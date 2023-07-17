package sr.speciation;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.core.Citation;
import beast.base.core.Description;


import beast.base.inference.parameter.RealParameter;
import sa.evolution.speciation.SABDParameterization;
import sa.evolution.speciation.SABirthDeathModel;
import sr.evolution.tree.SRTree;
import sr.evolution.sranges.StratigraphicRange;

import static sr.util.Tools.getEndSubrangeLength;
import static sr.util.Tools.getStartSubrangeLength;

/**
 * @author Alexandra Gavryushkina
 * @author Ugne Stolz
 */


@Description("A variant of the fossilized birth-death model under budding (asymmetric) speciation with stratigraphic ranges")
@Citation("Stadler T, Gavryushkina A, Warnock RCM, Drummond AJ, Heath TA (2017) \n" +
        "The fossilized birth-death model under different types of speciation")
@Citation("Gavryushkina A, Warnock RCM, Drummond AJ, Heath TA, Stadler T (2017) \n" +
        "Bayesian total-evidence dating under the fossilized birth-death model with stratigraphic ranges.")
public class SRangesBirthDeathModel extends SABirthDeathModel {

    public Input<SABDParameterization> parameterizationInput = new Input<>("parameterization", "The parameterization to use.");
    final public Input<Boolean> integrateOverRangesInput = new Input<Boolean>("integrateOverRanges",
            "If true, there should be no samples in between range bounds and these samples are integrated over." +
                    "Otherwise, no integration happens. You should not provide in between samples, " +
                    "when the data they provide is minimal and sampling rate inconsistent with the rest of the tree.", true);

    //     Next two inputs only relevant for transmission applications when we want to include intermediate samples in
//     the range and also set a period before or after sampling where we integrate over the sampling
    final public Input<Boolean> startPeriodInput = new Input<Boolean>("startPeriod",
            "If true, integrate over sampling between empty start sample and first non-empty sample of the range.",
            false);

    final public Input<Boolean> endPeriodInput = new Input<Boolean>("endPeriod",
            "If true, integrate over sampling between last non-empty sample of the range and empty end sample at the end.",
            false);


    @Override
    public double q(double t, double c1, double c2) {
        double v = Math.exp(-c1 * t);
        return 4 * v / Math.pow(v*(1-c2) + (1+c2), 2.0);
    }

    @Override
    public double log_q(double t, double c1, double c2) {
        return Math.log(q(t,c1,c2));
    }

    private double q_tilde(double t, double c1, double c2) {
        return Math.sqrt(Math.exp(-t*(lambda + mu + psi))*q(t,c1,c2));
    }

    private double log_q_tilde(double t, double c1, double c2) {
        return 0.5*(-t*(lambda + mu + psi) + log_q(t,c1,c2));
    }

    private Node findAncestralRangeLastNode(Node node) {
        Node parent = node.getParent();
        if (node.isDirectAncestor()){
            parent = parent.getParent();
            node = node.getParent();
        }
        if (parent == null) {
            return parent;
        } else {
            if (parent.isFake()) {
                return parent;
            } else if (parent.getChild(0) == node) {
                return findAncestralRangeLastNode(parent);
            } else {
                return null;
            }
        }
    }

    protected void updateParameters() {
        if (parameterizationInput.get()==null){
            super.updateParameters();
        } else {
            lambda = parameterizationInput.get().lambda();
            mu = parameterizationInput.get().mu();
            psi = parameterizationInput.get().psi();
            if (!conditionOnRootInput.get()){
                origin = parameterizationInput.get().origin();
            }  else {
                origin = Double.POSITIVE_INFINITY;
            }

            r = removalProbability.get().getValue();
            if (rhoProbability.get() != null ) {
                rho = rhoProbability.get().getValue();
            } else {
                rho = 0.;
            }
            c1 = Math.sqrt((lambda - mu - psi) * (lambda - mu - psi) + 4 * lambda * psi);
            c2 = -(lambda - mu - 2*lambda*rho - psi) / c1;
        }
    }


    boolean integrateOverRanges;
    boolean useStartPeriod;
    boolean useEndPeriod;
    @Override
    public void initAndValidate() {
        if (!integrateOverRangesInput.get() && (startPeriodInput.get() || endPeriodInput.get())) {
            throw new IllegalArgumentException("You can only use startPeriod and endPeriod when " +
                    "integrateOverRanges is true.");
        }
        super.initAndValidate();
        integrateOverRanges = integrateOverRangesInput.get();
        useStartPeriod = startPeriodInput.get();
        useEndPeriod = endPeriodInput.get();


    }

    @Override
    public double calculateLogP()
    {
        SRTree tree = (SRTree) treeInput.get();
        int nodeCount = tree.getNodeCount();
        updateParameters();
        if (lambdaExceedsMu && lambda <= mu) {
            return Double.NEGATIVE_INFINITY;
        }

        if (lambda < 0 || mu < 0 || psi < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        double x0 = origin;
        double x1=tree.getRoot().getHeight();

        if (x0 < x1 ) {
            return Double.NEGATIVE_INFINITY;
        }

        if (!conditionOnRootInput.get()){
            logP = log_q(x0, c1, c2);
        } else {
            if (tree.getRoot().isFake()){   //when conditioning on the root we assume the process
                //starts at the time of the first branching event and
                //that means that the root can not be a sampled ancestor
                return Double.NEGATIVE_INFINITY;
            } else {
                logP = log_q(x1, c1, c2);
            }
        }

        if (conditionOnSamplingInput.get()) {
            logP -= log_oneMinusP0(x0, c1, c2);
        }

        if (conditionOnRhoSamplingInput.get()) {
            if (conditionOnRootInput.get()) {
                logP -= Math.log(lambda) + log_oneMinusP0Hat(x1, c1, c2)+ log_oneMinusP0Hat(x1, c1, c2);
            }  else {
                logP -= log_oneMinusP0Hat(x0, c1, c2);
            }
        }

        for (int i = 0; i < nodeCount; i++) {
            Node node = tree.getNode(i);
            if (node.isLeaf()) {
                if  (!node.isDirectAncestor())  {
                    Node fossilParent = node.getParent();
                    if (node.getHeight() > 0.000000000005 || rho == 0.) {

                        if (((SRTree)tree).belongToSameSRange(i, fossilParent.getNr())) {
                            logP += Math.log(psi) - log_q_tilde(node.getHeight(), c1, c2) + log_p0s(node.getHeight(), c1, c2);
                        } else {
                            logP += Math.log(psi) - log_q(node.getHeight(), c1, c2) + log_p0s(node.getHeight(), c1, c2);
                        }
                    } else {
                        logP += Math.log(rho);
                    }
                }
            } else {
                if (node.isFake()) {
                    if (r == 1) {
                        System.out.println("r = 1 but there are sampled ancestors in the tree");
                        System.exit(0);
                    }
                    logP += Math.log(psi) + Math.log(1 - r);
                    Node parent = node.getParent();
                    Node child = node.getNonDirectAncestorChild();
                    Node DAchild = node.getDirectAncestorChild();
                    if (parent != null && tree.belongToSameSRange(parent.getNr(),DAchild.getNr())) {
                        logP += - log_q_tilde(node.getHeight(), c1, c2) + log_q(node.getHeight(), c1, c2);
                    }
                    if (child != null && tree.belongToSameSRange(i,child.getNr())) {
                        logP += - log_q(node.getHeight(), c1, c2) +  log_q_tilde(node.getHeight(), c1, c2);
                    }
                } else {
                    logP += Math.log(lambda) + log_q(node.getHeight(), c1, c2);
                }
            }
        }

        // integrate over fossils in the range. This seems to suggest that we take out the psi in the previous equations
        for (StratigraphicRange range:(tree).getSRanges()) {
            Node first =  tree.getNode(range.getNodeNrs().get(0));
            if (!range.isSingleFossilRange()) {
                double tFirst = first.getHeight();
                int rangeSize =  range.getNodeNrs().size();
                if (!integrateOverRanges && rangeSize > 2){
                    if (useStartPeriod)
                        logP += (psi*(1-r))*(getStartSubrangeLength(range, tree));
                    if (useEndPeriod)
                        logP += (psi*(1-r))*(getEndSubrangeLength(range, tree));
                } else {
                    double tLast = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size()-1)).getHeight();
                    logP += (psi*(1-r))*(tFirst - tLast);
                }

//                if (useStartPeriod){
//                    logP += psi*(getStartSubrangeLength(range, tree));
//                } else if (useEndPeriod) {
//                    logP += psi*(getEndSubrangeLength(range, tree));
//                } else {
//                    double tLast = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size()-1)).getHeight();
//                    logP += psi*(tFirst - tLast);
//                }
//                double tLast = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size()-1)).getHeight();

            }
            Node ancestralLast = findAncestralRangeLastNode(first);
            if (ancestralLast != null) {
                double tOld = ancestralLast.getHeight();
                double tYoung = first.getHeight();
                logP += Math.log(1-q(tYoung, c1, c2)/q_tilde(tYoung, c1, c2)*q_tilde(tOld, c1, c2)/q(tOld, c1, c2));
            }
        }
        return logP;
    }

}

