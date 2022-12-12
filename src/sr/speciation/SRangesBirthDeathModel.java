package sr.speciation;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.core.Citation;
import beast.base.core.Description;


import sa.evolution.speciation.SABirthDeathModel;
import sr.evolution.tree.SRTree;
import sr.evolution.sranges.StratigraphicRange;

/**
 * @author Alexandra Gavryushkina
 */


@Description("A variant of the fossilized birth-death model under budding (asymmetric) speciation with stratigraphic ranges")
@Citation("Stadler T, Gavryushkina A, Warnock RCM, Drummond AJ, Heath TA (2017) \n" +
        "The fossilized birth-death model under different types of speciation")
@Citation("Gavryushkina A, Warnock RCM, Drummond AJ, Heath TA, Stadler T (2017) \n" +
        "Bayesian total-evidence dating under the fossilized birth-death model with stratigraphic ranges.")
public class SRangesBirthDeathModel extends SABirthDeathModel {

//    final public Input<SRTree> treeInput = new Input<>("srTree",
//            "beast.tree on which this operation is performed",
//            Input.Validate.REQUIRED);

    private double q_tilde(double t, double c1, double c2) {
        return Math.sqrt(Math.exp(t*(lambda + mu + psi))*q(t,c1,c2));
    }

    private double log_q_tilde(double t, double c1, double c2) {
        return 0.5*(t*(lambda + mu + psi) + log_q(t,c1,c2));
    }

    private double log_lambda_times_int_limits_p(double tOld, double tYoung, double c1, double c2) {
        return Math.log((lambda+mu+psi-c1)*(tOld - tYoung) + 2*Math.log(Math.exp(-c1*tYoung)*(1-c2)+(1+c2)) -
                2*Math.log(Math.exp(-c1*tOld)*(1-c2)+(1+c2))) - Math.log(2);
    }

    private Node findAncestralRangeLastNode(Node node) {
        Node parent = node.getParent();
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

        //double x0 = tree.getRoot().getHeight() + origToRootDistance;
        double x0 = origin;

//        if (taxonInput.get() != null) { //TODO rewrite this part for SA sRanges
//
//            if (taxonAge > origin) {
//                return Double.NEGATIVE_INFINITY;
//            }
//            double logPost = 0.0;
//
//            if (conditionOnSamplingInput.get()) {
//                logPost -= Math.log(oneMinusP0(x0, c1, c2));
//            }
//
//            if (conditionOnRhoSamplingInput.get()) {
//                logPost -= Math.log(oneMinusP0Hat(x0, c1, c2));
//            }
//
//            if (SATaxonInput.get().getValue() == 0) {
//                logPost += Math.log(1 - oneMinusP0(taxonAge, c1, c2));
//            } else {
//                logPost += Math.log(oneMinusP0(taxonAge, c1, c2));
//            }
//
//            return logPost;
//        }

        double x1=tree.getRoot().getHeight();

        if (x0 < x1 ) {
            return Double.NEGATIVE_INFINITY;
        }

        if (!conditionOnRootInput.get()){
            logP = -log_q(x0, c1, c2);
        } else {
            if (tree.getRoot().isFake()){   //when conditioning on the root we assume the process
                //starts at the time of the first branching event and
                //that means that the root can not be a sampled ancestor
                return Double.NEGATIVE_INFINITY;
            } else {
                logP = -log_q(x1, c1, c2);
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
                    if (node.getHeight() > 0.000000000005 || rho == 0.) {
                        Node fossilParent = node.getParent();
                        if (((SRTree)tree).belongToSameSRange(i, fossilParent.getNr())) {
                            logP += Math.log(psi) + log_q_tilde(node.getHeight(), c1, c2) + log_p0s(node.getHeight(), c1, c2);
                        } else {
                            logP += Math.log(psi) + log_q(node.getHeight(), c1, c2) + log_p0s(node.getHeight(), c1, c2);
                        }
                    } else {
                        logP += Math.log(4*rho);
                    }
                }
            } else {
                if (node.isFake()) {
                    logP += Math.log(psi);
                    Node parent = node.getParent();
                    Node child = node.getNonDirectAncestorChild();
                    if (parent != null && ((SRTree)tree).belongToSameSRange(parent.getNr(),i)) {
                        logP += log_q_tilde(node.getHeight(), c1, c2) - log_q(node.getHeight(), c1, c2);
                    }
                    if (child != null && ((SRTree)tree).belongToSameSRange(i,child.getNr())) {
                        logP += log_q(node.getHeight(), c1, c2)-  log_q_tilde(node.getHeight(), c1, c2);
                    }
                } else {
                    logP += Math.log(lambda) - log_q(node.getHeight(), c1, c2);
                }
            }
        }

        for (StratigraphicRange range:((SRTree)tree).getSRanges()) {
            Node first =  tree.getNode(range.getNodeNrs().get(0));
            if (!range.isSingleFossilRange()) {
                double tFirst =first.getHeight();
                double tLast = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size()-1)).getHeight();
                logP += psi*(tFirst - tLast);
            }
            Node ancestralLast = findAncestralRangeLastNode(first);
            if (ancestralLast != null) {
                double tOld = ancestralLast.getHeight();
                double tYoung = first.getHeight();
                logP += log_lambda_times_int_limits_p(tOld, tYoung, c1, c2);
            }
        }

        return logP;
    }

}

