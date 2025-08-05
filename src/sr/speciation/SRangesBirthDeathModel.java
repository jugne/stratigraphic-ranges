package sr.speciation;

import beast.base.evolution.tree.Node;
import beast.base.core.Citation;
import beast.base.core.Description;


import sa.evolution.speciation.SABirthDeathModel;
import sr.evolution.tree.SRTree;
import sr.evolution.sranges.StratigraphicRange;

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



    @Override
    public double calculateLogP()
    {
        SRTree tree = (SRTree) combinedTree.getTree();
        int nodeCount = tree.getNodeCount();
        updateParameters();
        if (lambdaExceedsMu && lambda <= mu) {
            return Double.NEGATIVE_INFINITY;
        }

        if (lambda < 0 || mu < 0 || psi < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        double x0 = origin;
        double x1 = combinedTree.getHeightOfNode(tree.getRoot().getNr());

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
                    if (combinedTree.getHeightOfNode(i) > 0.000000000005 || rho == 0.) {

                        if (((SRTree)tree).belongToSameSRange(i, fossilParent.getNr())) {
                            logP += Math.log(psi) - log_q_tilde(combinedTree.getHeightOfNode(i), c1, c2) + log_p0s(combinedTree.getHeightOfNode(i), c1, c2);
                        } else {
                            logP += Math.log(psi) - log_q(combinedTree.getHeightOfNode(i), c1, c2) + log_p0s(combinedTree.getHeightOfNode(i), c1, c2);
                        }
                    } else {
                        logP += Math.log(rho);
                    }
                }
            } else {
                if (node.isFake()) {
                    logP += Math.log(psi);
                    Node parent = node.getParent();
                    Node child = node.getNonDirectAncestorChild();
                    Node DAchild = node.getDirectAncestorChild();
                    if (parent != null && ((SRTree)tree).belongToSameSRange(parent.getNr(),DAchild.getNr())) {
                        logP += - log_q_tilde(combinedTree.getHeightOfNode(i), c1, c2) + log_q(combinedTree.getHeightOfNode(i), c1, c2);
                    }
                    if (child != null && ((SRTree)tree).belongToSameSRange(i,child.getNr())) {
                        logP += - log_q(combinedTree.getHeightOfNode(i), c1, c2) +  log_q_tilde(combinedTree.getHeightOfNode(i), c1, c2);
                    }
                } else {
                    logP += Math.log(lambda) + log_q(combinedTree.getHeightOfNode(i), c1, c2);
                }
            }
        }

        // integrate over fossils in the range. This seems to suggest that we take out the psi in the previous equations
        for (StratigraphicRange range:((SRTree)tree).getSRanges()) {
            Node first =  tree.getNode(range.getNodeNrs().get(0));
            if (!range.isSingleFossilRange()) {
                double tFirst = combinedTree.getHeightOfNode(first.getNr());
                double tLast = combinedTree.getHeightOfNode(range.getNodeNrs().get(range.getNodeNrs().size()-1));
                logP += psi*(tFirst - tLast);
            }
            Node ancestralLast = findAncestralRangeLastNode(first);
            if (ancestralLast != null) {
                double tOld = combinedTree.getHeightOfNode(ancestralLast.getNr());
                double tYoung = combinedTree.getHeightOfNode(first.getNr());
                logP += Math.log(1-q(tYoung, c1, c2)/q_tilde(tYoung, c1, c2)*q_tilde(tOld, c1, c2)/q(tOld, c1, c2));
            }
        }
        return logP;
    }

}

