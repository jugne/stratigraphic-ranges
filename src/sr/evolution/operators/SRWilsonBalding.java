package sr.evolution.operators;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;

import java.util.ArrayList;

/**
 * Implements the Wilson-Balding proposal for the sRange tree.
 * If the tree has a flag of being Infection Interval tree then additional moves are allowed which attach infection interval (sRange)
 * as a right branch at the bifurcation time equal to the start of the interval (firstOccurrence time)
 * @author Alexandra Gavryushkina
 * @author Ugne Stolz
 */
public class SRWilsonBalding extends SRTreeOperator {

    @Override
    public void initAndValidate() {

    }

    /**
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        SRTree tree = (SRTree) InputUtil.get(treeInput, this);
        Boolean proposeIITrees = tree.isIITree();

        double oldMinAge, newMinAge, newRange, oldRange, newAge, fHastingsRatio, dimensionCoefficient, orientationCoefficient;
        int newDimension, oldDimension;

        int nodeCount = tree.getNodeCount();

        ArrayList<Integer> allowableNodeIndices = new ArrayList<>();
        ArrayList<Integer> sRangeInternalNodeNrs = tree.getSRangesInternalNodeNrs();

        // choose a random node avoiding root and leaves that are direct ancestors
        for (int index=0; index<nodeCount; index++) {
            Node node = tree.getNode(index);
            //the node is not the root, it is not a sampled ancestor on a zero branch, it is not an internal node of a
            // stratigraphic range and is not the child of the root-origin node in case of IITrees
            if (!node.isRoot() && !node.isDirectAncestor() && !sRangeInternalNodeNrs.contains(node.getNr()) &&
                    (!proposeIITrees || !node.getParent().isRoot())) {
                allowableNodeIndices.add(index);
            }
        }

        Node i;

        int allowableNodeCount = allowableNodeIndices.size();

        if (allowableNodeCount == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        i=tree.getNode(allowableNodeIndices.get(Randomizer.nextInt(allowableNodeCount)));

        Node iP = i.getParent();
        Node CiP;
        if (iP.getLeft().getNr() == i.getNr()) {
            CiP = iP.getRight();
        } else {
            CiP = iP.getLeft();
        }

        // make sure that there is at least one candidate edge to attach node iP to
        if ((iP.getParent() == null || (proposeIITrees && iP.getParent().getParent()==null))  && CiP.getHeight() <= i.getHeight()) {
            return Double.NEGATIVE_INFINITY;
        }

        // choose another random node to insert i above or to attach i to this node if it is a leaf
        Node j;
        Node jP;

        final int leafNodeCount = tree.getLeafNodeCount();

        if (leafNodeCount != tree.getExternalNodes().size()) {
            System.out.println("node counts are incorrect. NodeCount = " + nodeCount + " leafNodeCount = " + leafNodeCount + " external node count = " + tree.getExternalNodes().size());
        }

        // make sure that the target branch <jP, j> or target leaf j is above the subtree being moved

        int nodeNumber;
        double newParentHeight;
        boolean attachingToLeaf;
        boolean adjacentEdge;
        boolean originRootBranch;
        //boolean adjacentLeaf;
        do {
            adjacentEdge = false;
            //adjacentLeaf = false;
            nodeNumber = Randomizer.nextInt(nodeCount + leafNodeCount);
            if (nodeNumber < nodeCount) {
                j = tree.getNode(nodeNumber);
                jP = j.getParent();
                if (jP != null)
                    newParentHeight = jP.getHeight();
                else newParentHeight = Double.POSITIVE_INFINITY;
                if (!CiP.isDirectAncestor())
                    adjacentEdge = (CiP.getNr() == j.getNr() || iP.getNr() == j.getNr());
                originRootBranch = proposeIITrees && (jP == null || (jP.isRoot() && jP.getHeight()==j.getHeight()));
                attachingToLeaf = false;
            } else {
                j = tree.getExternalNodes().get(nodeNumber - nodeCount);
                jP = j.getParent();
                newParentHeight = j.getHeight();
                attachingToLeaf = true;
                originRootBranch = false;
                //adjacentLeaf = (iP.getNr() == j.getNr());
            }
        } while (j.isDirectAncestor() || (newParentHeight <= i.getHeight()) || (i.getNr() == j.getNr()) ||
                adjacentEdge || originRootBranch /*|| adjacentLeaf */);


        if (attachingToLeaf && iP.getNr() == j.getNr()) {
            System.out.println("Proposal failed because j = iP");
            return Double.NEGATIVE_INFINITY;
        }

        if (jP != null && jP.getNr() == i.getNr()) {
            System.out.println("Proposal failed because jP = i. Heights of i = " + i.getHeight() + " Height of jP = " + jP.getHeight());
            return Double.NEGATIVE_INFINITY;
        }



        //oldDimension = nodeCount - tree.getDirectAncestorNodeCount() - 1;
        oldDimension = allowableNodeCount;
        orientationCoefficient = 1.0;
        StratigraphicRange pruningRange = null;
        StratigraphicRange attachingRange = null;

        //classify the type of move being performed before changing the tree structure
        boolean pruningFromSA = CiP.isDirectAncestor();
//        boolean pruningFromSRange = !CiP.isDirectAncestor() && sRangeInternalNodeNrs.contains(iP.getNr());
        boolean randomAttach = false;
        boolean pruningFromSRange = false;
        if (!CiP.isDirectAncestor()) {
            pruningRange = tree.getSharedRange(iP.getNr(), CiP.getNr());
            pruningFromSRange = pruningRange != null;
        }
        boolean randomPrune = !pruningFromSA && !pruningFromSRange && !(tree.isIITree() && i.isFake() && i.getHeight() == iP.getHeight());
        boolean attachingToSRange = false;
        if (!attachingToLeaf && jP != null){
            attachingRange = tree.getSharedRange(jP.getNr(),j.getNr());
            attachingToSRange = attachingRange!=null;
        }

        boolean observedTransmissionAttach = false;

        newRange=1;
        //Hastings numerator calculation + newAge of iP
        if (attachingToLeaf) {
            newAge = j.getHeight();
        } else {
            if (jP != null) {
                if (!tree.isIITree() || !i.isFake() || j.getHeight() >= i.getHeight() || Randomizer.nextBoolean()) {
                    newMinAge = Math.max(i.getHeight(), j.getHeight());
                    newRange = jP.getHeight() - newMinAge;
                    newAge = newMinAge + (Randomizer.nextDouble() * newRange);
                } else {
                    newAge = i.getHeight();
                    observedTransmissionAttach = true;
                }

            } else {
                if (!tree.isIITree() || !i.isFake() || j.getHeight() >= i.getHeight() ||  Randomizer.nextBoolean()) {
                    double randomNumberFromExponential;
                    randomNumberFromExponential = Randomizer.nextExponential(1);
                    newRange = Math.exp(randomNumberFromExponential);
                    newAge = j.getHeight() + randomNumberFromExponential;
                } else {
                    newAge = i.getHeight();
                    observedTransmissionAttach = true;
                }
            }
            // if there were two choices whether to attach iP at a random height between the height of j and
            // the height of jP (or infinity) or at the height of i (when creating observed transmission),
            // then adjust the hastings by the probability of each choice of 1/2.
            // As instead of the density, we have its reciprocal, we multiply by 2.
            if (tree.isIITree() && i.isFake() && j.getHeight() < i.getHeight()) {
                newRange=newRange*2;
            }
        }

        Node PiP = iP.getParent();

        oldRange=1;

        //Hastings denominator calculation
        if (!CiP.isDirectAncestor()) {
            if (!tree.isIITree() || !i.isFake() || i.getHeight() != iP.getHeight()) {
                oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
                if (PiP != null) {
                    oldRange = PiP.getHeight() - oldMinAge;
                } else {
                    oldRange = Math.exp(iP.getHeight() - oldMinAge);
                    //oldRange = x0 - oldMinAge;
                }
            }
            // if node i is the beginning of an infection interval, then the reverse move
            // has two choices. In both cases (iP is a normal bifurcation node or
            // it is a bifurcation node that is an observed transmission) we need to adjust the density of the reverse move
            // by the probability of choosing one or the other which is 1/2.
            // As instead of the density, we have its reciprocal, we multiply by 2.
            if (tree.isIITree() && i.isFake() && CiP.getHeight() < i.getHeight()) {
                oldRange = oldRange * 2;
            }
        }

        //update
        if (iP.getNr() != j.getNr() && CiP.getNr() != j.getNr()) { // special case 1: iP = j when pruning from sampled ancestor iP and attaching to the branch above (PiP, iP)
                                                                    // special case 2: CiP = j when pruning from a branch (PiP, CiP) and attaching to a leaf CiP
                                                                    // In both cases the internal tree structure does not change, only the height of iP
            iP.removeChild(CiP); //remove <iP, CiP>

            if (PiP != null) {
                boolean left = PiP.getLeft().getNr() == iP.getNr();
                Node anotherChild;
                if (left) {
                    anotherChild = PiP.getRight();
                } else {
                    anotherChild = PiP.getLeft();
                }
                PiP.removeChild(iP);   // remove <PiP,iP>
                CiP.setParent(PiP);
                if (left) {
                    PiP.setLeft(CiP);
                    PiP.setRight(anotherChild);
                } else {
                    PiP.setLeft(anotherChild);
                    PiP.setRight(CiP);
                }  // add <PiP, CiP> at the same orientation as <PiP, iP>
                PiP.makeDirty(Tree.IS_FILTHY);
                CiP.makeDirty(Tree.IS_FILTHY);
            } else {
                CiP.setParent(null); // completely remove <iP, CiP>
                tree.setRootOnly(CiP);
            }

            boolean jLeft= (jP != null) && jP.getLeft().getNr() == j.getNr();

            if (jP != null) {
                Node CjP;
                if (jLeft) {
                    CjP = jP.getRight();
                } else {
                    CjP = jP.getLeft();
                }
                jP.removeChild(j);  // remove <jP, j>
                iP.setParent(jP);
                if (jLeft) {
                    jP.setLeft(iP);
                    jP.setRight(CjP);
                } else {
                    jP.setLeft(CjP);
                    jP.setRight(iP);
                } // add <jP, iP> at the same orientation as <jP, j>

                jP.makeDirty(Tree.IS_FILTHY);
            } else {
                iP.setParent(null); // completely remove <PiP, iP>
                tree.setRootOnly(iP);
            }
            j.setParent(iP);
            if (!attachingToSRange && !attachingToLeaf && !observedTransmissionAttach)
                randomAttach = true;
            if (attachingToSRange || observedTransmissionAttach ||
                    (!attachingToLeaf && Randomizer.nextBoolean())) {
                iP.setLeft(j);
                iP.setRight(i);
            } else {
                iP.setLeft(i);
                iP.setRight(j);
            }
            iP.makeDirty(Tree.IS_FILTHY);
            j.makeDirty(Tree.IS_FILTHY);
        } else {
            if (iP.getNr() == j.getNr()) {
                if (!attachingToSRange && !observedTransmissionAttach)
                    randomAttach = true;
                if (attachingToSRange || observedTransmissionAttach || Randomizer.nextBoolean()) { //in special case 1: when attaching to the range
                                                                     //make i right
                                                                     //otherwise choose randomly
                    iP.setLeft(CiP);
                    iP.setRight(i);
                } else {
                    iP.setLeft(i);
                    iP.setRight(CiP);
                }
            }
            if (CiP.getNr() == j.getNr()) {//in special case 2: always make i left
                    iP.setLeft(i);
                    iP.setRight(CiP);
            }
        }
        iP.setHeight(newAge);

        // remove or add nodes to the ranges and calculate the orientation coefficient
        if (pruningFromSRange) {
            pruningRange.removeNodeNr(tree, iP.getNr());
        }
        if (attachingToSRange) {
            attachingRange.addNodeNrAfter(tree, jP.getNr(), iP.getNr());
        }
        if (randomAttach){
            orientationCoefficient *= 2.0;
        }
        if (randomPrune){
            orientationCoefficient *= 0.5;
        }

        newDimension = 0;
        sRangeInternalNodeNrs = tree.getSRangesInternalNodeNrs();

        for (int index=0; index<nodeCount; index++) {
            Node node = tree.getNode(index);
            //the node is not the root, it is not a sampled ancestor on a zero branch,
            // it is not an internal node of a stratigraphic range
            // it is not a child of the root in case of II trees
            if (!node.isRoot() && !node.isDirectAncestor() && !sRangeInternalNodeNrs.contains(node.getNr()) &&
                    (!proposeIITrees || !node.getParent().isRoot()))
                newDimension++;
        }
        dimensionCoefficient = (double) oldDimension / newDimension;

        fHastingsRatio = Math.abs(orientationCoefficient * dimensionCoefficient * newRange / oldRange);

        return Math.log(fHastingsRatio);

    }

}
