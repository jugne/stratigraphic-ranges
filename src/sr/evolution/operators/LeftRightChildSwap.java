package sr.evolution.operators;

import beast.base.evolution.tree.Node;
import beast.base.inference.util.InputUtil;
import sr.evolution.tree.SRTree;
import beast.base.util.Randomizer;

import java.util.ArrayList;

/**
 * @author Alexandra Gavryushkina
 */
public class LeftRightChildSwap extends SRTreeOperator {

    @Override
    public void initAndValidate() {
    }

    /**
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        SRTree tree = (SRTree) InputUtil.get(treeInput, this);


        int nodeCount = tree.getNodeCount();

        Node node = null;

        // choose a random node that is a bifurcation node and neither of
        // the children are range branches
        for (int i=0; i<5; i++) {
            node = tree.getNode(Randomizer.nextInt(nodeCount));
            if (!node.isLeaf() && !node.isFake() && !tree.belongToSameSRange(node.getNr(), node.getLeft().getNr()) &&
            !(node.getHeight()==node.getRight().getHeight())) {
                break;
            }
            node = null;
        }


        if (node == null) {
            ArrayList<Integer> allowableNodeIndices = new ArrayList<>();

            for (int index=0; index<nodeCount; index++) {
                Node candidateNode = tree.getNode(index);
                //the node is not a leaf or sampled ancestor, the node is not fake, none of its children
                // belongs to the same sRange as the node, node is not an observed transmission
                if (!candidateNode.isLeaf() && !candidateNode.isFake() &&
                        !tree.belongToSameSRange(index, candidateNode.getLeft().getNr()) &&
                        !(candidateNode.getHeight()==candidateNode.getRight().getHeight()))
                    allowableNodeIndices.add(index);
            }

            int allowableNodeCount = allowableNodeIndices.size();

            if (allowableNodeCount == 0) {
                return Double.NEGATIVE_INFINITY;
            }

            node=tree.getNode(allowableNodeIndices.get(Randomizer.nextInt(allowableNodeCount)));
        }


        Node left = node.getLeft();
        Node right = node.getRight();

        node.setLeft(right);
        node.setRight(left);

        return 0.0;
    }
}
