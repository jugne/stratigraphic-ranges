package sr.evolution.operators;

import beast.base.evolution.tree.Node;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;
import beast.base.util.Randomizer;

import java.util.ArrayList;


public class MixedSpeciationSwap extends SRTreeOperator {

    @Override
    public void initAndValidate() {
    }

    /**
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
    	
    	SRTree tree = treeInput.get();

        // choose a random node avoiding leaves and nodes that are direct ancestors or belong to a range
        int nodeCount = tree.getNodeCount();
        SRNode node = null;

        for (int i=0; i<5; i++) {
            node = (SRNode) tree.getNode(Randomizer.nextInt(nodeCount));
            if (!node.isLeaf() && !node.isFake() && !tree.belongToSameSRange(node.getNr(), node.getLeft().getNr())) {
                break;
            }
            node = null;
        }


        if (node == null) {
            ArrayList<Integer> allowableNodeIndices = new ArrayList<>();

            for (int index=0; index<nodeCount; index++) {
                Node candidateNode = tree.getNode(index);
                //the node is not a leaf or sampled ancestor, the node is not fake, none of its children
                // belongs to the same srange as node
                if (!candidateNode.isLeaf() && !candidateNode.isFake() &&
                        !tree.belongToSameSRange(index, candidateNode.getLeft().getNr()))
                    allowableNodeIndices.add(index);
            }

            int allowableNodeCount = allowableNodeIndices.size();

            if (allowableNodeCount == 0) {
                return Double.NEGATIVE_INFINITY;
            }

            node = (SRNode) tree.getNode(allowableNodeIndices.get(Randomizer.nextInt(allowableNodeCount)));
        }
        
        node.setBudding(!node.isBudding());

        return 0.0;
    }
}
