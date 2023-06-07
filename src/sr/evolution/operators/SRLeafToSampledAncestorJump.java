package sr.evolution.operators;



import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;

import java.util.*;

/**
 * @author Alexandra Gavryushkina
 * @author Ugne Stolz
 */

@Description("Implements a narrow move between trees of different dimensions (number of nodes in trees)." +
        "It takes a random sampled node which is either a leaf with the younger sibling" +
        "or a sampled internal node. In the first case, the leaf becomes a sampled internal node by replacing its " +
        "parent (the height of the leaf remains unchanged). In the second case, the sampled internal node becomes " +
        "a leaf by inserting a new parent node at a height which is uniformly chosen on the interval " +
        " between the sampled node height and its old parent height.")
public class SRLeafToSampledAncestorJump extends SRTreeOperator {

    public Input<IntegerParameter> categoriesInput = new Input<IntegerParameter>("rateCategories", "rate category per branch");

    public Input<RealParameter> rInput =
            new Input<RealParameter>("removalProbability", "The probability of an individual to be removed from the process immediately after the sampling");

    @Override
    public void initAndValidate() {
    }

    @Override
    public double proposal() {

        double newHeight, newRange, oldRange, orientationCoefficient;
        StratigraphicRange sameRange = null;
        orientationCoefficient = 1;
        int categoryCount = 1;
        if (categoriesInput.get() != null) {
            categoryCount = categoriesInput.get().getUpper() - categoriesInput.get().getLower() +1;
        }

        SRTree tree = (SRTree) treeInput.get();

		Integer[] fitLeafNodeNrs = getFitToMoveNodeNrs(tree);
		if (fitLeafNodeNrs.length == 0)
			return Double.NEGATIVE_INFINITY;

		SRNode leaf = (SRNode) tree.getNode(fitLeafNodeNrs[Randomizer.nextInt(fitLeafNodeNrs.length)]);
        SRNode parent = (SRNode) leaf.getParent();
        if (Math.abs(leaf.getHeight()-parent.getHeight())>0 && Math.abs(leaf.getHeight()-parent.getHeight())<0.00000005){
            System.out.println("");
        }

        if (leaf.isDirectAncestor()) {
            oldRange = 1;
            if (parent.isRoot()) {
                final double randomNumber = Randomizer.nextExponential(1);
                if (randomNumber<0.00000005)
                    System.out.println("");
                newHeight = parent.getHeight() + randomNumber;
                newRange = Math.exp(randomNumber);
            } else {
                newRange = parent.getParent().getHeight() - parent.getHeight();
                newHeight = parent.getHeight() + Randomizer.nextDouble() * newRange;
            }

            if (categoriesInput.get() != null) {
                int index = leaf.getNr();
                int newValue = Randomizer.nextInt(categoryCount) + categoriesInput.get().getLower(); // from 0 to n-1, n must > 0,
                categoriesInput.get().setValue(index, newValue);
            }

            SRNode otherChild = (SRNode) getOtherChild(parent,leaf);
            sameRange = tree.getSharedRange(otherChild.getNr(),leaf.getNr());
            if (sameRange==null && Randomizer.nextBoolean()){
                parent.removeAllChildren(true);
                parent.addChild(leaf);
                parent.addChild(otherChild);
            } else {
                parent.removeAllChildren(true);
                parent.addChild(otherChild);
                parent.addChild(leaf);
                if (sameRange!=null)
                    sameRange.removeNodeNr(tree, parent.getNr());
            }
            parent.makeAllDirty(tree.IS_FILTHY);
            if (sameRange==null)
                orientationCoefficient*=2.;
        } else {
            newRange = 1;
            SRNode otherChild = (SRNode) getOtherChild(parent,leaf);
            //make sure that the branch where a new sampled node to appear is not above that sampled node
            if (otherChild.getHeight() >= leaf.getHeight() )
                return Double.NEGATIVE_INFINITY;

            if (parent.isRoot()) {
                oldRange = Math.exp(parent.getHeight() - leaf.getHeight());
            } else {
                int siblingNr =  getOtherChild(parent, leaf).getNr();
                int grandParentNr =   parent.getParent().getNr();
                if (tree.belongToSameSRange(siblingNr,grandParentNr))
                    return Double.NEGATIVE_INFINITY;

                oldRange = parent.getParent().getHeight() - leaf.getHeight();
            }
            newHeight = leaf.getHeight();
            if  (categoriesInput.get() != null) {
                int index = leaf.getNr();
                categoriesInput.get().setValue(index, -1);
            }

            sameRange = tree.getSharedRange(otherChild.getNr(),leaf.getNr());
            StratigraphicRange leafRange = tree.getRangeOfNode(leaf);
            if (parent.getChild(0).getNr()!=otherChild.getNr()){
                parent.removeAllChildren(true);
                parent.addChild(otherChild);
                parent.addChild(leaf);
            }
            if (sameRange == null) {
                orientationCoefficient *= 0.5;
            } else {
                System.out.println("There is a bug in the SRLeafToSampledAncestorJump operator. " +
                        "Please report this to the developers!");
                System.exit(1);
            }
            leafRange.removeNodeNr(tree, parent.getNr());
            parent.makeAllDirty(tree.IS_FILTHY);
        }
        parent.setHeight(newHeight);

        //make sure that either there are no direct ancestors or r<1
        if ((rInput.get() != null) && (tree.getDirectAncestorNodeCount() > 0 && rInput.get().getValue() == 1))  {
            return Double.NEGATIVE_INFINITY;
        }

        return Math.log(orientationCoefficient*newRange/oldRange);
    }

	private Integer[] getFitToMoveNodeNrs(SRTree tree) {
		Set<Integer> fitNodeNrs = new HashSet<Integer>();
        for (StratigraphicRange range : tree.getSRanges()){
            if (range.isSingleFossilRange())
                fitNodeNrs.add(range.getNodeNrs().get(0));
        }
		return fitNodeNrs.toArray(new Integer[0]);
	}
}
