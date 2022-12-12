package sr.evolution.operators;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import sr.evolution.tree.SRTree;

/**
 * copy from TreeOperator.
 */
abstract public class SRTreeOperator extends Operator {

    final public Input<SRTree> treeInput = new Input<>("tree",
            "beast.tree on which this operation is performed",
            Input.Validate.REQUIRED);
    final public Input<Boolean> markCladesInput = new Input<>("markClades",
            "Mark all ancestors of nodes changed by the operator as changed," +
            " up to the MRCA of all nodes changed by the operator.",
            false);

    /**
     * @param parent the parent
     * @param child  the child that you want the sister of
     * @return the other child of the given parent.
     */
    protected Node getOtherChild(final Node parent, final Node child) {
        if (parent.getLeft().getNr() == child.getNr()) {
            return parent.getRight();
        } else {
            return parent.getLeft();
        }
    }

    /**
     * replace child with another node
     *
     * @param parent parent node
     * @param child child node
     * @param replacement node that will replace the child
     */
    public void replace(final Node parent, final Node child, final Node replacement) {
        parent.removeChild(child);
        parent.addChild(replacement);
        parent.makeDirty(Tree.IS_FILTHY);
        replacement.makeDirty(Tree.IS_FILTHY);
    }

}
