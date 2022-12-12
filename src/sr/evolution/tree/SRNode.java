package sr.evolution.tree;

import beast.base.evolution.tree.Node;

import java.util.TreeMap;

/**
 * Created by gavryusa on 04/05/17.
 */
public class SRNode extends Node {

    /**
     * @return (deep) copy of node
     */
    public SRNode copy() {
        final SRNode node = new SRNode();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.metaData = new TreeMap<>(metaData);
        node.parent = null;
        node.setID(getID());

        for (final Node child : getChildren()) {
            node.addChild(child.copy());
        }
        return node;
    } // copy

    @Override
    public int sort()  {
        throw new RuntimeException("Do not sort ordered trees. Calculation stopped.");
    }
}
