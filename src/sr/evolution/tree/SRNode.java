package sr.evolution.tree;

import beast.base.core.BEASTInterface;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;

import java.util.TreeMap;

/**
 * Created by gavryusa on 04/05/17.
 * @author Alexandra Gavryushkina
 * @author Ugne Stolz
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
        Log.warning("Do not sort ordered trees. Tree not sorted");
        return 0;
    }

    /**
     * @return beast.tree in Newick format, with length and meta data
     *         information. Unlike toNewick(), here Nodes are numbered, instead of
     *         using the node labels.
     *         If there are internal nodes with non-null IDs then their numbers are also printed.
     *         Also, all internal nodes are labelled if printInternalNodeNumbers
     *         is set true. This is useful for example when storing a State to file
     *         so that it can be restored.
     */

    public String toShortNewickForLog(final boolean printInternalNodeNumbers) {
        final StringBuilder buf = new StringBuilder();

        if (!isLeaf()) {
            buf.append("(");
            boolean isFirst = true;
            for (Node child : getChildren()) {
                if (isFirst)
                    isFirst = false;
                else
                    buf.append(",");
                buf.append(((SRNode) child).toShortNewickForLog(printInternalNodeNumbers));
            }
            buf.append(")");
        }

        if (isLeaf() || getID() != null || printInternalNodeNumbers) {
            buf.append(getNr()+1);
        }

        buf.append(getNewickMetaData());
        buf.append(":").append(getNewickLengthMetaData()).append(getLength());
        return buf.toString();
    }
}
