package sr.evolution.sranges;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import sr.evolution.tree.SRTree;

import java.util.ArrayList;
import java.util.List;

/**
 *@author Alexandra Gavryushkina
 * @author Ugne Stolz
 */
public class StratigraphicRange extends BEASTObject {

    public final Input<Taxon> taxonFirstOccurrenceInput = new Input<Taxon>("firstOccurrence", "A BEAST taxon object that corresponds to the first " +
            "occurrence of the taxon");
    public final Input<Taxon> taxonLastOccurrenceInput = new Input<Taxon>("lastOccurrence", "A BEAST taxon object that corresponds to the last " +
            "occurrence of the taxon");


    String firstOccurrenceID=null;

    String lastOccurrenceID=null;

    boolean isSingleFossilRange=false;

    /**
     * The list of nodes that belong to the range.
     * The nodes in the list should go in the descending order with respect to the height.
     * The node at position 0 always corresponds to the first occurrence of the fossil.
     * The last node in the list always corresponds to the last occurrence of the fossil.
     * The intermediate nodes (if any) are always branching nodes.
     * For a single fossil range the first and the last occurrences coincide and there is only a single.
     * node in the list.
     */
    private List<Integer> nodes = new ArrayList<>();  //

    private List<Integer> directAncestorNodes = new ArrayList<>();  //

    @Override
    public void initAndValidate() {
        if (taxonFirstOccurrenceInput.get() != null || taxonLastOccurrenceInput.get() != null) {
            if (taxonFirstOccurrenceInput == null || taxonLastOccurrenceInput.get() == null) {
                throw new RuntimeException("Either firstOccurrence or lastOccurence is not specified. Either specify both or none");
            }
            firstOccurrenceID = taxonFirstOccurrenceInput.get().getID();
            lastOccurrenceID = taxonLastOccurrenceInput.get().getID();
            isSingleFossilRange = (firstOccurrenceID.equals(lastOccurrenceID));
        }
    }

    public boolean containsNodeNr(SRTree tree, int nodeNr) {
        int nr = nodeNr;
        if (tree.getNode(nodeNr).isFake())
            nr = tree.getNode(nodeNr).getDirectAncestorChild().getNr();
        return nodes.contains(nr);
    }

    public void addNodeNrAfter(SRTree tree, int nodeAfterNr, int nodeNr) {
        int afterNr = nodeAfterNr;
        int nr = nodeNr;
        if (tree.getNode(afterNr).isFake())
            afterNr = tree.getNode(afterNr).getDirectAncestorChild().getNr();
        if (tree.getNode(nodeNr).isFake())
            nr = tree.getNode(nodeNr).getDirectAncestorChild().getNr();
        int i= nodes.indexOf(afterNr)+1;
        nodes.add(i,nr);
    }

    public void removeNodeNr(SRTree tree, int nodeNr) {
        nodes.remove((Integer)nodeNr);
    }

    public void removeAllNodeNrs() {
        nodes.clear();
    }

    public List<Integer> getNodeNrs() {
        return nodes;
    }

    /**
     * sets the node that corresponds to the the first occurrence of the range in the tree at position 0
     * the single fossil is treated as the first occurrence
     * @param nodeNr
     */
    public void setFirstOccurrenceNodeNr(SRTree tree, int nodeNr) {
        int nr = nodeNr;
        if (tree.getNode(nodeNr).isFake())
            nr = tree.getNode(nodeNr).getDirectAncestorChild().getNr();
        if (nodes.isEmpty()) {
            nodes.add(nr);
        } else {
            nodes.set(0,nr);
        }
    }

    /**
     * adds the node that corresponds to the last occurrence of the range in the tree at the last position
     * in the array. A single fossil is treated as the first and the last occurrence.
     *
     * @param nodeNr
     */
    public void setLastOccurrenceNodeNr(SRTree tree, int nodeNr) {
        int nr = nodeNr;
        if (tree.getNode(nodeNr).isFake())
            nr = tree.getNode(nodeNr).getDirectAncestorChild().getNr();
        if (isSingleFossilRange()) {
            if (nodes.isEmpty()) {
                nodes.add(nr);
            } else {
                nodes.set(0,nr);
            }
            return;
        } else {
            if (nodes.isEmpty()) {
                nodes.add(null);
            }
        }
        nodes.add(nr);
    }

    public void makeSingleFossilRange() {
        if (firstOccurrenceID != null && lastOccurrenceID != null && !firstOccurrenceID.equals(lastOccurrenceID))
            throw new RuntimeException("First and last occurrences of taxon "+ getID() + "are not equal but an attempt" +
                    "to mark it as a single fossil range was made");
        if (firstOccurrenceID == null && lastOccurrenceID != null) firstOccurrenceID = lastOccurrenceID;
        if (lastOccurrenceID == null && firstOccurrenceID != null) lastOccurrenceID = firstOccurrenceID;
        isSingleFossilRange = true;
        if (nodes.get(0)==null){
            nodes.remove(0);
        }
    }

    public boolean isSingleFossilRange() {
        if (!isSingleFossilRange && taxonFirstOccurrenceInput.get() != null && taxonLastOccurrenceInput.get() != null &&
                taxonFirstOccurrenceInput.get().getID().equals(taxonLastOccurrenceInput.get().getID())) {
            isSingleFossilRange = true;
        }
        if (!isSingleFossilRange && firstOccurrenceID != null && lastOccurrenceID != null &&
                firstOccurrenceID.equals(lastOccurrenceID)) {
            isSingleFossilRange = true;
        }
        return isSingleFossilRange;
    }

    public List<Integer> getInternalNodeNrs(SRTree tree) {
        List<Integer> internalNodeNrs = new ArrayList<>();
        for (int i=1; i< nodes.size(); i++) {
            internalNodeNrs.add(nodes.get(i));
            if (tree.getNode(i).isDirectAncestor()){
                internalNodeNrs.add(tree.getNode(i).getParent().getNr());
            }
        }
        return internalNodeNrs;
    }

    public void setFirstOccurrenceID(String ID) {
        if (taxonFirstOccurrenceInput.get() != null && ! taxonFirstOccurrenceInput.get().getID().equals(ID)) {
            throw new RuntimeException(ID + " was attempted to be assigned as the name of the first occurrence taxon " +
                    "for a stratigraphic range for which the first occurrence taxon input is specified and has name " +
                    taxonFirstOccurrenceInput.get().getID());
        }
        firstOccurrenceID = ID;
    }

    public void setLastOccurrenceID(String ID) {
        if (taxonLastOccurrenceInput.get() != null && ! taxonLastOccurrenceInput.get().getID().equals(ID)) {
            throw new RuntimeException(ID + " was attempted to be assigned as the name of the last occurrence taxon " +
                    "for a stratigraphic range for which the last occurrence taxon input is specified and has name " +
                    taxonFirstOccurrenceInput.get().getID());
        }
        lastOccurrenceID = ID;
    }

    public String getFirstOccurrenceID() {
        if (firstOccurrenceID == null && taxonFirstOccurrenceInput.get() != null) {
            firstOccurrenceID = taxonFirstOccurrenceInput.get().getID();
        }
        return firstOccurrenceID;
    }

    public String getLastOccurrenceID() {
        if (lastOccurrenceID == null && taxonLastOccurrenceInput.get() != null) {
            lastOccurrenceID = taxonLastOccurrenceInput.get().getID();
        }
        return lastOccurrenceID;
    }

    public void addNodeNr(SRTree tree, int nodeNr) {
        int nr = nodeNr;
        if (tree.getNode(nodeNr).isFake())
            nr = tree.getNode(nodeNr).getDirectAncestorChild().getNr();
        nodes.add(nr);
    }
}
