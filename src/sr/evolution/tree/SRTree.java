package sr.evolution.tree;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import sr.evolution.sranges.StratigraphicRange;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexandra Gavryushkina
 * A labeled oriented tree on stratigraphic ranges under budding speciation.
 * At every internal node, the branch leading to the left child represents the ancestral species
 * and the branch leading to the right child represents the descendant species.
 * Branches starting with sampled ancestors inheret the orientation from the parental branches:
 * if a node is fake then the non-direct ancestor child gets the same left/right orientation as the fake node.
 * If the root is a fake node then the non-direct ancestor child is always left.
 */
public class SRTree extends Tree implements TreeInterface {

    public Input<List<StratigraphicRange>> stratigraphicRangeInput = new Input<>("stratigraphicRange", "all stratigraphic ranges", new ArrayList<>());

    public Input<Tree> treeInput = new Input<>("tree", "tree to start with");

    protected ArrayList<StratigraphicRange> sRanges;
    protected ArrayList<StratigraphicRange> storedSRanges;

    @Override
    public void initAndValidate() {
        if (treeInput.get() != null) {
            assignFromWithoutID(treeInput.get());
        }

        super.initAndValidate();
//        addOrientationMetadata(root.getNr());
    }

    protected void initSRanges() {
        if (stratigraphicRangeInput.get().size() != 0) {
            sRanges = (ArrayList) stratigraphicRangeInput.get();
            List<Node> externalNodes = getExternalNodes();
            List<Node> unusedNodes = getExternalNodes();
            for (StratigraphicRange range:sRanges) {
                range.removeAllNodeNrs();
                for (Node node:externalNodes) {
                    if(node.getID().equals(range.getFirstOccurrenceID()) && !range.isSingleFossilRange()) {
                        if (!node.isDirectAncestor()) {
                            throw new RuntimeException("The first occurrence always has to be a sampled ancestor but " +
                                    range.getFirstOccurrenceID() + " is not a sampled ancestor. Something went wrong in " +
                                    "initializing the stratigraphic range tree."  );
                        }
                        range.setFirstOccurrenceNodeNr(this, node.getNr());
                        unusedNodes.remove(node);
                    }
                    if (node.getID().equals(range.getLastOccurrenceID())) {
                        range.setLastOccurrenceNodeNr(this, node.getNr());
                        unusedNodes.remove(node);
                    }
                }
                range.initAndValidate();
            }
            for (Node n : unusedNodes){
                StratigraphicRange tmpRange = new StratigraphicRange();
                tmpRange.setFirstOccurrenceNodeNr(this, n.getNr());
                tmpRange.setFirstOccurrenceID(n.getID());
                tmpRange.setLastOccurrenceNodeNr(this, n.getNr());
                tmpRange.setLastOccurrenceID(n.getID());
                tmpRange.initAndValidate();
                sRanges.add(tmpRange);
            }

        } else {
            sRanges = new ArrayList<>();
            ArrayList<StratigraphicRange> firstRanges = new ArrayList<>();
            ArrayList<StratigraphicRange> lastRanges = new ArrayList<>();
            for (Node node:getExternalNodes()) {
                String ID = node.getID();
                String IDwithoutPrefix = ID;
                String prefix = "";
                int i=ID.lastIndexOf('_');
                if (i>0) {
                    IDwithoutPrefix = ID.substring(0,i);
                    prefix = ID.substring(i+1);
                }
                if (prefix.equals("first")) {
                    if (!node.isDirectAncestor()) {
                        throw new RuntimeException("The first occurrence always has to be a sampled ancestor but " +
                                node.getID() + " is not a sampled ancestor.");
                    }
                    boolean found = false;
                    for (StratigraphicRange candidateRange:lastRanges) {
                        if (candidateRange.getID().equals(IDwithoutPrefix)) {
                            candidateRange.setFirstOccurrenceID(ID);
                            candidateRange.setFirstOccurrenceNodeNr(this, node.getNr());
                            sRanges.add(candidateRange);
                            lastRanges.remove(candidateRange);
                            found=true;
                            break;
                        }
                    }
                    if (!found) {
                        StratigraphicRange range = new StratigraphicRange();
                        range.setID(IDwithoutPrefix);
                        range.setFirstOccurrenceID(ID);
                        if (node.isDirectAncestor()) {
                            range.setFirstOccurrenceNodeNr(this, node.getNr());
                        } else {
                            range.setFirstOccurrenceNodeNr(this, node.getNr());
                        }
                        firstRanges.add(range);
                    }
                } else {
                    boolean found = false;
                    for (StratigraphicRange candidateRange:firstRanges) {
                        if (candidateRange.getID().equals(IDwithoutPrefix)) {
                            if (!prefix.equals("last")) {
                                throw new RuntimeException("Taxa " + candidateRange.getFirstOccurrenceID() + " and " +
                                        ID  + " are found in the tree. If " + ID + " is the last occurrence then add " +
                                        "_last at the end.");
                            }
                            candidateRange.setLastOccurrenceID(ID);
                            candidateRange.setLastOccurrenceNodeNr(this, node.getNr());
                            sRanges.add(candidateRange);
                            firstRanges.remove(candidateRange);
                            found=true;
                            break;
                        }
                    }
                    if (!found) {
                        StratigraphicRange range = new StratigraphicRange();
                        range.setID(IDwithoutPrefix);
                        range.setLastOccurrenceID(ID);
                        range.setLastOccurrenceNodeNr(this, node.getNr());
                        lastRanges.add(range);
                    }
                }
            }
            if (!firstRanges.isEmpty()) {
                throw new RuntimeException("There are taxa with first occurrence only " + firstRanges.toString() +". " +
                        "Single fossil ranges can not have _first at the end." );
            }
            for (StratigraphicRange range:lastRanges) {
                range.makeSingleFossilRange();
            }
            sRanges.addAll(lastRanges);
        }

        initStoredRanges();
    }

    public void initStoredRanges() {
        storedSRanges = new ArrayList<>();
        for (int i=0; i<sRanges.size(); i++) {
            StratigraphicRange range_src = sRanges.get(i);
            StratigraphicRange range_sink = new StratigraphicRange();
            ArrayList<Integer> nodeNrs_src = (ArrayList) range_src.getNodeNrs();
            for (int j=0; j<nodeNrs_src.size(); j++) {
                range_sink.addNodeNr(this, nodeNrs_src.get(j));
            }
            range_sink.setFirstOccurrenceID(range_src.getFirstOccurrenceID());
            range_sink.setLastOccurrenceID(range_src.getLastOccurrenceID());
            storedSRanges.add(range_sink);
        }

    }

    /**
     * copy of all values from existing tree *
     */
    @Override
    public void assignFrom(final StateNode other) {
        final Tree tree = (Tree) other;
        final Node[] nodes = new Node[tree.getNodeCount()];//tree.getNodesAsArray();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            nodes[i] = newNode();
        }
        setID(tree.getID());
        //index = tree.index;
        root = nodes[tree.getRoot().getNr()];
        root.assignFrom(nodes, tree.getRoot());
        root.setParent(null);
        nodeCount = tree.getNodeCount();
        internalNodeCount = tree.getInternalNodeCount();
        leafNodeCount = tree.getLeafNodeCount();
        initArrays();
        initSRanges();
//        addOrientationMetadata(tree.getRoot().getNr());
//        addOrientationMetadata(root.getNr());

    }

    /**
     * as assignFrom, but only copy tree structure *
     */
    @Override
    public void assignFromFragile(final StateNode other) {
        // invalidate cache
        postCache = null;

        final SRTree tree = (SRTree) other;
//        addOrientationMetadata(tree.getRoot().getNr());
        if (m_nodes == null) {
            initArrays();
        }
        root = m_nodes[tree.getRoot().getNr()];
        final Node[] otherNodes = tree.m_nodes;
        final int rootNr = root.getNr();
        assignFrom(0, rootNr, otherNodes);
        root.setHeight(otherNodes[rootNr].getHeight());
        root.setParent(null);
        if (otherNodes[rootNr].getLeft() != null) {
            root.setLeft(m_nodes[otherNodes[rootNr].getLeft().getNr()]);

        } else {
            root.setLeft(null);
        }
        if (otherNodes[rootNr].getRight() != null) {
            root.setRight(m_nodes[otherNodes[rootNr].getRight().getNr()]);
        } else {
            root.setRight(null);
        }
        assignFrom(rootNr + 1, nodeCount, otherNodes);
        if(stratigraphicRangeInput.get()!= null) {
            initSRanges();
        }
//        addOrientationMetadata(root.getNr());
    }

    /**
     * helper to assignFromFragile *
     */
    private void assignFrom(final int start, final int end, final Node[] otherNodes) {
        for (int i = start; i < end; i++) {
            Node sink = m_nodes[i];
            Node src = otherNodes[i];
            sink.setHeight(src.getHeight());
            sink.setParent(m_nodes[src.getParent().getNr()]);
            if (src.getLeft() != null) {
                sink.setLeft(m_nodes[src.getLeft().getNr()]);
                if (src.getRight() != null) {
                    sink.setRight(m_nodes[src.getRight().getNr()]);
                } else {
                    sink.setRight(null);
                }
            }
        }
    }

    /**
     * StateNode implementation *
     */
    @Override
    protected void store() {
//        addOrientationMetadata(root.getNr());
        storeNodes(0, nodeCount);
        storedRoot = m_storedNodes[root.getNr()];
        for (StratigraphicRange range_src:sRanges) {
            int index = sRanges.indexOf(range_src);
            StratigraphicRange range_sink = storedSRanges.get(index);
            range_sink.removeAllNodeNrs();
            for (int i=0; i< range_src.getNodeNrs().size(); i++) {
                int nodeNr = range_src.getNodeNrs().get(i);
                range_sink.addNodeNr(this, nodeNr);
            }
        }
    }

    /**
     * Stores nodes with index i, for start <= i < end
     * (i.e. including start but not including end)
     *
     * @param start the first index to be stored
     * @param end   nodes are stored up to but not including this index
     */
    private void storeNodes(final int start, final int end) {
        // Use direct members for speed (we are talking 5-7% or more from total time for large trees :)
        for (int i = start; i < end; i++) {
            final SRNode sink = (SRNode)m_storedNodes[i];
            final SRNode src = (SRNode)m_nodes[i];
            sink.setHeight(src.getHeight());

            if ( src.getParent() != null ) {
                sink.setParent(m_storedNodes[src.getParent().getNr()]);
            } else {
                // currently only called in the case of sampled ancestor trees
                // where root node is not always last in the list
                sink.setParent(null);
            }

            final List<Node> children = sink.getChildrenMutable();
            final List<Node> srcChildren = src.getChildren();

            if( children.size() == srcChildren.size() ) {
                // shave some more time by avoiding list clear and add
                for (int k = 0; k < children.size(); ++k) {
                    final SRNode srcChild = (SRNode)srcChildren.get(k);
                    // don't call addChild, which calls  setParent(..., true);
                    final Node c = m_storedNodes[srcChild.getNr()];
                    c.setParent(sink);
                    children.set(k, c);
                }
            } else {
                children.clear();
                //sink.removeAllChildren(false);
                for (final Node srcChild : srcChildren) {
                    // don't call addChild, which calls  setParent(..., true);
                    final Node c = m_storedNodes[srcChild.getNr()];
                    c.setParent(sink);
                    children.add(c);
                    //sink.addChild(c);
                }
            }
        }
    }

    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        super.fromXML(node);
        orientateTree();
    }

    @Override
    public void restore() {

        super.restore();

/*        // necessary for sampled ancestor trees
        nodeCount = m_storedNodes.length;

        final Node[] tmp = m_storedNodes;
        m_storedNodes = m_nodes;
        m_nodes = tmp;
        root = m_nodes[storedRoot.getNr()];

        // necessary for sampled ancestor trees,
        // we have the nodes, no need for expensive recursion
        leafNodeCount = 0;
        for( Node n : m_nodes ) {
            leafNodeCount += n.isLeaf() ? 1 : 0;
        }

        //leafNodeCount = root.getLeafNodeCount();

        hasStartedEditing = false;

        for( Node n : m_nodes ) {
            n.se
            n.isDirty = Tree.IS_CLEAN;
        }

        postCache = null;*/

        ArrayList<StratigraphicRange> tmp_ranges = storedSRanges;
        storedSRanges = sRanges;
        sRanges=tmp_ranges;
//        orientateTree();
    }

    /**
     * Orientates the tree according to stored (donor-recipient) metadata.
     */
    public void orientateTree() {
        orientateNodeChildren(getRoot().getNr());
    }

    /**
     * Orientate node children depending on stored metadata.
     *
     * @param subtreeRootNr the none number
     */
    private void orientateNodeChildren(int subtreeRootNr) {
        Node subTreeRoot = this.getNode(subtreeRootNr);
        if (!subTreeRoot.isLeaf()) {
            if(subTreeRoot.getChild(0).metaDataString == null)
                return;
            if ((!subTreeRoot.isFake() && !subTreeRoot.getChild(0).metaDataString.contains("orientation=donor"))
                    || (subTreeRoot.isFake() && subTreeRoot.getChild(1).getHeight() != subTreeRoot.getHeight())) {
                Node left = subTreeRoot.getChild(1);
                Node right = subTreeRoot.getChild(0);

                subTreeRoot.removeAllChildren(false);

//                subTreeRoot.setLeft(left);
//                subTreeRoot.setRight(right);
				subTreeRoot.addChild(left);
				subTreeRoot.addChild(right);
            }

            orientateNodeChildren(subTreeRoot.getChild(0).getNr());
            orientateNodeChildren(subTreeRoot.getChild(1).getNr());
        }
    }

    /**
     * Add orientation metadata to each node. Left (0) child is always a donor.
     * Right (1) child is always a recipient. Allows for: - tree state to be stored
     * and restored from file even when BEAST applies sorting. - metadata can be
     * used to color the output tree lineages for donors and recipients.
     *
     * @param subtreeRootNr the node number
     */
    public void addOrientationMetadata(int subtreeRootNr) {
        Node subRoot = this.getNode(subtreeRootNr);
        if (subRoot.isRoot()) {
            subRoot.metaDataString = "orientation=donor";
        }

        if (!subRoot.isLeaf()) {
            if (subRoot.isFake()) {
                subRoot.getLeft().metaDataString = subRoot.metaDataString;
                subRoot.getRight().metaDataString = subRoot.metaDataString;
            } else if (subRoot.getChildCount()==1){
                subRoot.getLeft().metaDataString = subRoot.metaDataString;
            } else {
                subRoot.getLeft().metaDataString = "orientation=donor";
                subRoot.getRight().metaDataString = "orientation=recipient";
            }

            addOrientationMetadata(subRoot.getLeft().getNr());
            if(subRoot.getChildCount()!=1){
                addOrientationMetadata(subRoot.getRight().getNr());
            }
        }
    }


    // SRange methods:

    public ArrayList<StratigraphicRange> getSRanges() {
        return  sRanges;
    }

    public ArrayList<Integer> getSRangesInternalNodeNrs() { // nodes that do not represent the first occurrences, does not include nodes
        // from single fossil range
        ArrayList<Integer> internalNodeNrs = new ArrayList<>();
        for (StratigraphicRange range: sRanges) {
            internalNodeNrs.addAll(range.getInternalNodeNrs(this));
        }
        return internalNodeNrs;
    }

    public StratigraphicRange sRangesContainsID(String taxonName) {
        for (StratigraphicRange range:sRanges) {
            if (range.getFirstOccurrenceID().equals(taxonName)) {
                return range;
            }
            if (range.getLastOccurrenceID().equals(taxonName)) {
                return range;
            }
        }
        return null;
    }

    public StratigraphicRange getRangeOfNode(Node node) {
        int nodeNr = node.getNr();
        if (node.isFake())
            nodeNr = node.getDirectAncestorChild().getNr();
        for (StratigraphicRange candidate_range:sRanges) {
            if (candidate_range.containsNodeNr(this, nodeNr)) {
                return candidate_range;
            }
        }
        return null;
    }

    public StratigraphicRange getSharedRange(int node1Nr, int node2Nr){
        if (m_nodes[node1Nr].isFake())
            node1Nr = m_nodes[node1Nr].getDirectAncestorChild().getNr();
        if (m_nodes[node2Nr].isFake())
            node2Nr = m_nodes[node2Nr].getDirectAncestorChild().getNr();
        for (StratigraphicRange range:sRanges) {
            if (range.containsNodeNr(this, node1Nr) && range.containsNodeNr(this, node2Nr)) {
                return range;
            }
        }
        return null;
    }


    public boolean belongToSameSRange(int node1Nr, int node2Nr) {
        if (getSharedRange(node1Nr, node2Nr)!=null)
            return true;
        return false;
    }

    @Override
    public void log(long sample, PrintStream out) {
        Tree tree = (Tree) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        //tree.getRoot().sort();
        final String newick = tree.getRoot().toShortNewick(false);
        out.print(newick);
        out.print(";");
    }

}
