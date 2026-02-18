package sr.evolution.tree;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.StateNode;
import sr.evolution.sranges.StratigraphicRange;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexandra Gavryushkina
 * @author Ugne Stolz
 * A labeled oriented tree on stratigraphic ranges under budding speciation.
 * At every internal node contains metadata on its orientation.
 * Branch leading to the left child represents the ancestral species
 * and the branch leading to the right child represents the descendant species.
 * Sampled ancestors are always attached from left, with the appropriate metadata recorded.
 */
public class SRTree extends Tree implements TreeInterface {

    public Input<List<StratigraphicRange>> stratigraphicRangeInput = new Input<>("stratigraphicRange", "all stratigraphic ranges", new ArrayList<>());

    public Input<Tree> treeInput = new Input<>("tree", "tree to start with");

    protected ArrayList<StratigraphicRange> sRanges;
    protected ArrayList<StratigraphicRange> storedSRanges;

    /**
     * Initializes and validates the object, assigns the tree if provided,
     * and initializes the stratigraphic ranges.
     */
    @Override
    public void initAndValidate() {
        if (treeInput.get() != null) {
            assignFromWithoutID(treeInput.get());
        }

        super.initAndValidate();
    }

    /**
     * Initializes the stratigraphic ranges based on the input or
     * inferred from the provided tree structure.
     */
    protected void initSRanges() {
        if (stratigraphicRangeInput.get().size() != 0) {
            sRanges = (ArrayList<StratigraphicRange>) stratigraphicRangeInput.get();
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
//                    if (!node.isDirectAncestor()) {
//                        throw new RuntimeException("The first occurrence always has to be a sampled ancestor but " +
//                                node.getID() + " is not a sampled ancestor.");
//                    }
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
            if (!lastRanges.isEmpty()) {
                throw new RuntimeException("There are taxa with last occurrence only " + lastRanges.toString() );
            }
            for (StratigraphicRange range:firstRanges) {
                range.makeSingleFossilRange();
            }
            sRanges.addAll(firstRanges);
        }

        initStoredRanges();
    }

    /**
     * Initializes the stored stratigraphic ranges,
     * based on the current stratigraphic ranges.
     */
    public void initStoredRanges() {
        storedSRanges = new ArrayList<>();
        for (int i=0; i<sRanges.size(); i++) {
            StratigraphicRange range_src = sRanges.get(i);
            StratigraphicRange range_sink = new StratigraphicRange();
            ArrayList<Integer> nodeNrs_src = (ArrayList<Integer>) range_src.getNodeNrs();
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
    }

    /**
     * as assignFrom, but only copy tree structure *
     */
    @Override
    public void assignFromFragile(final StateNode other) {
        // invalidate cache
        postCache = null;

        final Tree tree = (Tree) other;
        sr.util.Tools.orientateNodeChildren(tree.getRoot().getNr(), tree);
        if (m_nodes == null) {
            initArrays();
        }
        root = m_nodes[tree.getRoot().getNr()];
        final Node[] otherNodes = tree.getNodesAsArray();
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
            for (StratigraphicRange range : this.getSRanges()) {
                int firstNr = range.getNodeNrs().get(0);
                Node first = this.getNode(firstNr);
                first = first.isDirectAncestor() ? first.getParent() : first;
                Node last = this.getNode(range.getNodeNrs().get(range.getNodeNrs().size() - 1));
                while (!range.isSingleFossilRange() && firstNr != last.getNr()) {
                        if (first.isLeaf())
                            throw new RuntimeException("Error when restoring. Please contact the package developers if this occurs!");
                        int nr = first.getLeft().isFake() ? first.getLeft().getDirectAncestorChild().getNr() : first.getLeft().getNr();
                        if (nr == last.getNr())
                            break;
                        range.addNodeNrAfter(this, first.getNr(), nr);
                        firstNr = nr;
                        first = this.getNode(nr).isDirectAncestor() ? this.getNode(nr).getParent() : this.getNode(nr);
                }
            }
        }
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
     * Stores the current state of the tree and stratigraphic ranges.
     */
    @Override
    protected void store() {
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

    /**
     * Populates the tree from XML file data.
     * @param node XML node containing tree data.
     */
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        final String newick = node.getTextContent();
        final TreeParser parser = new TreeParser();
        try {
            parser.thresholdInput.setValue(1e-10, parser);
        } catch (Exception e1) {
            e1.printStackTrace();
        }
        try {
            parser.offsetInput.setValue(0, parser);
            setRoot(parser.parseNewick(newick));
        } catch (Exception e) {
            e.printStackTrace();
        }
        initArrays();
        orientateTree();
//        for (StratigraphicRange range : this.getSRanges()) {
//            int firstNr = range.getNodeNrs().get(0);
//            Node first = this.getNode(range.getNodeNrs().get(0));
//            first = first.isDirectAncestor() ? first.getParent() : first;
//            Node last = this.getNode(range.getNodeNrs().get(range.getNodeNrs().size() - 1));
//            while (!range.isSingleFossilRange() && firstNr !=last.getNr()){
//                if (first.isLeaf())
//                    System.out.println();
//                int nr = first.getLeft().isFake() ? first.getLeft().getDirectAncestorChild().getNr(): first.getLeft().getNr();
//                range.addNodeNrAfter(this, first.getNr(), nr);
//                firstNr = nr;
//                first = this.getNode(nr).isDirectAncestor() ? this.getNode(nr).getParent(): this.getNode(nr);
//            }
//        }
    }

    @Override
    public String toString() {
        StackTraceElement[] ste = Thread.currentThread().getStackTrace();
        if (ste[2].getMethodName().equals("toXML")) {
            addOrientationMetadata();
        }
        return root.toString();
    }

    /**
     * Restores the tree and stratigraphic ranges from the stored state.
     */
    @Override
    public void restore() {
        super.restore();

        ArrayList<StratigraphicRange> tmp_ranges = storedSRanges;
        storedSRanges = sRanges;
        sRanges = tmp_ranges;
    }

    /**
     * Orientates the tree according to stored (ancestor/donor - descendant/recipient) metadata.
     */
    public void orientateTree() {
        sr.util.Tools.orientateNodeChildren(getRoot().getNr(), this);
    }

    /**
     * Returns the stratigraphic ranges of the tree.
     * @return the stratigraphic ranges of the tree.
     */
    public ArrayList<StratigraphicRange> getSRanges() {
        return  sRanges;
    }
    public Integer getNonSingleSRangeCount() {
        int count = 0;
        for (StratigraphicRange range: sRanges) {
            if (!range.isSingleFossilRange()) {
                count++;
            }
        }
        return  count;
    }

    /**
     * @return the list of internal node numbers in the stratigraphic ranges.
     */
    public ArrayList<Integer> getSRangesInternalNodeNrs() { // nodes that do not represent the first occurrences, does not include nodes
        // from single fossil range
        ArrayList<Integer> internalNodeNrs = new ArrayList<>();
        for (StratigraphicRange range: sRanges) {
            internalNodeNrs.addAll(range.getInternalNodeNrs(this));
        }
        return internalNodeNrs;
    }

    /**
     *
     * @param taxonName
     * @return  the stratigraphic range containing the taxon with the given name,
     *          or null if no such range exists.
     */
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

    /**
     * Retrieves the stratigraphic range to which a node belongs.
     */
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

    /**
     * @param node1Nr
     * @param node2Nr
     * @return  the stratigraphic range to which both nodes belong,
     *          or null if no such range exists.
     */
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


    /**
     *
     * @param node1Nr
     * @param node2Nr
     * @return  true if both nodes belong to the same stratigraphic range,
     *          false otherwise.
     */
    public boolean belongToSameSRange(int node1Nr, int node2Nr) {
        if (getSharedRange(node1Nr, node2Nr)!=null)
            return true;
        return false;
    }

    /**
     * Adds the orientation (donor-recipient) metadata.
     */
    public void addOrientationMetadata() {
        addOrientationMetadataNode(getRoot().getNr());
    }

    /**
     * Add orientation metadata to each node. Left (0) child is always a ancestor species.
     * Right (1) child is always a descendant species. Allows for: - tree state to be stored
     * and restored from file even when BEAST applies sorting. - metadata can be
     * used to color the output tree lineages for donors and recipients.
     *
     * @param subtreeRootNr the node number
     */
    private void addOrientationMetadataNode(int subtreeRootNr) {
        Node subRoot = this.getNode(subtreeRootNr);
        if (subRoot.isRoot()) {
            subRoot.metaDataString = "orientation=ancestor";
        }

        if (!subRoot.isLeaf()) {
            if (subRoot.isFake()) {
                subRoot.getLeft().metaDataString = subRoot.metaDataString;
                subRoot.getRight().metaDataString = subRoot.metaDataString;
            } else if (subRoot.getChildCount()==1){
                subRoot.getLeft().metaDataString = subRoot.metaDataString;
            } else {
                subRoot.getLeft().metaDataString = "orientation=ancestor";
                subRoot.getRight().metaDataString = "orientation=descendant";
            }

            addOrientationMetadataNode(subRoot.getLeft().getNr());
            if(subRoot.getChildCount()!=1){
                addOrientationMetadataNode(subRoot.getRight().getNr());
            }
        }
    }

    /**
     * @param sample chain sample number
     * @param out     log stream
     */
    @Override
    public void log(long sample, PrintStream out) {
        SRTree tree = (SRTree) getCurrent();
        out.print("tree STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        //tree.getRoot().sort();
        final String newick = ((SRNode) tree.getRoot()).toShortNewickForLog(false);
        out.print(newick);
        out.print(";");
    }

    /**
     * deep copy, returns a completely new tree
     *
     * @return a deep copy of this beast.tree.
     */
    @Override
    public Tree copy() {
        SRTree tree = new SRTree();
        tree.setID(getID());
        tree.index = index;
        tree.root = root.copy();
        tree.nodeCount = nodeCount;
        tree.internalNodeCount = internalNodeCount;
        tree.leafNodeCount = leafNodeCount;
        return tree;
    }

}
