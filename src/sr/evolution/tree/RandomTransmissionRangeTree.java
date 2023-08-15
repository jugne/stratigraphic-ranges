package sr.evolution.tree;

import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.evolution.tree.coalescent.RandomTree;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.HeapSort;
import beast.base.util.Randomizer;
import sr.evolution.sranges.StratigraphicRange;

import java.util.*;

public class RandomTransmissionRangeTree extends SRTree implements StateNodeInitialiser  {

    public Input<List<TaxonSet>> taxonsetsInput = new Input<>("patientTaxonSets",
            "a separate list of taxa for samples collected from the same patient", new ArrayList<>());

    final public Input<Alignment> taxaInput = new Input<>("taxa", "set of taxa to initialise tree specified by alignment");

    final public Input<RealParameter> birthRate = new Input<>("birthRate",
            "Tree prior birth rate to initialize");

    int nrOfTaxa;

    Set<String> taxa;

    @Override
    public void initAndValidate() {

        taxa = new LinkedHashSet<>();
        if (taxaInput.get() != null) {
            taxa.addAll(taxaInput.get().getTaxaNames());
        } else {
            taxa.addAll(m_taxonset.get().asStringList());
        }

        nrOfTaxa = taxa.size();

        initStateNodes();
        super.initAndValidate();
        initStoredRanges();
        this.store();
    }
    int nextNodeNr;

    /**
     * Simulates a coalescent tree, given a taxon list.
     *
     * @param taxa         the set of taxa to simulate a coalescent tree between
     * @param demoFunction the demographic function to use
     */
    public void simulateTree(final Set<String> taxa, final PopulationFunction demoFunction) {
        if (taxa.size() == 0)
            return;

        String msg = "Failed to generate a random tree (probably a bug).";
        for (int attempts = 0; attempts < 1000; ++attempts) {

            nextNodeNr = nrOfTaxa;
            final Set<Node> candidates = new LinkedHashSet<>();
            int i = 0;
            for (String taxon : taxa) {
                final Node node = newNode();
                node.setNr(i);
                node.setID(taxon);
                node.setHeight(0.0);
                candidates.add(node);
                i += 1;
            }

            if (m_initial.get() != null) {
                processCandidateTraits(candidates, m_initial.get().m_traitList.get());
            } else {
                processCandidateTraits(candidates, m_traitList.get());
            }

            root = simulateCoalescent(candidates, demoFunction);

            return;

        }
        throw new RuntimeException(msg);
    }

    /**
     * Apply traits to a set of nodes.
     * @param candidates List of nodes
     * @param traitSets List of TraitSets to apply
     */
    private void processCandidateTraits(Set<Node> candidates, List<TraitSet> traitSets) {
        for (TraitSet traitSet : traitSets) {
            for (Node node : candidates) {
                node.setMetaData(traitSet.getTraitName(), traitSet.getValue(node.getID()));
            }
        }
    }
    private Node simulateCoalescent(final Set<Node> candidates, final PopulationFunction demoFunction) {
        final List<Node> remainingCandidates = new ArrayList<>();
        List<Node> newNodes = new ArrayList<>();
        List<String> firstOccurrenceTaxonNames = new ArrayList<>();
        List<Node> lastOccurrenceNodes = new ArrayList<>();

        for(Node node:candidates) {
            String taxonName = node.getID();
            sRanges= (ArrayList) stratigraphicRangeInput.get();
            StratigraphicRange range = sRangesContainsID(taxonName);
            if (range != null && !range.isSingleFossilRange()) {
                if (range.getFirstOccurrenceID().equals(taxonName)) {
                    final Node newNode = newNode();
                    newNode.setHeight(node.getHeight());
                    newNode.setNr(nextNodeNr++);
                    newNode.setLeft(node);
                    node.setParent(newNode);
                    newNodes.add(newNode);
                    remainingCandidates.add(newNode);

                } else {
                    lastOccurrenceNodes.add(node);
                    firstOccurrenceTaxonNames.add(range.getFirstOccurrenceID());
                }
            } else {
                remainingCandidates.add(node);
            }
        }

        for (int i=0; i< lastOccurrenceNodes.size(); i++) {
            Node lastOccurrenceNode = lastOccurrenceNodes.get(i);
            for (Node node:newNodes) {
                if (node.getLeft().getID().equals(firstOccurrenceTaxonNames.get(i))){
                    Node firstOccur = node.getLeft();
                    node.setLeft(lastOccurrenceNode);
                    lastOccurrenceNode.setParent(node);
                    node.setRight(firstOccur);
                    candidates.remove(lastOccurrenceNode);
                }
            }
        }

        if (remainingCandidates.size() == 0) {
            throw new IllegalArgumentException("empty nodes set");
        }

        final List<Node> rootNode = simulateCoalescent(remainingCandidates, demoFunction, 0.0);
        if (rootNode.size() == 1) {
            return rootNode.get(0);
        }

        throw new RuntimeException("failed to generate a random tree!");
    }


    public List<Node> simulateCoalescent(final List<Node> nodes, final PopulationFunction demographic, double currentHeight) {
        // If only one node, return it
        // continuing results in an infinite loop
        if (nodes.size() == 1)
            return nodes;

        final double[] heights = new double[nodes.size()];
        for (int i = 0; i < nodes.size(); i++) {
            heights[i] = nodes.get(i).getHeight();
        }
        final int[] indices = new int[nodes.size()];
        HeapSort.sort(heights, indices);

        // node list
        nodeList.clear();
        activeNodeCount = 0;
        for (int i = 0; i < nodes.size(); i++) {
            nodeList.add(nodes.get(indices[i]));
        }
        setCurrentHeight(currentHeight);

        // get at least two tips
        while (getActiveNodeCount() < 2) {
            currentHeight = getMinimumInactiveHeight();
            setCurrentHeight(currentHeight);
        }

        // simulate coalescent events
        double nextCoalescentHeight = currentHeight
                + PopulationFunction.Utils.getSimulatedInterval(demographic, getActiveNodeCount(), currentHeight);

        while ((nodeList.size() > 1)) {

            if (nextCoalescentHeight >= getMinimumInactiveHeight()) {
                currentHeight = getMinimumInactiveHeight();
                setCurrentHeight(currentHeight);
            } else {
                currentHeight = coalesceTwoActiveNodes(nextCoalescentHeight);
            }

            if (nodeList.size() > 1) {
                // get at least two tips
                while (getActiveNodeCount() < 2) {
                    currentHeight = getMinimumInactiveHeight();
                    setCurrentHeight(currentHeight);
                }
                nextCoalescentHeight = currentHeight
                        + PopulationFunction.Utils.getSimulatedInterval(demographic, getActiveNodeCount(),
                        currentHeight);
            }
        }

        return nodeList;
    }
    /**
     * @return the height of youngest inactive node.
     */
    private double getMinimumInactiveHeight() {
        if (activeNodeCount < nodeList.size()) {
            return (nodeList.get(activeNodeCount)).getHeight();
        } else
            return Double.POSITIVE_INFINITY;
    }

    /**
     * Set the current height.
     * @param height
     */
    private void setCurrentHeight(final double height) {
        while (getMinimumInactiveHeight() <= height) {
            activeNodeCount += 1;
        }
    }

    /**
     * @return the number of active nodes (equate to lineages)
     */
    private int getActiveNodeCount() {
        return activeNodeCount;
    }


    /**
     * Coalesce two nodes in the active list. This method removes the two
     * (randomly selected) active nodes and replaces them with the new node at
     * the top of the active list.
     * @param height
     * @return
     */
    private double coalesceTwoActiveNodes(double height) {
        final int node1 = Randomizer.nextInt(activeNodeCount);
        int node2 = node1;
        while (node2 == node1) {
            node2 = Randomizer.nextInt(activeNodeCount);
        }

        final Node left = nodeList.get(node1);
        final Node right = nodeList.get(node2);

        final Node newNode = newNode();
        newNode.setNr(nextNodeNr++);   // multiple tries may generate an excess of nodes assert(nextNodeNr <= nrOfTaxa*2-1);
        newNode.setHeight(height);
        newNode.setLeft(left);
        left.setParent(newNode);
        newNode.setRight(right);
        right.setParent(newNode);

        nodeList.remove(left);
        nodeList.remove(right);

        activeNodeCount -= 2;

        nodeList.add(activeNodeCount, newNode);

        activeNodeCount += 1;

        if (getMinimumInactiveHeight() < height) {
            throw new RuntimeException(
                    "This should never happen! Somehow the current active node is older than the next inactive node!\n"
                            + "One possible solution you can try is to increase the population size of the population model.");
        }
        return height;
    }


    final private ArrayList<Node> nodeList = new ArrayList<>();
    private int activeNodeCount = 0;


    @Override
    public void initStateNodes() {
        if (taxaInput.get() != null) {
            taxa.addAll(taxaInput.get().getTaxaNames());
        } else {
            taxa.addAll(m_taxonset.get().asStringList());
        }

        final RealParameter birthRateParameter = birthRate.get();
        final Double lambda = (birthRateParameter == null) ? 1.0 : birthRateParameter.getValue();
        final Double initialPopSize = 1.0 / lambda; // scales coalescent tree height inverse to birth rate
        final RealParameter popSize = new RealParameter(initialPopSize.toString());
        final ConstantPopulation pf = new ConstantPopulation();
        pf.setInputValue("popSize", popSize);

        simulateTree(taxa, pf);


        nodeCount = 2 * taxa.size() - 1;
        internalNodeCount = taxa.size() - 1;
        leafNodeCount = taxa.size();

        HashMap<String,Integer> taxonToNR = null;
        // preserve node numbers where possible
        if (m_initial.get() != null) {
            if( leafNodeCount == m_initial.get().getLeafNodeCount() ) {
                // dont ask me how the initial tree is rubbish  (i.e. 0:0.0)
                taxonToNR = new HashMap<>();
                for (Node n : m_initial.get().getExternalNodes()) {
                    taxonToNR.put(n.getID(), n.getNr());
                }
            }
        } else {
            taxonToNR = new HashMap<>();
            String[] taxa = getTaxaNames();
            for(int k = 0; k < taxa.length; ++k) {
                taxonToNR.put(taxa[k], k);
            }
        }
        // multiple simulation tries may produce an excess of nodes with invalid nr's. reset those.
//        setNodesNrs(root, 0, new int[1], taxonToNR);

        initArrays();
        if (!taxonsetsInput.get().isEmpty())
            samePatientSamplingInit();

        if (m_initial.get() != null) {
            m_initial.get().assignFromWithoutID(this);
        }
//        addOrientationMetadata(this.root.getNr());
    }

    private void samePatientSamplingInit() {
        final RealParameter birthRateParameter = birthRate.get();
        final Double lambda = (birthRateParameter == null) ? 1.0 : birthRateParameter.getValue();
        final Double initialPopSize = 1.0 / lambda; // scales coalescent tree height inverse to birth rate
        final RealParameter popSize = new RealParameter(initialPopSize.toString());
        final ConstantPopulation pf = new ConstantPopulation();
        pf.setInputValue("popSize", popSize);

        final RandomTree rnd = new RandomTree();
        rnd.setInputValue("taxa", taxaInput.get());
        if (m_initial.get().hasDateTrait())
            rnd.setInputValue("trait", m_initial.get().getDateTrait());

        rnd.setInputValue("populationModel", pf);
        rnd.setInputValue("populationModel", pf);
        rnd.initAndValidate();

        /////////////

        int nSets = taxonsetsInput.get().size();
        int nSamples = nrOfTaxa;
        List<List<Node>> taxonsets = new ArrayList<List<Node>>(nSets);
        List<Node> singlePatientSamples = new ArrayList<Node>();
        double maxHeight = Double.NEGATIVE_INFINITY;

        // If there are same patient sampling over time supplied
        if (nSets != 0) {
            for (int i = 0; i < nSets; i++)
                taxonsets.add(new ArrayList<Node>());

            for (int idx = 0; idx < nSamples; idx++) {
                Node leaf = this.getNode(idx);
                if (leaf.getParent() != null)
                    leaf.getParent().removeAllChildren(false);
                int i = 0;
                for (TaxonSet t : taxonsetsInput.get()) {
                    List<String> tmp = t.asStringList();
                    if (tmp.contains(leaf.getID())) {
                        taxonsets.get(i).add(leaf);
                        break;
                    } else {
                        i += 1;
                    }
                }
                singlePatientSamples.add(leaf);
                if (leaf.getHeight() > maxHeight)
                    maxHeight = leaf.getHeight();
            }

            for (int ii = 0; ii < nSets; ii++) {
                taxonsets.get(ii).sort(Comparator.comparing(Node::getHeight));
            }
        }

        int id = nSamples;
        for (int s = 0; s < nSets; s++) {
            Node[] tmp = new Node[taxonsets.get(s).size()];
            tmp = taxonsets.get(s).toArray(tmp);
            int nNodes = tmp.length;
            for (int n = 0; n < nNodes - 1; n++) {
                Node parent = new Node();
                parent.setID(Integer.toString(id));
                parent.setNr(id);
                id += 1;

                parent.addChild(tmp[n]);
                parent.addChild(tmp[n + 1]);
                singlePatientSamples.remove(tmp[n]);
                singlePatientSamples.remove(tmp[n + 1]);
                parent.setHeight(tmp[n + 1].getHeight());
                tmp[n] = null;
                tmp[n + 1] = parent;
            }
            singlePatientSamples.add(tmp[nNodes - 1]);
            if (tmp[nNodes - 1].getHeight() > maxHeight)
                maxHeight = tmp[nNodes - 1].getHeight();
        }

        while (singlePatientSamples.size() > 1) {
            Node node1 = singlePatientSamples.get(Randomizer.nextInt(singlePatientSamples.size()));
            singlePatientSamples.remove(node1);
            Node node2 = singlePatientSamples.get(Randomizer.nextInt(singlePatientSamples.size()));
            singlePatientSamples.remove(node2);

            double deltaT = Randomizer.nextExponential(
                    (singlePatientSamples.size() + 2) * (singlePatientSamples.size() + 1) * 0.5 * 1.0 / initialPopSize);
            double coalTime = maxHeight + deltaT;

            Node parent = new Node();
            parent.setID(Integer.toString(id));
            parent.setNr(id);
            id += 1;
            parent.setHeight(coalTime);
            parent.addChild(node1);
            parent.addChild(node2);

            singlePatientSamples.add(parent);
            maxHeight = coalTime;
        }

        Node root = singlePatientSamples.get(0);

        copyTreeStructure(new Tree(root), this);
    }

    // copy the structure of the source tree to the destination tree
    // preserving the leaf node names and numbers
    private void copyTreeStructure(final Tree src, final Tree dst) {
        final Node[] dstNodes = dst.getNodesAsArray();
        final Node[] srcNodes = src.getNodesAsArray();

        final Map<String, Integer> srcTipNumbers = new HashMap<>();

        final int nodeCount = src.getNodeCount();
        final int leafNodeCount = src.getLeafNodeCount();

        // Clear the children of all internal nodes in the destination tree
        for (int nodeNumber = leafNodeCount; nodeNumber < nodeCount; nodeNumber++)
            dstNodes[nodeNumber].removeAllChildren(false);

        // Record the node number of all leaves in the source tree
        for (int nodeNumber = 0; nodeNumber < leafNodeCount; nodeNumber++) {
            final String srcName = srcNodes[nodeNumber].getID();
            srcTipNumbers.put(srcName, nodeNumber);
        }

        // Set the heights of all nodes to match the source height
        for (int nodeNumber = 0; nodeNumber < nodeCount; nodeNumber++) {
            final Node dstNode = dstNodes[nodeNumber];
            Node srcNode;

            // find the corresponding node from the source tree
            if (nodeNumber < leafNodeCount) { // if this is a leaf node
                final String speciesName = dstNode.getID();
                System.out.println(speciesName);
                final int srcTipNumber = srcTipNumbers.get(speciesName);

                srcNode = srcNodes[srcTipNumber];
            } else { // if this is an internal node
                srcNode = srcNodes[nodeNumber];
            }

            // Copy height
            dstNode.setHeight(srcNode.getHeight());

            // Clear and copy metadata
            final Set<String> dstMetaDataNames = dstNode.getMetaDataNames();
            final Set<String> srcMetaDataNames = srcNode.getMetaDataNames();
            final Set<String> srcLengthMetaDataNames = srcNode.getLengthMetaDataNames();

            for (String metaDataName : dstMetaDataNames)
                dstNode.removeMetaData(metaDataName);

            for (String metaDataName : srcMetaDataNames)
                dstNode.setMetaData(metaDataName, srcNode.getMetaData(metaDataName));

            for (String lengthMetaDataName : srcLengthMetaDataNames)
                dstNode.setMetaData(lengthMetaDataName, srcNode.getLengthMetaData(lengthMetaDataName));

            // if this is not the root node, also set the parent and child
            // connections to match the source
            if (nodeNumber != nodeCount - 1) {
                boolean rightChild = false;
                Node siblingNode = null;
                final int parentNumber = srcNode.getParent().getNr();
                final Node parentNode = dstNodes[parentNumber];
                if (srcNode.getParent().getChild(0) == srcNode)
                    rightChild = true;
                if (rightChild && parentNode.getChildCount() > 0) {
                    siblingNode = parentNode.getChild(0);
                    parentNode.removeChild(siblingNode);
                }
                dstNode.setParent(parentNode);
                parentNode.addChild(dstNode);

                if (siblingNode != null) {
                    siblingNode.setParent(parentNode);
                    parentNode.addChild(siblingNode);
                }

            }
        }
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {

    }
}
