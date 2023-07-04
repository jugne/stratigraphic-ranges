package sr.evolution.tree;

import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.util.HeapSort;
import beast.base.util.Randomizer;
import sr.evolution.sranges.StratigraphicRange;

import java.util.*;

/**
 * @author Alexandra Gavryushkina with the main part copied from RandomTree
 */
public class RandomSRangeTree extends SRTree implements StateNodeInitialiser {

    final public Input<Alignment> taxaInput = new Input<>("taxa", "set of taxa to initialise tree specified by alignment");

    final public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel", "population function for generating coalescent", Input.Validate.REQUIRED);

    // total nr of taxa
    int nrOfTaxa;

    List<Integer>[] children;

    Set<String> taxa;

    // number of the next internal node, used when creating new internal nodes
    int nextNodeNr;

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

    @SuppressWarnings({"rawtypes", "unchecked"})
    private void swap(final List list, final int i, final int j) {
        final Object tmp = list.get(i);
        list.set(i, list.get(j));
        list.set(j, tmp);
    }

    @SuppressWarnings("unchecked")
    @Override
    public void initStateNodes() {
        if (taxaInput.get() != null) {
            taxa.addAll(taxaInput.get().getTaxaNames());
        } else {
            taxa.addAll(m_taxonset.get().asStringList());
        }


        final PopulationFunction popFunction = populationFunctionInput.get();

        simulateTree(taxa, popFunction);

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
        setNodesNrs(root, 0, new int[1], taxonToNR);

        initArrays();

        if (m_initial.get() != null) {
            m_initial.get().assignFromWithoutID(this);
        }
//        addOrientationMetadata(this.root.getNr());
    }

    private int setNodesNrs(final Node node, int internalNodeCount, int[] n, Map<String,Integer> initial) {
        if( node.isLeaf() )  {
            if( initial != null ) {
                node.setNr(initial.get(node.getID()));
            } else {
                node.setNr(n[0]);
                n[0] += 1;
            }
        } else {
            for (final Node child : node.getChildren()) {
                internalNodeCount = setNodesNrs(child, internalNodeCount, n, initial);
            }
            node.setNr(nrOfTaxa + internalNodeCount);
            internalNodeCount += 1;
        }
        return internalNodeCount;
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {

    }

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

}
