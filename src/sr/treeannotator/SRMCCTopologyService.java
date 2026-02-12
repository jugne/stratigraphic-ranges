package sr.treeannotator;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.services.TopologySettingService;
import sr.evolution.tree.SRTree;

import java.io.IOException;
import java.io.PrintStream;

/**
 * SR-specific MCC topology service that uses relationship-based credibility
 * instead of clade-based credibility for finding the Maximum Credibility tree.
 *
 * Implements the approach from section 1.1.3 where the MCC tree is the one
 * that maximizes the product of relationship probabilities.
 *
 * @author Ugne Stolz
 */
public class SRMCCTopologyService implements TopologySettingService {

    protected TreeAnnotator.TreeSet treeSet;
    protected int totalTreesUsed;
    protected RelationshipSystem relationshipSystem;

    @Override
    public Tree setTopology(TreeAnnotator.TreeSet treeSet, PrintStream progressStream, TreeAnnotator annotator) throws IOException {
        progressStream.println("Finding maximum relationship credibility tree for SR trees...");
        this.treeSet = treeSet;
        this.totalTreesUsed = annotator.getTotalTreesUsed();
        return summarizeSRTrees(false, progressStream);
    }

    /**
     * Finds the MCC tree by:
     * 1. Collecting all relationships from all trees
     * 2. Computing posterior probabilities for each relationship
     * 3. Finding the tree with maximum product (or sum of log) of relationship probabilities
     *
     * @param useSumCredibility If true, use sum of probabilities; if false, use sum of log probabilities
     * @param progressStream Stream for progress output
     * @return The MCC tree
     * @throws IOException If tree reading fails
     */
    protected Tree summarizeSRTrees(boolean useSumCredibility, PrintStream progressStream) throws IOException {
        // Phase 1: Collect all relationships and their counts
        progressStream.println("Phase 1: Collecting relationships from " + totalTreesUsed + " trees...");
        relationshipSystem = new RelationshipSystem();

        treeSet.reset();
        int counter = 0;
        while (treeSet.hasNext()) {
            Tree tree = treeSet.next();
            if (tree instanceof SRTree) {
                relationshipSystem.add((SRTree) tree);
            } else {
                throw new IOException("Tree is not an SR tree. This service requires SR trees.");
            }
            counter++;

            if (counter % 1000 == 0) {
                progressStream.print(".");
                progressStream.flush();
            }
        }
        progressStream.println();
        progressStream.println("Collected relationships from " + counter + " trees.");

        // Calculate posterior probabilities
        relationshipSystem.calculatePosteriorProbabilities(totalTreesUsed);

        // Optional: Print summary of relationships
        if (progressStream != System.err) {
            progressStream.println("\nRelationship Summary:");
            progressStream.println(relationshipSystem.getSummary());
        }

        // Phase 2: Find tree with highest credibility
        progressStream.println("\nPhase 2: Finding maximum credibility tree...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        Tree bestTree = null;
        double bestScore = Double.NEGATIVE_INFINITY;

        int reported = 0;
        counter = 0;
        treeSet.reset();

        while (treeSet.hasNext()) {
            Tree tree = treeSet.next();

            if (!(tree instanceof SRTree)) {
                throw new IOException("Tree is not an SR tree. This service requires SR trees.");
            }

            double score = scoreTree((SRTree) tree, useSumCredibility);

            if (score > bestScore) {
                bestTree = tree;
                bestScore = score;
            }

            // Progress reporting
            while (reported < 61 && 1000.0 * reported < 61000.0 * (counter + 1) / totalTreesUsed) {
                progressStream.print("*");
                reported++;
                progressStream.flush();
            }
            counter++;
        }

        progressStream.println();
        progressStream.println();

        if (useSumCredibility) {
            progressStream.println("Highest Sum Relationship Credibility: " + bestScore);
        } else {
            progressStream.println("Highest Log Relationship Credibility: " + bestScore);
        }

        if (bestTree != null) {
            bestTree.initAndValidate();
        }

        return bestTree;
    }

    /**
     * Scores a tree based on its relationship credibility.
     *
     * @param tree The SR tree to score
     * @param useSumCredibility If true, use sum of probabilities; if false, use sum of log probabilities
     * @return The credibility score
     */
    public double scoreTree(SRTree tree, boolean useSumCredibility) {
        if (useSumCredibility) {
            return relationshipSystem.getSumRelationshipCredibility(tree);
        } else {
            return relationshipSystem.getLogRelationshipCredibility(tree);
        }
    }

    @Override
    public String getServiceName() {
        return "SR-MCC";
    }

    @Override
    public String getDescription() {
        return "Maximum relationship credibility tree for stratigraphic-range trees";
    }

    /**
     * Gets the relationship system used for this analysis.
     *
     * @return The relationship system
     */
    public RelationshipSystem getRelationshipSystem() {
        return relationshipSystem;
    }
}
