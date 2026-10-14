package sr.treeannotator;

import sr.evolution.tree.SRTree;

/**
 * Common interface for the methods used to summarize a posterior sample of SR trees and to
 * find / annotate a maximum credibility tree.
 *
 * Two implementations exist:
 * <ul>
 *     <li>{@link RelationshipSystem} - relationship-based credibility (ancestry + orientation).</li>
 *     <li>{@link CladeSystem} - clade-based credibility using the new SR clade definition.</li>
 * </ul>
 *
 * @author Alexandra Gavryushkina
 */
public interface TreeSummarizer {

    /**
     * Collects the features (relationships or clades) from a tree.
     *
     * @param tree           the SR tree to process
     * @param collectHeights if true, also collect node heights for statistics
     */
    void add(SRTree tree, boolean collectHeights);

    /**
     * Computes posterior probabilities for all collected features.
     *
     * @param totalTrees the total number of trees in the posterior sample
     */
    void calculatePosteriorProbabilities(int totalTrees);

    /**
     * @return the log credibility of a tree (sum of log feature probabilities).
     */
    double getLogCredibility(SRTree tree);

    /**
     * @return the sum credibility of a tree (sum of feature probabilities).
     */
    double getSumCredibility(SRTree tree);

    /**
     * @return a human-readable summary of all collected features.
     */
    String getSummary();

    /**
     * Annotates the MCC tree with feature probabilities and height statistics.
     *
     * @param tree                       the MCC tree to annotate
     * @param includeRelationshipDetails if true, include detailed taxa annotations
     */
    void annotateMCCTree(SRTree tree, boolean includeRelationshipDetails);
}
