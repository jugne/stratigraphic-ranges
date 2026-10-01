package treeannotator;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.tree.Tree;
import beast.base.parser.NexusParser;
import org.junit.Test;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;
import sr.treeannotator.AncestryRelationship;
import sr.treeannotator.OrientationRelationship;
import sr.treeannotator.RelationshipSystem;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * Relationship-based counterpart of {@link MccCladeSystemTest}.
 *
 * <p>Reads the same real posterior tree sample, builds the SR {@link RelationshipSystem}
 * (ancestry and orientation relationships) and the maximum-credibility tree, mirroring the
 * pipeline of {@link sr.treeannotator.SRTreeAnnotator} when the relationship-based summarizer
 * is used (i.e. {@code clades=false}):
 *
 * <ol>
 *   <li>read trees from the NEXUS file and convert them to {@link SRTree}s,</li>
 *   <li>discard a 10% burn-in,</li>
 *   <li>add every remaining tree to the relationship system <b>with heights</b>
 *       ({@code add(tree, true)} - the "keep heights" path used for annotation),</li>
 *   <li>compute posterior relationship probabilities,</li>
 *   <li>pick the tree with the highest log relationship credibility as the MCC tree and
 *       annotate it.</li>
 * </ol>
 *
 * <p>The input ({@code test_data/morph+dna_co_ncc_for_mcc_test.trees}) contains 11 trees, so a
 * 10% burn-in discards 1 tree and the relationship system is built from the remaining 10.</p>
 *
 * <p>The test prints the constructed relationship system, the MCC tree and the per-tree log
 * relationship credibilities. True values are calculated by hand.</p>
 *
 * @author Alexandra Gavryushkina
 */
public class MccRelationshipSystemTest {

    /** Trees file produced by running examples/morph+dna_co_ncc_for_mcc_test.xml. */
    private static final String TREES_FILE = "morph+dna_co_ncc_for_mcc_test.trees";

    /** Burn-in as a percentage of the sample, matching SRTreeAnnotator's default. */
    private static final int BURNIN_PERCENTAGE = 10;

    @Test
    public void testRelationshipSystemAndMccTree() throws Exception {
        // --- read trees (SRTreeAnnotator.readTrees) ----------------------------------
        List<SRTree> trees = readTrees(locateTreesFile());
        assertEquals("expected 11 trees in the sample", 11, trees.size());

        // --- burn-in: 10% of 11 trees -> discard 1, analyse 10 -----------------------
        int burninCount = (BURNIN_PERCENTAGE * trees.size()) / 100;
        List<SRTree> analyzedTrees = trees.subList(burninCount, trees.size());
        int totalTreesUsed = analyzedTrees.size();
        assertEquals("10% burn-in should discard exactly 1 tree", 1, burninCount);
        assertEquals("10 trees should remain after burn-in", 10, totalTreesUsed);

        // --- build the relationship system, keeping heights (add(tree, true)) --------
        RelationshipSystem system = new RelationshipSystem();
        for (SRTree tree : analyzedTrees) {
            system.add(tree, true);
        }
        system.calculatePosteriorProbabilities(totalTreesUsed);

        // --- log relationship credibility of every analysed tree; MCC = the highest
        //     (SRTreeAnnotator step 2). Scores are deterministic (each relationship posterior
        //     is an exact multiple of 1/10), so they can be pinned down exactly. -------------
        SRTree mccTree = null;
        double bestScore = Double.NEGATIVE_INFINITY;
        double[] logCredibilities = new double[totalTreesUsed];
        for (int i = 0; i < totalTreesUsed; i++) {
            SRTree tree = analyzedTrees.get(i);
            double score = system.getLogRelationshipCredibility(tree);
            logCredibilities[i] = score;
            if (score > bestScore) {
                bestScore = score;
                mccTree = tree;
            }
        }
        assertNotNull("an MCC tree should be selected", mccTree);

        // Capture the MCC tree's identity BEFORE annotation adds relationship metadata:
        //  - its posterior sample label (e.g. "STATE_500000"), and
        //  - its oriented, node-numbered Newick (node numbers are 1-based and map to the
        //    Translate block of the trees file).
        String mccState = mccTree.getID();
        String mccNewick = ((SRNode) mccTree.getRoot()).toShortNewickForLog(false);

        system.annotateMCCTree(mccTree, true);

        // --- report the constructed relationship system ------------------------------
        System.out.println("=== MCC relationship system (10 trees, " + BURNIN_PERCENTAGE
                + "% burn-in, heights kept) ===");
        System.out.println("MCC tree sample:               " + mccState);
        System.out.println("MCC tree (oriented, node-numbered) Newick:");
        System.out.println("  " + mccNewick);
        System.out.println("Highest log relationship credibility: " + bestScore);
        System.out.println("Log relationship credibilities of the " + totalTreesUsed + " analysed trees:");
        for (int i = 0; i < totalTreesUsed; i++) {
            System.out.println("  [" + i + "] " + analyzedTrees.get(i).getID() + " : " + logCredibilities[i]);
        }
        System.out.println("Ancestry relationships:    " + system.getAncestryMap().size());
        System.out.println("Orientation relationships: " + system.getOrientationMap().size());
        System.out.println();
        System.out.println(system.getSummary());

        // --- sanity checks on the overall shape --------------------------------------
        assertTrue("the relationship system should not be empty",
                system.getAncestryMap().size() + system.getOrientationMap().size() > 0);
        // every probability is count/10, so it must lie in (0, 1]
        for (AncestryRelationship r : system.getAncestryMap().values()) {
            assertTrue(r.getProbability() > 0.0 && r.getProbability() <= 1.0);
            assertEquals(r.getCount() / (double) totalTreesUsed, r.getProbability(), 1e-12);
        }
        for (OrientationRelationship r : system.getOrientationMap().values()) {
            assertTrue(r.getProbability() > 0.0 && r.getProbability() <= 1.0);
            assertEquals(r.getCount() / (double) totalTreesUsed, r.getProbability(), 1e-12);
        }
        assertTrue("MCC tree should have a finite credibility score", Double.isFinite(bestScore));

        // =============================================================================
        //  TRUE VALUES - FILL IN BY HAND
        // =============================================================================

        // -----------------------------------------------------------------------------
        // (1) Size of the relationship system.
        //     Set to a non-negative number to activate the check; -1 = not yet filled in.
        // -----------------------------------------------------------------------------
        int expectedAncestryCount = 12;
        int expectedOrientationCount = 46;

        if (expectedAncestryCount >= 0) {
            assertEquals("number of distinct ancestry relationships",
                    expectedAncestryCount, system.getAncestryMap().size());
        }
        if (expectedOrientationCount >= 0) {
            assertEquals("number of distinct orientation relationships",
                    expectedOrientationCount, system.getOrientationMap().size());
        }

        // -----------------------------------------------------------------------------
        // (2) The MCC tree: which posterior sample was selected (its STATE label) and its
        //     oriented, node-numbered topology (node numbers map to the Translate block).
        //     Leave the string empty to skip the check.
        // -----------------------------------------------------------------------------
        String expectedMccState = "_800000";
        if (!expectedMccState.isEmpty()) {
            assertEquals("MCC tree should be the expected posterior sample",
                    expectedMccState, mccState);
        }

        // -----------------------------------------------------------------------------
        // (3) Log relationship credibility of the MCC tree.
        //     Leave as Double.NaN to skip the check.
        // -----------------------------------------------------------------------------
        double expectedBestScore = -12.24689;
        if (!Double.isNaN(expectedBestScore)) {
            assertEquals("log relationship credibility of the MCC tree",
                    expectedBestScore, bestScore, 1e-5);
        }

        // -----------------------------------------------------------------------------
        // (4) Log relationship credibility of every analysed tree, in sample order.
        //     Leave the array empty to skip the check; it is only applied when its length
        //     matches the number of analysed trees.
        // -----------------------------------------------------------------------------
        double[] expectedLogCredibilities = {
                -14.7318, -15.42495, -14.32634, -14.7318, -14.03865, -15.42495, -13.63319, -12.24689, -14.03865, -14.03865
        };
        if (expectedLogCredibilities.length == totalTreesUsed) {
            System.out.println("asserting the crediblities:");
            for (int i = 0; i < totalTreesUsed; i++) {
                assertEquals("log relationship credibility for tree " + i
                                + " (" + analyzedTrees.get(i).getID() + ")",
                        expectedLogCredibilities[i], logCredibilities[i], 1e-5);
            }
        }

        // -----------------------------------------------------------------------------
        // (5) Individual ancestry relationships (A, T): the first occurrence (or the only
        //     occurrence) of taxon A is a direct ancestor of the MRCA of T.
        //     Taxon names are BASE names, i.e. without the "_first"/"_last" suffix
        //     (e.g. "Spheniscus_urbinai", not "Spheniscus_urbinai_first").
        //     Add one line per relationship you have counted by hand; the expected count is
        //     out of 10 analysed trees.
        // -----------------------------------------------------------------------------
        //
        //   assertAncestry(system, "Pygoscelis_grandis", set("Aptenodytes_forsteri"), <count>, totalTreesUsed);
        //   assertAncestry(system, "Spheniscus_urbinai",
        //           set("Spheniscus_demersus", "Spheniscus_megaramphus"), <count>, totalTreesUsed);
        //
        // TODO: fill in the true ancestry relationships here.

        // -----------------------------------------------------------------------------
        // (6) Individual orientation relationships (T1 => T2): T1 descends from the
        //     ancestral (left) lineage of the split, T2 from the descendant (right) lineage.
        //     ORDER MATTERS - (T1, T2) and (T2, T1) are different relationships.
        //     Add one line per relationship you have counted by hand.
        // -----------------------------------------------------------------------------
        //
        //   assertOrientation(system, set("Eudyptes_filholi"), set("Spheniscus_demersus"),
        //           <count>, totalTreesUsed);
        //   assertOrientation(system, set("Marplesornis_novaezealandiae"),
        //           set("Aptenodytes_forsteri", "Pygoscelis_grandis"), <count>, totalTreesUsed);
        //
        // TODO: fill in the true orientation relationships here.

        // -----------------------------------------------------------------------------
        // (7) Relationships that must NOT be present. In particular, a LAST occurrence of a
        //     multi-occurrence range never produces an ancestry relationship (only the first
        //     occurrence, or a singleton, does), and a split whose ancestral side contains
        //     only the last occurrence of the enclosing range produces no orientation
        //     relationship.
        // -----------------------------------------------------------------------------
        //
        //   assertNoAncestry(system, "<taxon>", set("..."));
        //   assertNoOrientation(system, set("..."), set("..."));
        //
        // TODO: fill in the relationships that must be absent here.
    }

    // ------------------------------------------------------------------
    //  Assertion helpers for hand-computed relationships
    // ------------------------------------------------------------------

    /** Asserts that ancestry relationship (A, T) was seen in exactly {@code count} of the trees. */
    private void assertAncestry(RelationshipSystem system, String ancestorTaxon,
                                Set<String> descendantTaxa, int count, int totalTrees) {
        AncestryRelationship key = new AncestryRelationship(ancestorTaxon, descendantTaxa);
        AncestryRelationship actual = system.getAncestryMap().get(key);
        assertNotNull("ancestry relationship " + key + " should exist", actual);
        assertEquals("count of ancestry relationship " + key, count, actual.getCount());
        assertEquals("probability of ancestry relationship " + key,
                count / (double) totalTrees, actual.getProbability(), 1e-10);
    }

    /** Asserts that orientation relationship (T1 => T2) was seen in exactly {@code count} trees. */
    private void assertOrientation(RelationshipSystem system, Set<String> ancestralTaxa,
                                   Set<String> descendantTaxa, int count, int totalTrees) {
        OrientationRelationship key = new OrientationRelationship(ancestralTaxa, descendantTaxa);
        OrientationRelationship actual = system.getOrientationMap().get(key);
        assertNotNull("orientation relationship " + key + " should exist", actual);
        assertEquals("count of orientation relationship " + key, count, actual.getCount());
        assertEquals("probability of orientation relationship " + key,
                count / (double) totalTrees, actual.getProbability(), 1e-10);
    }

    /** Asserts that ancestry relationship (A, T) does not occur in any of the analysed trees. */
    private void assertNoAncestry(RelationshipSystem system, String ancestorTaxon,
                                  Set<String> descendantTaxa) {
        AncestryRelationship key = new AncestryRelationship(ancestorTaxon, descendantTaxa);
        assertTrue("ancestry relationship " + key + " must NOT exist",
                !system.getAncestryMap().containsKey(key));
    }

    /** Asserts that orientation relationship (T1 => T2) does not occur in any analysed tree. */
    private void assertNoOrientation(RelationshipSystem system, Set<String> ancestralTaxa,
                                     Set<String> descendantTaxa) {
        OrientationRelationship key = new OrientationRelationship(ancestralTaxa, descendantTaxa);
        assertTrue("orientation relationship " + key + " must NOT exist",
                !system.getOrientationMap().containsKey(key));
    }

    /** Convenience set builder, as in {@link WithinRangeCladeTest}. */
    private Set<String> set(String... taxa) {
        return new TreeSet<>(Arrays.asList(taxa));
    }

    // ------------------------------------------------------------------
    //  Helpers (identical to MccCladeSystemTest)
    // ------------------------------------------------------------------

    /**
     * Reads SR trees from a NEXUS file using the same conversion path as
     * {@link sr.treeannotator.SRTreeAnnotator#readTrees()}: parse with {@link NexusParser},
     * then wrap each parsed tree in an {@link SRTree} and orientate it.
     */
    private List<SRTree> readTrees(File file) throws Exception {
        List<SRTree> trees = new ArrayList<>();
        NexusParser parser = new NexusParser();
        parser.parseFile(file);
        for (Tree tree : parser.trees) {
            if (tree instanceof SRTree) {
                trees.add((SRTree) tree);
            } else {
                // NOTE: unlike SRTreeAnnotator.readTrees(), we declare the multi-occurrence
                // ranges explicitly. SRTree.initSRanges() can otherwise only reconstruct ranges
                // by splitting taxon IDs on the last '_', which misclassifies binomial names
                // (e.g. "Spheniscus_demersus" -> range "Spheniscus", occurrence "demersus") and
                // throws "taxa with last occurrence only". Declaring the first/last pairs takes
                // the pre-set-ranges branch, where every remaining tip auto-becomes a singleton
                // range - matching how the source XML defines the ranges.
                SRTree srTree = new SRTree();
                srTree.setInputValue("stratigraphicRange", multiOccurrenceRanges());
                srTree.assignFrom(tree);
                srTree.orientateTree();
                trees.add(srTree);
            }
        }
        return trees;
    }

    /**
     * The multi-occurrence stratigraphic ranges (first/last occurrence pairs) declared in
     * examples/morph+dna_co_ncc_for_mcc_test.xml. A fresh set is built per tree because
     * {@link SRTree#assignFrom} mutates the range objects (node-number bookkeeping).
     */
    private List<StratigraphicRange> multiOccurrenceRanges() {
        List<StratigraphicRange> sranges = new ArrayList<>();
        sranges.add(range("Spheniscus_urbinai_first", "Spheniscus_urbinai_last"));
        sranges.add(range("Spheniscus_megaramphus_first", "Spheniscus_megaramphus_last"));
        sranges.add(range("Pygoscelis_grandis_first", "Pygoscelis_grandis_last"));
        return sranges;
    }

    private StratigraphicRange range(String firstOcc, String lastOcc) {
        StratigraphicRange sr = new StratigraphicRange();
        sr.setInputValue("firstOccurrence", new Taxon(firstOcc));
        sr.setInputValue("lastOccurrence", new Taxon(lastOcc));
        return sr;
    }

    /**
     * Resolves the trees file regardless of the working directory the test is launched from
     * (project root or module root).
     */
    private File locateTreesFile() {
        String[] candidates = {
                "test_data/" + TREES_FILE,
                TREES_FILE,
                "../stratigraphic-ranges/test_data/" + TREES_FILE,
        };
        for (String c : candidates) {
            File f = new File(c);
            if (f.exists()) {
                return f.getAbsoluteFile();
            }
        }
        throw new IllegalStateException("Could not find " + TREES_FILE
                + " (looked in: test_data/, ., ../stratigraphic-ranges/test_data/)");
    }
}
