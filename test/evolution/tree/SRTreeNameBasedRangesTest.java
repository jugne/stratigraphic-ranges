package evolution.tree;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.Test;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRTree;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

/**
 * Stratigraphic range detection from taxon names (no explicit stratigraphicRange input),
 * the path used when relogging trees read from file.
 */
public class SRTreeNameBasedRangesTest {

    private static SRTree fromNewick(String newick) {
        Tree initial = new TreeParser(newick, false);
        SRTree tree = new SRTree();
        tree.assignFrom(initial);
        return tree;
    }

    private static void assertThrowsContaining(String newick, String messagePart) {
        try {
            fromNewick(newick);
            fail("Expected RuntimeException containing: " + messagePart);
        } catch (RuntimeException e) {
            assertTrue("unexpected message: " + e.getMessage(), e.getMessage().contains(messagePart));
        }
    }

    @Test
    public void unsuffixedTaxonIsSingleOccurrence() {
        SRTree tree = fromNewick("((A:1.0,B_last:0.5):1.0,B_first:0.0):0.5");
        assertEquals(2, tree.getSRanges().size());
        StratigraphicRange a = tree.sRangesContainsID("A");
        assertNotNull(a);
        assertTrue(a.isSingleFossilRange());
        StratigraphicRange b = tree.sRangesContainsID("B_first");
        assertNotNull(b);
        assertFalse(b.isSingleFossilRange());
        assertEquals("B_last", b.getLastOccurrenceID());
    }

    @Test
    public void underscoreInSingleOccurrenceNameIsKept() {
        SRTree tree = fromNewick("((Aptenodytes_forsteri:1.0,B_last:0.5):1.0,B_first:0.0):0.5");
        StratigraphicRange a = tree.sRangesContainsID("Aptenodytes_forsteri");
        assertNotNull(a);
        assertTrue(a.isSingleFossilRange());
        assertEquals("Aptenodytes_forsteri", a.getID());
    }

    @Test
    public void firstOnlyTaxonIsSingleOccurrence() {
        SRTree tree = fromNewick("((A_first:1.0,B_last:0.5):1.0,B_first:0.0):0.5");
        StratigraphicRange a = tree.sRangesContainsID("A_first");
        assertNotNull(a);
        assertTrue(a.isSingleFossilRange());
    }

    @Test
    public void lastOnlyTaxonIsRejected() {
        assertThrowsContaining("((A:1.0,B_last:0.5):1.0,C:0.2):0.5", "last occurrence only: [B]");
    }

    /** Builds the SRTree the way the XML parser does: through the "tree" input and initAndValidate. */
    private static SRTree fromNewickAsStartingTree(String newick) {
        Tree initial = new TreeParser(newick, false);
        SRTree tree = new SRTree();
        tree.initByName("tree", initial);
        return tree;
    }

    @Test
    public void nonSampledAncestorFirstOccurrenceIsRejectedForStartingTree() {
        // B_first is a regular leaf at height 0.7, not a sampled ancestor
        String[] newicks = {
                "((A:1.0,B_last:0.5):1.0,B_first:0.7):0.5",
                // same, with _first encountered before _last
                "((A:1.0,B_first:0.7):1.0,B_last:0.5):0.5"};
        for (String newick : newicks) {
            try {
                fromNewickAsStartingTree(newick);
                fail("Expected RuntimeException for " + newick);
            } catch (RuntimeException e) {
                assertTrue("unexpected message: " + e.getMessage(),
                        e.getMessage().contains("B_first is not a sampled ancestor"));
            }
        }
    }

    @Test
    public void nonSampledAncestorFirstOccurrenceIsToleratedOnAssignFrom() {
        // feast's TreeLogFileState assigns a placeholder ultrametric tree before the logged trees;
        // ranges are still built by name, only the sampled-ancestor validation is skipped
        SRTree tree = fromNewick("((A:1.0,B_last:0.5):1.0,B_first:0.7):0.5");
        assertEquals(2, tree.getSRanges().size());
        StratigraphicRange b = tree.sRangesContainsID("B_first");
        assertNotNull(b);
        assertFalse(b.isSingleFossilRange());
        assertEquals(2, b.getNodeNrs().size());
    }

    @Test
    public void duplicateRangeNamesAreRejected() {
        assertThrowsContaining("((B:1.0,B_last:0.5):1.0,B_first:0.0):0.5", "share the name B");
    }
}