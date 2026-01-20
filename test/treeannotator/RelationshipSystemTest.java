package treeannotator;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import junit.framework.TestCase;
import org.junit.Test;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRTree;
import sr.treeannotator.AncestryRelationship;
import sr.treeannotator.OrientationRelationship;
import sr.treeannotator.RelationshipSystem;

import java.util.*;

/**
 * Test for RelationshipSystem that collects ancestry and orientation relationships
 * from SR trees for tree annotation.
 *
 * This test uses simple scenarios that can be verified by hand:
 *
 * TEST 1: Singleton taxa (A, B, C) - 10 trees
 * ============================================
 *
 * Trees 1-5: A is sampled ancestor of {B, C}, B ancestral to C
 *   Newick: (A:0.0,(B:1.0,C:1.0):1.0):1.0
 *   - Ancestry: (A, {B,C})
 *   - Orientation: ({B}, {C})
 *
 * Trees 6-7: A is sampled ancestor of {B, C}, C ancestral to B
 *   Newick: (A:0.0,(C:1.0,B:1.0):1.0):1.0
 *   - Ancestry: (A, {B,C})
 *   - Orientation: ({C}, {B})
 *
 * Trees 8-9: A is regular tip, ((A,B),C), A ancestral to B
 *   Newick: ((A:1.0,B:1.0):1.0,C:1.0):1.0
 *   - Orientation: ({A,B}, {C}) at root
 *   - Orientation: ({A}, {B}) at A-B node
 *
 * Tree 10: A is regular tip, ((B,A),C), B ancestral to A
 *   Newick: ((B:1.0,A:1.0):1.0,C:1.0):1.0
 *   - Orientation: ({A,B}, {C}) at root
 *   - Orientation: ({B}, {A}) at A-B node
 *
 * EXPECTED POSTERIORS:
 *   Ancestry (A, {B,C}): 7/10 = 0.70
 *   Orientation ({B}, {C}): 5/10 = 0.50
 *   Orientation ({C}, {B}): 2/10 = 0.20
 *   Orientation ({A,B}, {C}): 3/10 = 0.30
 *   Orientation ({A}, {B}): 2/10 = 0.20
 *   Orientation ({B}, {A}): 1/10 = 0.10
 *
 *
 * TEST 2: Stratigraphic range (F_first, F_last) + singletons (G, H) - 5 trees
 * ===========================================================================
 *
 * Trees 1-3: F_first is SA, F_last is SA of {G, H}, G ancestral to H
 *   - Ancestry: (F, {F,G,H}) from F_first
 *   - NO ancestry from F_last (it's a last occurrence!)
 *   - Orientation: ({G}, {H})
 *
 * Trees 4-5: F_first is SA, F_last is tip (not SA), ((F_last,G),H)
 *   - Ancestry: (F, {F,G,H}) from F_first
 *   - Orientation: ({F,G}, {H}) at root
 *   - Orientation: ({F}, {G}) at F_last-G node
 *
 * EXPECTED POSTERIORS:
 *   Ancestry (F, {F,G,H}): 5/5 = 1.00
 *   Orientation ({G}, {H}): 3/5 = 0.60
 *   Orientation ({F,G}, {H}): 2/5 = 0.40
 *   Orientation ({F}, {G}): 2/5 = 0.40
 *
 */
public class RelationshipSystemTest extends TestCase {

    private static final double TOLERANCE = 1e-10;

    /**
     * Test 1: Singleton taxa only - verifies basic ancestry and orientation collection.
     */
    @Test
    public void testSingletonTaxa() throws Exception {
        RelationshipSystem system = new RelationshipSystem();

        // Trees 1-5: A is SA of {B,C}, B ancestral to C
        // Newick: (A:0.0,(B:1.0,C:1.0):1.0):1.0
        String newick_A_SA_B_anc_C = "(A:0.0,(B:1.0,C:1.0):1.0):1.0";
        for (int i = 0; i < 5; i++) {
            SRTree tree = createSingletonTree(newick_A_SA_B_anc_C);
            system.add(tree);
        }

        // Trees 6-7: A is SA of {B,C}, C ancestral to B
        // Newick: (A:0.0,(C:1.0,B:1.0):1.0):1.0
        String newick_A_SA_C_anc_B = "(A:0.0,(C:1.0,B:1.0):1.0):1.0";
        for (int i = 0; i < 2; i++) {
            SRTree tree = createSingletonTree(newick_A_SA_C_anc_B);
            system.add(tree);
        }

        // Trees 8-9: A is tip, ((A,B),C), A ancestral to B
        String newick_A_tip_A_anc_B = "((A:1.0,B:1.0):1.0,C:1.0):1.0";
        for (int i = 0; i < 2; i++) {
            SRTree tree = createSingletonTree(newick_A_tip_A_anc_B);
            system.add(tree);
        }

        // Tree 10: A is tip, ((B,A),C), B ancestral to A
        String newick_A_tip_B_anc_A = "((B:1.0,A:1.0):1.0,C:1.0):1.0";
        SRTree tree10 = createSingletonTree(newick_A_tip_B_anc_A);
        system.add(tree10);

        // Calculate posteriors
        system.calculatePosteriorProbabilities(10);

        // Print summary for manual verification
        System.out.println("=== TEST 1: Singleton Taxa ===");
        System.out.println(system.getSummary());

        // Verify ancestry relationships
        Map<AncestryRelationship, AncestryRelationship> ancestryMap = system.getAncestryMap();

        // (A, {B,C}) should have probability 0.7
        AncestryRelationship ancA_BC = new AncestryRelationship("A", new TreeSet<>(Arrays.asList("B", "C")));
        AncestryRelationship foundAnc = ancestryMap.get(ancA_BC);
        assertNotNull("Ancestry (A, {B,C}) should exist", foundAnc);
        assertEquals("Ancestry (A, {B,C}) count", 7, foundAnc.getCount());
        assertEquals("Ancestry (A, {B,C}) probability", 0.7, foundAnc.getProbability(), TOLERANCE);

        // Verify orientation relationships
        Map<OrientationRelationship, OrientationRelationship> orientMap = system.getOrientationMap();

        // ({B}, {C}) should have probability 0.5
        OrientationRelationship orient_B_C = new OrientationRelationship(
            new TreeSet<>(Arrays.asList("B")),
            new TreeSet<>(Arrays.asList("C"))
        );
        OrientationRelationship foundOrient_B_C = orientMap.get(orient_B_C);
        assertNotNull("Orientation ({B}, {C}) should exist", foundOrient_B_C);
        assertEquals("Orientation ({B}, {C}) count", 5, foundOrient_B_C.getCount());
        assertEquals("Orientation ({B}, {C}) probability", 0.5, foundOrient_B_C.getProbability(), TOLERANCE);

        // ({C}, {B}) should have probability 0.2
        OrientationRelationship orient_C_B = new OrientationRelationship(
            new TreeSet<>(Arrays.asList("C")),
            new TreeSet<>(Arrays.asList("B"))
        );
        OrientationRelationship foundOrient_C_B = orientMap.get(orient_C_B);
        assertNotNull("Orientation ({C}, {B}) should exist", foundOrient_C_B);
        assertEquals("Orientation ({C}, {B}) count", 2, foundOrient_C_B.getCount());
        assertEquals("Orientation ({C}, {B}) probability", 0.2, foundOrient_C_B.getProbability(), TOLERANCE);

        // ({A,B}, {C}) should have probability 0.3
        OrientationRelationship orient_AB_C = new OrientationRelationship(
            new TreeSet<>(Arrays.asList("A", "B")),
            new TreeSet<>(Arrays.asList("C"))
        );
        OrientationRelationship foundOrient_AB_C = orientMap.get(orient_AB_C);
        assertNotNull("Orientation ({A,B}, {C}) should exist", foundOrient_AB_C);
        assertEquals("Orientation ({A,B}, {C}) count", 3, foundOrient_AB_C.getCount());
        assertEquals("Orientation ({A,B}, {C}) probability", 0.3, foundOrient_AB_C.getProbability(), TOLERANCE);

        // ({A}, {B}) should have probability 0.2
        OrientationRelationship orient_A_B = new OrientationRelationship(
            new TreeSet<>(Arrays.asList("A")),
            new TreeSet<>(Arrays.asList("B"))
        );
        OrientationRelationship foundOrient_A_B = orientMap.get(orient_A_B);
        assertNotNull("Orientation ({A}, {B}) should exist", foundOrient_A_B);
        assertEquals("Orientation ({A}, {B}) count", 2, foundOrient_A_B.getCount());
        assertEquals("Orientation ({A}, {B}) probability", 0.2, foundOrient_A_B.getProbability(), TOLERANCE);

        // ({B}, {A}) should have probability 0.1
        OrientationRelationship orient_B_A = new OrientationRelationship(
            new TreeSet<>(Arrays.asList("B")),
            new TreeSet<>(Arrays.asList("A"))
        );
        OrientationRelationship foundOrient_B_A = orientMap.get(orient_B_A);
        assertNotNull("Orientation ({B}, {A}) should exist", foundOrient_B_A);
        assertEquals("Orientation ({B}, {A}) count", 1, foundOrient_B_A.getCount());
        assertEquals("Orientation ({B}, {A}) probability", 0.1, foundOrient_B_A.getProbability(), TOLERANCE);

        System.out.println("TEST 1 PASSED: All posteriors match expected values!\n");
    }

    /**
     * Test 2: Stratigraphic range - verifies that only FIRST occurrence creates ancestry.
     */
    @Test
    public void testStratigraphicRange() throws Exception {
        RelationshipSystem system = new RelationshipSystem();

        // Trees 1-3: F_first is SA at root, F_last is SA of {G,H}, G ancestral to H
        // Structure: (F_first:0.0, (F_last:0.0, (G:1.0, H:1.0):1.0):1.0):1.0
        // F_first is SA of everything, F_last is SA of {G,H}
        String newick_Flast_SA = "(F_first:0.0,(F_last:0.0,(G:1.0,H:1.0):1.0):1.0):1.0";
        for (int i = 0; i < 3; i++) {
            SRTree tree = createRangeTree(newick_Flast_SA, "F_first", "F_last");
            system.add(tree);
        }

        // Trees 4-5: F_first is SA, F_last is regular tip in ((F_last,G),H)
        // Structure: (F_first:0.0, ((F_last:1.0, G:1.0):1.0, H:1.0):1.0):1.0
        String newick_Flast_tip = "(F_first:0.0,((F_last:1.0,G:1.0):1.0,H:1.0):1.0):1.0";
        for (int i = 0; i < 2; i++) {
            SRTree tree = createRangeTree(newick_Flast_tip, "F_first", "F_last");
            system.add(tree);
        }

        // Calculate posteriors
        system.calculatePosteriorProbabilities(5);

        // Print summary for manual verification
        System.out.println("=== TEST 2: Stratigraphic Range ===");
        System.out.println(system.getSummary());

        // Verify ancestry relationships
        Map<AncestryRelationship, AncestryRelationship> ancestryMap = system.getAncestryMap();

        // (F, {F,G,H}) from F_first should have probability 1.0 (all 5 trees)
        AncestryRelationship ancF_FGH = new AncestryRelationship("F", new TreeSet<>(Arrays.asList("F", "G", "H")));
        AncestryRelationship foundAncF = ancestryMap.get(ancF_FGH);
        assertNotNull("Ancestry (F, {F,G,H}) should exist", foundAncF);
        assertEquals("Ancestry (F, {F,G,H}) count", 5, foundAncF.getCount());
        assertEquals("Ancestry (F, {F,G,H}) probability", 1.0, foundAncF.getProbability(), TOLERANCE);

        // IMPORTANT: There should be NO ancestry relationship from F_last!
        // F_last is a last occurrence, so even when it's a SA, it should NOT create ancestry
        // If the bug existed, we would see (F, {G,H}) with count 3
        AncestryRelationship ancF_GH = new AncestryRelationship("F", new TreeSet<>(Arrays.asList("G", "H")));
        AncestryRelationship foundAncF_GH = ancestryMap.get(ancF_GH);
        assertNull("Ancestry (F, {G,H}) should NOT exist (F_last is last occurrence)", foundAncF_GH);

        // Verify orientation relationships
        Map<OrientationRelationship, OrientationRelationship> orientMap = system.getOrientationMap();

        // ({G}, {H}) should have probability 0.6 (trees 1-3)
        OrientationRelationship orient_G_H = new OrientationRelationship(
            new TreeSet<>(Arrays.asList("G")),
            new TreeSet<>(Arrays.asList("H"))
        );
        OrientationRelationship foundOrient_G_H = orientMap.get(orient_G_H);
        assertNotNull("Orientation ({G}, {H}) should exist", foundOrient_G_H);
        assertEquals("Orientation ({G}, {H}) count", 3, foundOrient_G_H.getCount());
        assertEquals("Orientation ({G}, {H}) probability", 0.6, foundOrient_G_H.getProbability(), TOLERANCE);

        // ({F,G}, {H}) should have probability 0.4 (trees 4-5)
        OrientationRelationship orient_FG_H = new OrientationRelationship(
            new TreeSet<>(Arrays.asList("F", "G")),
            new TreeSet<>(Arrays.asList("H"))
        );
        OrientationRelationship foundOrient_FG_H = orientMap.get(orient_FG_H);
        assertNotNull("Orientation ({F,G}, {H}) should exist", foundOrient_FG_H);
        assertEquals("Orientation ({F,G}, {H}) count", 2, foundOrient_FG_H.getCount());
        assertEquals("Orientation ({F,G}, {H}) probability", 0.4, foundOrient_FG_H.getProbability(), TOLERANCE);

        // ({F}, {G}) should have probability 0.4 (trees 4-5, from F_last-G bifurcation)
        OrientationRelationship orient_F_G = new OrientationRelationship(
            new TreeSet<>(Arrays.asList("F")),
            new TreeSet<>(Arrays.asList("G"))
        );
        OrientationRelationship foundOrient_F_G = orientMap.get(orient_F_G);
        assertNotNull("Orientation ({F}, {G}) should exist", foundOrient_F_G);
        assertEquals("Orientation ({F}, {G}) count", 2, foundOrient_F_G.getCount());
        assertEquals("Orientation ({F}, {G}) probability", 0.4, foundOrient_F_G.getProbability(), TOLERANCE);

        System.out.println("TEST 2 PASSED: First occurrence creates ancestry, last occurrence does NOT!\n");
    }

    /**
     * Helper: Create an SRTree from Newick with singleton taxa (no ranges).
     * Each singleton taxon (A, B, C) needs to be registered as a single-fossil range.
     * For single-fossil ranges, firstOccurrence and lastOccurrence must be the SAME taxon.
     */
    private SRTree createSingletonTree(String newick) throws Exception {
        Tree treeInitial = new TreeParser(newick, false);

        // Create single-fossil ranges for each singleton taxon
        // For single-fossil range: firstOccurrence == lastOccurrence (same taxon)
        ArrayList<StratigraphicRange> sranges = new ArrayList<>();

        StratigraphicRange srA = new StratigraphicRange();
        Taxon taxonA = new Taxon("A");
        srA.setInputValue("firstOccurrence", taxonA);
        srA.setInputValue("lastOccurrence", taxonA);
        sranges.add(srA);

        StratigraphicRange srB = new StratigraphicRange();
        Taxon taxonB = new Taxon("B");
        srB.setInputValue("firstOccurrence", taxonB);
        srB.setInputValue("lastOccurrence", taxonB);
        sranges.add(srB);

        StratigraphicRange srC = new StratigraphicRange();
        Taxon taxonC = new Taxon("C");
        srC.setInputValue("firstOccurrence", taxonC);
        srC.setInputValue("lastOccurrence", taxonC);
        sranges.add(srC);

        SRTree tree = new SRTree();
        tree.setInputValue("stratigraphicRange", sranges);
        tree.assignFrom(treeInitial);
        tree.initAndValidate();
        return tree;
    }

    /**
     * Helper: Create an SRTree from Newick with one stratigraphic range (F) and singletons (G, H).
     */
    private SRTree createRangeTree(String newick, String firstOccurrence, String lastOccurrence) throws Exception {
        Tree treeInitial = new TreeParser(newick, false);

        ArrayList<StratigraphicRange> sranges = new ArrayList<>();

        // Create the stratigraphic range for F
        StratigraphicRange srF = new StratigraphicRange();
        Taxon taxonFirst = new Taxon(firstOccurrence);
        Taxon taxonLast = new Taxon(lastOccurrence);
        srF.setInputValue("firstOccurrence", taxonFirst);
        srF.setInputValue("lastOccurrence", taxonLast);
        sranges.add(srF);

        // Create single-fossil ranges for G and H (first == last for singletons)
        StratigraphicRange srG = new StratigraphicRange();
        Taxon taxonG = new Taxon("G");
        srG.setInputValue("firstOccurrence", taxonG);
        srG.setInputValue("lastOccurrence", taxonG);
        sranges.add(srG);

        StratigraphicRange srH = new StratigraphicRange();
        Taxon taxonH = new Taxon("H");
        srH.setInputValue("firstOccurrence", taxonH);
        srH.setInputValue("lastOccurrence", taxonH);
        sranges.add(srH);

        SRTree tree = new SRTree();
        tree.setInputValue("stratigraphicRange", sranges);
        tree.assignFrom(treeInitial);
        tree.initAndValidate();

        return tree;
    }
}
