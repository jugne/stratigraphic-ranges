package evolution.tree;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.Test;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.BranchRateLogger;
import sr.evolution.tree.SRTree;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class BranchRateLoggerTest {

    private static SRTree threeRangeTree() {
        String newick = "(((((A:3.4,2_last:0.0):1.0,2_first:0.0):0.7,(B:3.5,(3_last:1.7,3_first:0.0):0.8):1.6):0.55,1_last:0.0):0.85,1_first:0.0):0.5";
        Tree initial = new TreeParser(newick, false);
        ArrayList<StratigraphicRange> ranges = new ArrayList<>();
        for (String name : new String[]{"1", "2", "3"}) {
            StratigraphicRange r = new StratigraphicRange();
            r.setInputValue("firstOccurrence", new Taxon(name + "_first"));
            r.setInputValue("lastOccurrence", new Taxon(name + "_last"));
            ranges.add(r);
        }
        SRTree tree = new SRTree();
        tree.setInputValue("stratigraphicRange", ranges);
        tree.assignFrom(initial);
        return tree;
    }

    /** @return map childLabel -> [isRange, rangeName] */
    private static Map<String, String[]> logRows(SRTree tree) {
        BranchRateLogger logger = new BranchRateLogger();
        logger.initByName("tree", tree);
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream out = new PrintStream(bytes);
        logger.init(out);
        logger.log(0, out);
        out.flush();
        Map<String, String[]> rows = new HashMap<>();
        String[] lines = bytes.toString().split("\\r?\\n");
        for (int i = 1; i < lines.length; i++) {
            String[] cols = lines[i].split("\t", -1);
            rows.put(cols[0], new String[]{cols[3], cols[4]});
        }
        return rows;
    }

    @Test
    public void firstOccurrenceBranchIsNotPartOfRange() {
        Map<String, String[]> rows = logRows(threeRangeTree());
        // 1_first is the sampled ancestor at the root, so it has no branch and no row
        assertTrue(!rows.containsKey("1_first"));
        for (String first : new String[]{"2_first", "3_first"}) {
            assertTrue("missing row for " + first, rows.containsKey(first));
            assertEquals("isRange for " + first, "0", rows.get(first)[0]);
            assertEquals("rangeName for " + first, "", rows.get(first)[1]);
        }
    }

    @Test
    public void lastOccurrenceBranchIsPartOfRange() {
        Map<String, String[]> rows = logRows(threeRangeTree());
        for (String name : new String[]{"1", "2", "3"}) {
            String last = name + "_last";
            assertTrue("missing row for " + last, rows.containsKey(last));
            assertEquals("isRange for " + last, "1", rows.get(last)[0]);
            assertEquals("rangeName for " + last, name, rows.get(last)[1]);
        }
    }

    @Test
    public void singleOccurrenceTaxaAreNotRanges() {
        Map<String, String[]> rows = logRows(threeRangeTree());
        for (String tip : new String[]{"A", "B"}) {
            assertEquals("isRange for " + tip, "0", rows.get(tip)[0]);
            assertEquals("rangeName for " + tip, "", rows.get(tip)[1]);
        }
    }
}