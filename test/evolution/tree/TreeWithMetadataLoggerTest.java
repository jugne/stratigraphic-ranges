package evolution.tree;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.Test;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRTree;
import sr.evolution.tree.TreeWithMetadataLogger;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class TreeWithMetadataLoggerTest {

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

    private static String logOnce(boolean logRanges) {
        TreeWithMetadataLogger logger = new TreeWithMetadataLogger();
        logger.initByName("tree", threeRangeTree(), "logRanges", logRanges);
        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream out = new PrintStream(bytes);
        logger.log(0, out);
        out.flush();
        return bytes.toString();
    }

    @Test
    public void rangesLoggedByDefault() {
        String newick = logOnce(true);
        assertTrue(newick.contains("range=1"));
        assertTrue(newick.contains("orientation="));
    }

    @Test
    public void logRangesFalseSuppressesRangeMetadata() {
        String newick = logOnce(false);
        assertFalse("range metadata should be absent when logRanges=false: " + newick,
                newick.contains("range="));
        assertTrue(newick.contains("orientation="));
    }
}