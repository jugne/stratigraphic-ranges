package evolution.tree;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.Test;
import sr.evolution.tree.SRTree;

import static org.junit.Assert.assertTrue;

public class SRTreeStateXMLTest {

    @Test
    public void stateXMLCarriesOrientationMetadata() {
        Tree initial = new TreeParser("((A:1.0,B_last:0.5):1.0,B_first:0.0):0.5", false);
        SRTree tree = new SRTree();
        tree.assignFrom(initial);
        tree.setID("tree");
        String xml = tree.toXML();
        assertTrue(xml, xml.startsWith("<statenode id='tree'>"));
        assertTrue(xml, xml.contains("orientation=ancestor"));
        assertTrue(xml, xml.contains("orientation=descendant"));
    }
}