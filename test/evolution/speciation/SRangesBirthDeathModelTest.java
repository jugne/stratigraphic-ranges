package test.evolution.speciation;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import sr.evolution.tree.SRTree;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import junit.framework.TestCase;
import org.junit.jupiter.api.Test;
import sr.speciation.SRangesBirthDeathModel;
import sr.evolution.sranges.StratigraphicRange;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by gavryusa on 24/07/17.
 */
public class SRangesBirthDeathModelTest extends TestCase {

    @Test
    public void testLikelihood() throws Exception {

        ArrayList<String> taxa = new ArrayList<String>(Arrays.asList("1", "2", "3"));
        String newick = "(((((A:3.4,2_last:0.0):1.0,2_first:0.0):0.7,(B:3.5,(3_last:1.7,3_first:0.0):0.8):1.6):0.55,1_last:0.0):0.85,1_first:0.0):0.5";
        Tree tree_initial = new TreeParser(newick, false);

        StratigraphicRange sr1 = new StratigraphicRange();
        Taxon taxon1_first = new Taxon("1_first");
        Taxon taxon1_last = new Taxon("1_last");
        sr1.setInputValue("firstOccurrence", taxon1_first);
        sr1.setInputValue("lastOccurrence", taxon1_last);
        StratigraphicRange sr2 = new StratigraphicRange();
        Taxon taxon2_first = new Taxon("2_first");
        Taxon taxon2_last = new Taxon("2_last");
        sr2.setInputValue("firstOccurrence", taxon2_first);
        sr2.setInputValue("lastOccurrence", taxon2_last);
        StratigraphicRange sr3 = new StratigraphicRange();
        Taxon taxon3_first = new Taxon("3_first");
        Taxon taxon3_last = new Taxon("3_last");
        sr3.setInputValue("firstOccurrence", taxon3_first);
        sr3.setInputValue("lastOccurrence", taxon3_last);
        ArrayList<StratigraphicRange> sranges = new ArrayList<>();
        sranges.add(sr1);
        sranges.add(sr2);
        sranges.add(sr3);
        SRTree tree = new SRTree();
        tree.setInputValue("stratigraphicRange", sranges);
        tree.assignFrom(tree_initial);

        SRangesBirthDeathModel model = new SRangesBirthDeathModel();
        model.setInputValue("tree", tree);
        model.setInputValue("origin", new RealParameter("7.0"));
        model.setInputValue("birthRate", new RealParameter("1.5"));
        model.setInputValue("deathRate", new RealParameter("0.5"));
        model.setInputValue("samplingRate", new RealParameter("0.1"));
        model.setInputValue("removalProbability", new RealParameter("0.0"));
        model.setInputValue("rho", new RealParameter("0.5"));

        model.initAndValidate();

//        assertEquals(-33.57179092868063, model.calculateLogP(), 1e-14);
        assertEquals(-33.74668640318646, model.calculateLogP(), 1e-14);
    }
}
