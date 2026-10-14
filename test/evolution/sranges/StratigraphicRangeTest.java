package evolution.sranges;

import beast.base.evolution.alignment.Taxon;
import org.junit.Test;
import sr.evolution.sranges.StratigraphicRange;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

public class StratigraphicRangeTest {

    @Test
    public void onlyFirstOccurrenceGivesClearError() {
        StratigraphicRange range = new StratigraphicRange();
        range.taxonFirstOccurrenceInput.setValue(new Taxon("A_first"), range);
        try {
            range.initAndValidate();
            fail("Expected a RuntimeException when only firstOccurrence is given");
        } catch (RuntimeException e) {
            assertTrue(e.getMessage().contains("Either firstOccurrence or lastOccurence is not specified"));
        }
    }

    @Test
    public void onlyLastOccurrenceGivesClearError() {
        StratigraphicRange range = new StratigraphicRange();
        range.taxonLastOccurrenceInput.setValue(new Taxon("A_last"), range);
        try {
            range.initAndValidate();
            fail("Expected a RuntimeException when only lastOccurrence is given");
        } catch (RuntimeException e) {
            assertTrue(e.getMessage().contains("Either firstOccurrence or lastOccurence is not specified"));
        }
    }

    @Test
    public void bothOccurrencesInitialise() {
        StratigraphicRange range = new StratigraphicRange();
        range.taxonFirstOccurrenceInput.setValue(new Taxon("A_first"), range);
        range.taxonLastOccurrenceInput.setValue(new Taxon("A_last"), range);
        range.initAndValidate();
        assertEquals("A_first", range.getFirstOccurrenceID());
        assertEquals("A_last", range.getLastOccurrenceID());
        assertFalse(range.isSingleFossilRange());
    }

    @Test
    public void sameTaxonTwiceIsSingleFossilRange() {
        StratigraphicRange range = new StratigraphicRange();
        Taxon a = new Taxon("A");
        range.taxonFirstOccurrenceInput.setValue(a, range);
        range.taxonLastOccurrenceInput.setValue(a, range);
        range.initAndValidate();
        assertTrue(range.isSingleFossilRange());
    }
}