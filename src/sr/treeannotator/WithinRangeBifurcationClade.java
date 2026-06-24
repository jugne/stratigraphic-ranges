package sr.treeannotator;

import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

/**
 * Clade type 2: a bifurcation event within a range.
 *
 * Denoted (A, T1, T2), where {A} ∪ T is monophyletic (T = T1 ∪ T2), the MRCA of {A} ∪ T is
 * within the range of A, the taxa T1 are descendants of the ancestral (left) lineage that
 * continues the range of A, and the taxa T2 are descendants of the descendant (right) lineage.
 *
 * Orientation matters: T1 is ancestral, T2 is descendant.
 *
 * @author Alexandra Gavryushkina
 */
public class WithinRangeBifurcationClade implements SRClade {

    private final String rangeTaxon;          // A - the range within which the bifurcation occurs
    private final Set<String> ancestralTaxa;  // T1 - ancestral (left) lineage
    private final Set<String> descendantTaxa; // T2 - descendant (right) lineage
    private int count;
    private double probability;

    public WithinRangeBifurcationClade(String rangeTaxon, Set<String> ancestralTaxa, Set<String> descendantTaxa) {
        this.rangeTaxon = rangeTaxon;
        this.ancestralTaxa = new TreeSet<>(ancestralTaxa);
        this.descendantTaxa = new TreeSet<>(descendantTaxa);
        this.count = 0;
        this.probability = 0.0;
    }

    public String getRangeTaxon() {
        return rangeTaxon;
    }

    public Set<String> getAncestralTaxa() {
        return ancestralTaxa;
    }

    public Set<String> getDescendantTaxa() {
        return descendantTaxa;
    }

    @Override
    public String getType() {
        return "within_range_bifurcation";
    }

    @Override
    public int getCount() {
        return count;
    }

    @Override
    public void setCount(int count) {
        this.count = count;
    }

    @Override
    public void incrementCount() {
        this.count++;
    }

    @Override
    public double getProbability() {
        return probability;
    }

    @Override
    public void setProbability(double probability) {
        this.probability = probability;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        WithinRangeBifurcationClade that = (WithinRangeBifurcationClade) o;
        return Objects.equals(rangeTaxon, that.rangeTaxon) &&
               Objects.equals(ancestralTaxa, that.ancestralTaxa) &&
               Objects.equals(descendantTaxa, that.descendantTaxa);
    }

    @Override
    public int hashCode() {
        return Objects.hash(rangeTaxon, ancestralTaxa, descendantTaxa);
    }

    @Override
    public String toString() {
        return "(" + rangeTaxon + ", " + ancestralTaxa + ", " + descendantTaxa + ")";
    }
}
