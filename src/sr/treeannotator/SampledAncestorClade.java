package sr.treeannotator;

import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

/**
 * Clade type 3: an SA range or an SA singleton.
 *
 * Denoted (A, T), where the set of taxa T is monophyletic and the last occurrence of A
 * (or the singleton A) is a direct sampled ancestor of the MRCA of T.
 *
 * This combines two MRCA cases: the last occurrence of a sampled ancestor range, and a
 * singleton sampled ancestor.
 *
 * @author Alexandra Gavryushkina
 */
public class SampledAncestorClade implements SRClade {

    private final String ancestorTaxon;       // A - the sampled ancestor (last occurrence or singleton)
    private final Set<String> descendantTaxa; // T - monophyletic group below A
    private int count;
    private double probability;

    public SampledAncestorClade(String ancestorTaxon, Set<String> descendantTaxa) {
        this.ancestorTaxon = ancestorTaxon;
        this.descendantTaxa = new TreeSet<>(descendantTaxa);
        this.count = 0;
        this.probability = 0.0;
    }

    public String getAncestorTaxon() {
        return ancestorTaxon;
    }

    public Set<String> getDescendantTaxa() {
        return descendantTaxa;
    }

    @Override
    public String getType() {
        return "sampled_ancestor";
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
        SampledAncestorClade that = (SampledAncestorClade) o;
        return Objects.equals(ancestorTaxon, that.ancestorTaxon) &&
               Objects.equals(descendantTaxa, that.descendantTaxa);
    }

    @Override
    public int hashCode() {
        return Objects.hash(ancestorTaxon, descendantTaxa);
    }

    @Override
    public String toString() {
        return "(" + ancestorTaxon + ", " + descendantTaxa + ")";
    }
}
