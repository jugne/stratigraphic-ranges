package sr.treeannotator;

import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

/**
 * Clade type 1: a bifurcation event that is not within a range.
 *
 * Denoted (T1, T2), where T = T1 ∪ T2 is monophyletic and the MRCA of T is a bifurcation
 * event at which the taxa T1 are descendants of the ancestral (left) lineage and the taxa T2
 * are descendants of the descendant (right) lineage.
 *
 * Orientation matters: T1 is ancestral, T2 is descendant.
 *
 * @author Alexandra Gavryushkina
 */
public class BifurcationClade implements SRClade {

    private final Set<String> ancestralTaxa;  // T1 - ancestral (left) lineage
    private final Set<String> descendantTaxa; // T2 - descendant (right) lineage
    private int count;
    private double probability;

    public BifurcationClade(Set<String> ancestralTaxa, Set<String> descendantTaxa) {
        this.ancestralTaxa = new TreeSet<>(ancestralTaxa);
        this.descendantTaxa = new TreeSet<>(descendantTaxa);
        this.count = 0;
        this.probability = 0.0;
    }

    public Set<String> getAncestralTaxa() {
        return ancestralTaxa;
    }

    public Set<String> getDescendantTaxa() {
        return descendantTaxa;
    }

    @Override
    public String getType() {
        return "bifurcation";
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
        BifurcationClade that = (BifurcationClade) o;
        // Order matters: T1 is ancestral, T2 is descendant
        return Objects.equals(ancestralTaxa, that.ancestralTaxa) &&
               Objects.equals(descendantTaxa, that.descendantTaxa);
    }

    @Override
    public int hashCode() {
        return Objects.hash(ancestralTaxa, descendantTaxa);
    }

    @Override
    public String toString() {
        return "(" + ancestralTaxa + ", " + descendantTaxa + ")";
    }
}
