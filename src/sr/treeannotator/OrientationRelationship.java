package sr.treeannotator;

import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

/**
 * Represents an orientation relationship (T1, T2) where:
 * - T1 are taxa descending from the ancestral (left) lineage
 * - T2 are taxa descending from the descendant (right) lineage
 *
 * This can occur in two cases:
 * - Case 2a: Standard bifurcation where T = T1 ∪ T2 is monophyletic and the MRCA is a bifurcation event
 * - Case 2b: Bifurcation within a range where there exists taxon A such that {A} ∪ T is monophyletic
 *
 * This relationship captures the orientation of speciation events.
 * Note: Order matters - T1 is ancestral, T2 is descendant.
 *
 * @author Ugne Stolz
 */
public class OrientationRelationship {

    private final Set<String> ancestralTaxa;  // T1 - taxa on ancestral lineage
    private final Set<String> descendantTaxa; // T2 - taxa on descendant lineage
    private int count;
    private double probability;

    /**
     * Creates a new orientation relationship.
     *
     * @param ancestralTaxa Set of taxa descending from ancestral (left) lineage (T1)
     * @param descendantTaxa Set of taxa descending from descendant (right) lineage (T2)
     */
    public OrientationRelationship(Set<String> ancestralTaxa, Set<String> descendantTaxa) {
        // Use TreeSet for consistent ordering within each set
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

    public int getCount() {
        return count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public void incrementCount() {
        this.count++;
    }

    public double getProbability() {
        return probability;
    }

    public void setProbability(double probability) {
        this.probability = probability;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        OrientationRelationship that = (OrientationRelationship) o;

        // Order matters: T1 is ancestral, T2 is descendant
        return Objects.equals(ancestralTaxa, that.ancestralTaxa) &&
               Objects.equals(descendantTaxa, that.descendantTaxa);
    }

    @Override
    public int hashCode() {
        // Order matters in the hash
        return Objects.hash(ancestralTaxa, descendantTaxa);
    }

    @Override
    public String toString() {
        return "(" + ancestralTaxa + " => " + descendantTaxa + ")";
    }
}
