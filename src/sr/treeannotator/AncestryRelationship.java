package sr.treeannotator;

import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

/**
 * Represents an ancestry relationship (A, T) where:
 * - T is a monophyletic set of taxa
 * - The first occurrence of taxon A (or only occurrence if singleton) is a direct ancestor of the MRCA of T
 *
 * This relationship captures "A is ancestral to group T" and covers cases where a species range
 * or singleton is a sampled ancestor.
 *
 * @author Ugne Stolz
 */
public class AncestryRelationship {

    private final String ancestorTaxon;
    private final Set<String> descendantTaxa;
    private int count;
    private double probability;

    /**
     * Creates a new ancestry relationship.
     *
     * @param ancestorTaxon The ancestor taxon name (A)
     * @param descendantTaxa The set of descendant taxa (T)
     */
    public AncestryRelationship(String ancestorTaxon, Set<String> descendantTaxa) {
        this.ancestorTaxon = ancestorTaxon;
        // Use TreeSet for consistent ordering and equality checking
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

        AncestryRelationship that = (AncestryRelationship) o;

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