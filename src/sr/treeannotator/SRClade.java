package sr.treeannotator;

/**
 * Common interface for the three types of clades defined for stratigraphic-range (SR) trees.
 *
 * The clade definition is based on the MRCA of a set of taxa computed in the SA tree obtained
 * from the SR tree by omitting all nodes that correspond to the first occurrences of ranges
 * (singletons are kept). Depending on the type of node that is the MRCA, three clade types arise:
 *
 * <ol>
 *     <li>Bifurcation event not within a range: (T1, T2)</li>
 *     <li>Bifurcation event within a range: (A, T1, T2)</li>
 *     <li>SA range (last occurrence) or SA singleton: (A, T)</li>
 * </ol>
 *
 * @author Alexandra Gavryushkina
 */
public interface SRClade {

    /** A short label identifying the clade type, used for annotation. */
    String getType();

    int getCount();

    void setCount(int count);

    void incrementCount();

    double getProbability();

    void setProbability(double probability);
}
