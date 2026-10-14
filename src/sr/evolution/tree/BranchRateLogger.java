package sr.evolution.tree;


import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.branchratemodel.BranchRateModel;
import sr.evolution.sranges.StratigraphicRange;

import java.io.PrintStream;
import java.util.Locale;
import java.util.Objects;

import static sr.util.Tools.removeLastSubstring;

@Description("Per-branch logger: length, rate, and range membership (childLabel; 1/0 for isRange).")
public class BranchRateLogger extends BEASTObject implements Loggable {

    // ---- Inputs ----
    public final Input<SRTree> treeInput = new Input<>(
            "tree", "Tree to log per-branch statistics for.", Validate.REQUIRED);

    public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>(
            "branchrateModel",
            "Clock / branch-rate model to retrieve per-branch rates (optional).",
            Validate.OPTIONAL);


    // ---- State ----
    private SRTree tree;
    private BranchRateModel.Base clock;

    // Fixed metadata key per request
    private static final String RANGE_KEY = "range";

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        clock = branchRateModelInput.get();
    }

    @Override
    public void init(PrintStream out) {
        out.println("childLabel\tlength\trate\tisRange\trangeName");
    }

    @Override
    public void log(long sample, PrintStream out) {
        final int nodeCount = tree.getNodeCount();

            for (int i = 0; i < nodeCount; i++) {
                Node child = tree.getNode(i);
                Node parent = child.getParent();

                if (parent == null || child.isFake()) continue; // root has no branch

                double rate = (clock == null) ? Double.NaN : clock.getRateForBranch(child);
                if (parent.isFake()) {
                    rate = (clock == null) ? Double.NaN : clock.getRateForBranch(parent);
                    parent = parent.getParent(); // skip fake parent
                    if (parent == null) {
                        continue; // no parent, skip this child
                    }
                }

                double length = Math.max(0.0, parent.getHeight() - child.getHeight());


                String childLabel = child.isLeaf()
                        ? child.getID()
                        : Integer.toString(child.getNr());

                // the branch above the first occurrence is not part of the range
                String rangeId = null; // no range for this node
                StratigraphicRange range = tree.getRangeOfNode(child);
                if (range != null && !Objects.equals(childLabel, range.getFirstOccurrenceID())) {
                    rangeId = removeLastSubstring("_", range.getLastOccurrenceID());
                }

                boolean isRange = rangeId != null && !rangeId.isEmpty();

                out.printf(Locale.ROOT,
                        "%s\t%.12g\t%12g\t%d\t%s%n",
                        childLabel,
                        length,
                        rate,
                        isRange ? 1 : 0,
                        nz(rangeId));
            }
    }

    @Override
    public void close(PrintStream out) {
        // no-op
    }

    // ---- helpers ----
    private static String nz(String s) {
        return (s == null) ? "" : s;
    }
}
