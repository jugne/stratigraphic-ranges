package sr.util;

import beast.base.core.*;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.Parameter;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

@Description("Based on the SpeciesTreeLogger class, but without node sorting")
public class TreeWithMetadataLogger extends BEASTObject implements Loggable {
	final public Input<SRTree> srTeeInput = new Input<>("tree",
			"The range tree to be logged.", Input.Validate.REQUIRED);
    final public Input<BranchRateModel> clockModelInput = new Input<>("branchratemodel", "rate to be logged with branches of the tree");
    final public Input<Boolean> substitutionsInput = new Input<>("substitutions", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    final public Input<Integer> decimalPlacesInput = new Input<>("dp", "the number of decimal places to use writing branch lengths and rates, use -1 for full precision (default = full precision)", -1);
    final public Input<List<Function>> parameterInput = new Input<>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());

	final public Input<Boolean> logOrientationInput = new Input<>("logOrientation",
			"report if node is donor or recipient", true);

    boolean someMetaDataNeedsLogging;
    boolean substitutions = false;

    private DecimalFormat df;

    @Override
    public void initAndValidate() {
		if (parameterInput.get().size() == 0 && clockModelInput.get() == null && !logOrientationInput.get()) {
            someMetaDataNeedsLogging = false;
            return;
        }
		someMetaDataNeedsLogging = true;

        // without substitution model, reporting substitutions == reporting branch lengths 
        if (clockModelInput.get() != null) {
            substitutions = substitutionsInput.get();
        }

        int dp = decimalPlacesInput.get();

        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }
    }

    @Override
    public void init(PrintStream out) {
		SRTree srTree = srTeeInput.get();
        srTree.init(out);
    }

    @Override
    public void log(long nSample, PrintStream out) {
        // make sure we get the current version of the inputs
        SRTree srTree = (SRTree) srTeeInput.get().getCurrent();

        srTree.addOrientationMetadata();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
            if (metadata.get(i) instanceof StateNode) {
                metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
            }
        }
        BranchRateModel branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
//        tree.getRoot().sort();
		out.print(toNewick((SRNode) srTree.getRoot(), metadata, branchRateModel));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
    }

    /**
     * Appends a double to the given StringBuffer, formatting it using
     * the private DecimalFormat instance, if the input 'dp' has been
     * given a non-negative integer, otherwise just uses default
     * formatting.
     * @param buf
     * @param d
     */
    private void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
    }

	String toNewick(SRNode node, List<Function> metadataList, BranchRateModel branchRateModel) {
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick((SRNode) node.getLeft(), metadataList, branchRateModel));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick((SRNode) node.getRight(), metadataList, branchRateModel));
            }
            buf.append(")");
        } else {
			buf.append(node.getNr() + 1);
        }
        if (someMetaDataNeedsLogging) {
			if (node.getID() == null) {
				buf.append(node.getNr() + 1);
			}
            buf.append("[&");
            if (metadataList.size() > 0) {
                for (Function metadata : metadataList) {
                    buf.append(((BEASTObject)metadata).getID());
                    buf.append('=');
                    if (metadata instanceof Parameter<?>) {
                        Parameter<?> p = (Parameter<?>) metadata;
                        int dim = p.getMinorDimension1();
                        if (dim > 1) {
                            buf.append('{');
                            for (int i = 0; i < dim; i++) {
                                buf.append(p.getMatrixValue(node.getNr(), i));
                                if (i < dim - 1) {
                                    buf.append(',');
                                }
                            }
                            buf.append('}');
                        } else {
                            buf.append(metadata.getArrayValue(node.getNr()));
                        }
                    } else {
                        buf.append(metadata.getArrayValue(node.getNr()));
                    }
                    if (metadataList.indexOf(metadata) < metadataList.size() - 1) {
                        buf.append(",");
                    }
                }
                if (branchRateModel != null) {
                    buf.append(",");
                }
            }
            if (branchRateModel != null) {
                buf.append("rate=");
                appendDouble(buf, branchRateModel.getRateForBranch(node));
				if (logOrientationInput.get()) {
                    buf.append(",");
                }
            }

			if (logOrientationInput.get()) {
				buf.append(node.metaDataString);
            }

            buf.append(']');
        }

        buf.append(":");

        double nodeLength = node.getLength();

        if (substitutions) {
            appendDouble(buf, nodeLength * branchRateModel.getRateForBranch(node));
        } else {
            appendDouble(buf, nodeLength);
        }
        return buf.toString();
    }


    @Override
    public void close(PrintStream out) {
		SRTree tree = srTeeInput.get();
        tree.close(out);
    }

}

