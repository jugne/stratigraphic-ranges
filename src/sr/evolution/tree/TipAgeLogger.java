package sr.evolution.tree;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Ugne Stolz
 */
public class TipAgeLogger extends CalculationNode implements Loggable {
    public Input<SRTree> treeInput = new Input<>("tree",
            "sRange tree for range age logging.",
            Input.Validate.REQUIRED);

    public Input<Tree> simpleTreeInput = new Input<>("simpleTree",
            "tree for range age logging.",
            Input.Validate.REQUIRED);

//    public Input<Boolean> relogInput = new Input<>("relog",
//            "If true, this logger is run after the analysis completes. " +
//                    "Default false.",
//            false);



//    public Input<Boolean> onlyFirstInput = new Input<>("onlyFirst",
//            "If true, only the first descendant is logged " ,
//            Boolean.FALSE);
    HashMap<String, Integer> rangeIdMap = new HashMap<>();
    List<String> keys = new ArrayList<>();
    int nRanges = 0;
    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void init(PrintStream out) {
        final SRTree tree = treeInput.get();
        if (tree != null) {
            printToLog(tree, out);
        } else {
            Tree simpleTree = simpleTreeInput.get();
            if (simpleTree != null) {
                printToLog(simpleTree, out);
            }
        }
    }

    private void printToLog(Tree t, PrintStream out){
        for (Node n : t.getExternalNodes()){
            out.print(n.getID() + "\t");
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        final SRTree tree = treeInput.get();
        if (tree != null) {
            tree.orientateTree();
            for (Node n : tree.getExternalNodes()){
                out.print(n.getHeight() + "\t");
            }
        } else if (simpleTreeInput.get() != null) {
            Tree simpleTree = simpleTreeInput.get();
            if (simpleTree != null) {
                for (Node n : tree.getExternalNodes()){
                    out.print(n.getHeight() + "\t");
                }
            }
        }

    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }
}

