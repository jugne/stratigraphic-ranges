package sr.evolution.tree;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import sr.evolution.sranges.StratigraphicRange;
import sr.util.Tools;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import static sr.util.Tools.getSeparatingLengthAndNodeCount;
import static sr.util.Tools.removeLastSubstring;

/**
 * @author Ugne Stolz
 */
public class TipAgeLogger extends CalculationNode implements Loggable, Function {
    public Input<SRTree> treeInput = new Input<>("tree",
            "sRange tree for range age logging.",
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
        for (Node n : tree.getExternalNodes()){
            out.print(n.getID() + "\t");
//            keys.add(n.getID());
        }
//        for (String key :keys) {
//            out.print(key + "\t");
//        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        final SRTree tree = treeInput.get();
        tree.orientateTree();
        for (Node n : tree.getExternalNodes()){
            out.print(n.getHeight() + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        return treeInput.get().getDirectAncestorNodeCount();
    }

    @Override
    public double getArrayValue(int iDim) {
        return treeInput.get().getDirectAncestorNodeCount();
    }

}

