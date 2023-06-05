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

import static sr.util.Tools.removeLastSubstring;

/**
 * @author Alexandra Gavryushkina
 */
public class SpeciationLogger extends CalculationNode implements Loggable, Function {
    public Input<SRTree> treeInput = new Input<>("tree",
            "tree to report SA count for.",
            Input.Validate.REQUIRED);
    public Input<String> sepStringInput = new Input<>("sep",
            "separator string for ranges",
            "_");
    public Input<String> directStringInput = new Input<>("dir",
            "string to indicate speciation direction. " +
                    "It will separate the donor range and recipient species.",
            ">");

    public Input<Boolean> onlyFirstInput = new Input<>("onlyFirst",
            "If true, only the first descendant is logged " ,
            Boolean.TRUE);
    HashMap<String, Integer> speciationsCounter = new HashMap<>();
    List<String> keys = new ArrayList<>();
    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void init(PrintStream out) {
        final SRTree tree = treeInput.get();
        ArrayList<StratigraphicRange> ranges = tree.getSRanges();
        for (StratigraphicRange range : ranges) {
            for(Node leaf : tree.getExternalNodes()){
                String key = removeLastSubstring(sepStringInput.get(), range.getFirstOccurrenceID()) +
                        directStringInput.get() + removeLastSubstring(sepStringInput.get(), leaf.getID());
                if (speciationsCounter.get(key)==null){
                    speciationsCounter.put(key, 0);
                    keys.add(key);
                    out.print(key + "\t");
                }
//                    System.out.println("key = " + key);
//                speciationsCounter.putIfAbsent(key, 0);
//                keys.add(key);
//                out.print(key + "\t");
            }
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        final SRTree tree = treeInput.get();
        for (StratigraphicRange range : tree.getSRanges()) {
            for (Integer i : range.getInternalNodeNrs(tree)) {
                Node right = tree.getNode(i);
                if (!tree.getNode(i).isLeaf())
                    right = tree.getNode(i).getRight();
                if (onlyFirstInput.get() && !tree.getNode(i).isLeaf()) {
                    String key = removeLastSubstring(sepStringInput.get(), range.getFirstOccurrenceID()) +
                            directStringInput.get()  + removeLastSubstring(sepStringInput.get(), Tools.getFirstLeafID(right));
                    if (speciationsCounter.get(key)==null)
                        System.out.println("key = " + key);
                    speciationsCounter.put(key, 1);
                } else {
                    for (Node l : right.getAllLeafNodes()){
                        String key = removeLastSubstring(sepStringInput.get(),range.getFirstOccurrenceID()) +
                                directStringInput.get()  + removeLastSubstring(sepStringInput.get(),l.getID());
                        if (speciationsCounter.get(key)==null)
                            System.out.println("key = " + key);
                        speciationsCounter.put(key, 1);
                    }
                }
            }
        }
        for (String s : keys) {
            out.print(speciationsCounter.get(s) + "\t");
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

