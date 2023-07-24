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
public class SpeciationLogger extends CalculationNode implements Loggable, Function {
    public Input<SRTree> treeInput = new Input<>("tree",
            "sRange tree for speciation logging.",
            Input.Validate.REQUIRED);
    public Input<String> sepStringInput = new Input<>("sep",
            "separator string for ranges",
            "_");
    public Input<String> directStringInput = new Input<>("dir",
            "string to indicate speciation direction. " +
                    "It will separate the donor range and recipient species.",
            ">");

    public Input<Boolean> relogInput = new Input<>("relog",
            "If true, this logger is run after the analysis completes. " +
                    "Default false.",
            false);



//    public Input<Boolean> onlyFirstInput = new Input<>("onlyFirst",
//            "If true, only the first descendant is logged " ,
//            Boolean.FALSE);
    HashMap<String, double[]> speciationsCounter = new HashMap<>();
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
                    double[] vals = new double[2];
                    speciationsCounter.put(key, vals);
                    keys.add(key);
                    out.print(key + "\t");
                }
            }
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        final SRTree tree = treeInput.get();
        if (relogInput.get()){
            tree.orientateTree();
            tree.initSRanges();
        }
        for (StratigraphicRange range : tree.getSRanges()) {
            List<Integer> internalNodeNrs = range.getInternalNodeNrs(tree);
            if (relogInput.get() && !range.isSingleFossilRange()){
                SRNode start = (SRNode) tree.getNode(range.getNodeNrs().get(0));
                SRNode child = (SRNode) start.getParent().getLeft();
                int childNr = child.getNr();
                if (child.isFake())
                    childNr=child.getDirectAncestorChild().getNr();
                while (childNr != internalNodeNrs.get(0)){
                    range.addNodeNrAfter(tree, start.getNr(), child.getNr());
                    start = child;
                    child = (SRNode) start.getLeft();
                    childNr = child.getNr();
                    if (child.isFake())
                        childNr=child.getDirectAncestorChild().getNr();
                }
            }
            if (range.getInternalNodeNrs(tree).size() == 1)
                continue;
            for (Integer i : range.getInternalNodeNrs(tree)) {
                Node right = tree.getNode(i);
                if (!tree.getNode(i).isLeaf())
                    right = tree.getNode(i).getRight();
                if (right.isLeaf() && !right.getID().contains("last")){
                    String key = removeLastSubstring(sepStringInput.get(),range.getFirstOccurrenceID()) +
                            directStringInput.get()  + removeLastSubstring(sepStringInput.get(),right.getID());
                    double[] vals = getSeparatingLengthAndNodeCount(right.getParent(), right);
                    if (right.isDirectAncestor())
                        vals[1] = vals[1] - 1;
                    speciationsCounter.put(key, vals);
                }
                    for (Node l : right.getAllLeafNodes()){
                        if (!l.getID().contains("last")){
                            String key = removeLastSubstring(sepStringInput.get(),range.getFirstOccurrenceID()) +
                                    directStringInput.get()  + removeLastSubstring(sepStringInput.get(),l.getID());
                            double[] vals = getSeparatingLengthAndNodeCount(tree.getNode(i), l);
                            if (l.isDirectAncestor())
                                vals[1] = vals[1] - 1;
                            speciationsCounter.put(key, vals);
                        }
                    }
//                }
            }
        }
        for (String s : keys) {
//            out.print(speciationsCounter.get(s) + "\t");
            out.print(String.join(",", Arrays.stream(speciationsCounter.get(s))
                    .mapToObj(Double::toString)
                    .toArray(String[]::new))+ "\t");
            speciationsCounter.put(s, new double[]{0,0});
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

