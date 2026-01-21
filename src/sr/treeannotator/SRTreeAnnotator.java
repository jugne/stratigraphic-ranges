package sr.treeannotator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.Runnable;
import beast.base.parser.NexusParser;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import sr.evolution.tree.SRTree;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Tree annotator for stratigraphic-range (SR) trees using relationship-based
 * credibility instead of traditional clade-based credibility.
 *
 * Implements the approach from section 1.1.3 where:
 * - Relationships (ancestry and orientation) are collected instead of clades
 * - MCC tree maximizes the product of relationship probabilities
 *
 * @author Ugne Stolz
 */
@Description("Summarizes posterior tree samples for SR trees using relationship-based " +
        "credibility instead of traditional clade-based methods. Preserves orientation information.")
public class SRTreeAnnotator extends Runnable {

    final public Input<TreeFile> treesInput = new Input<>("trees",
            "Input tree file (Nexus or Newick format)",
            new TreeFile("[[none]]"));

    final public Input<OutFile> outputInput = new Input<>("out",
            "Output file for the MCC tree (Nexus format)",
            new OutFile("[[none]]"));

    final public Input<Integer> burnInPercentageInput = new Input<>("burnin",
            "Percentage of trees to discard as burn-in",
            10);

    final public Input<Boolean> useSumInput = new Input<>("sumProbabilities",
            "Use sum of probabilities instead of log product for scoring",
            false);

    final public Input<OutFile> summaryInput = new Input<>("summary",
            "Optional output file for relationship summary",
            new OutFile("[[none]]"));

    final public Input<Boolean> detailedInput = new Input<>("detailed",
            "Include detailed relationship annotations (ancestral/descendant taxa)",
            false);

    private String inputFileName;
    private String outputFileName;
    private String summaryFileName;
    private int burninPercentage;
    private boolean useSumCredibility;
    private boolean annotateRelationshipDetails;

    @Override
    public void initAndValidate() {
        // Validation will happen in run()
    }

    @Override
    public void run() throws Exception {
        // Get values from inputs
        inputFileName = treesInput.get().getPath();
        outputFileName = outputInput.get().getPath();
        burninPercentage = burnInPercentageInput.get();
        useSumCredibility = useSumInput.get();
        annotateRelationshipDetails = detailedInput.get();

        if (summaryInput.get() != null && !summaryInput.get().getName().equals("[[none]]")) {
            summaryFileName = summaryInput.get().getPath();
        } else {
            summaryFileName = null;
        }

        // Validate inputs
        if (inputFileName == null || inputFileName.equals("[[none]]")) {
            throw new IllegalArgumentException("Input tree file must be specified");
        }
        if (outputFileName == null || outputFileName.equals("[[none]]")) {
            throw new IllegalArgumentException("Output file must be specified");
        }

        annotate();
    }

    /**
     * Main annotation process.
     */
    public void annotate() throws Exception {
        Log.info("SR Tree Annotator");
        Log.info("Using relationship-based credibility (section 1.1.3)");
        Log.info("");

        // Read trees from file
        List<SRTree> trees = readTrees();

        if (trees.isEmpty()) {
            throw new IllegalArgumentException("No trees found in input file");
        }

        Log.info("Total trees in file: " + trees.size());

        // Apply burnin
        int burninCount = (burninPercentage * trees.size()) / 100;
        List<SRTree> analyzedTrees = trees.subList(burninCount, trees.size());
        int totalTreesUsed = analyzedTrees.size();

        Log.info("Burnin: " + burninPercentage + "% (" + burninCount + " trees)");
        Log.info("Trees to analyze: " + totalTreesUsed);
        Log.info("");

        // Phase 1: Collect relationships and compute probabilities
        Log.info("Step 1: Collecting relationships and attributes from trees...");
        RelationshipSystem relationshipSystem = new RelationshipSystem();

        int counter = 0;
        for (SRTree tree : analyzedTrees) {
            // Collect relationships WITH heights for annotation
            relationshipSystem.add(tree, true);
            counter++;
            if (counter % 100 == 0) {
                Log.info.print(".");
                if (counter % 1000 == 0) {
                    Log.info.print(" " + counter);
                }
            }
        }
        Log.info("");
        Log.info("Collected relationships from " + totalTreesUsed + " trees.");

        // Calculate posterior probabilities
        relationshipSystem.calculatePosteriorProbabilities(totalTreesUsed);

        // Write relationship summary to file if specified
        if (summaryFileName != null) {
            writeSummary(relationshipSystem, totalTreesUsed);
            Log.info("Relationship summary written to: " + summaryFileName);
        }

        // Step 2: Find MCC tree
        Log.info("Step 2: Finding maximum credibility tree...");
        Log.info("0              25             50             75            100");
        Log.info("|--------------|--------------|--------------|--------------|");

        SRTree bestTree = null;
        double bestScore = Double.NEGATIVE_INFINITY;

        counter = 0;
        int reported = 0;

        for (SRTree tree : analyzedTrees) {
            double score;
            if (useSumCredibility) {
                score = relationshipSystem.getSumRelationshipCredibility(tree);
            } else {
                score = relationshipSystem.getLogRelationshipCredibility(tree);
            }

            if (score > bestScore) {
                bestTree = tree;
                bestScore = score;
            }

            // Progress bar
            while (reported < 61 && 1000.0 * reported < 61000.0 * (counter + 1) / totalTreesUsed) {
                Log.info.print("*");
                reported++;
            }
            counter++;
        }

        Log.info("");
        Log.info("");

        if (useSumCredibility) {
            Log.info("Highest Sum Relationship Credibility: " + bestScore);
        } else {
            Log.info("Highest Log Relationship Credibility: " + bestScore);
        }

        // Step 3: Annotate MCC tree
        if (bestTree != null) {
            Log.info("\nStep 3: Annotating MCC tree with relationship probabilities and statistics...");
            relationshipSystem.annotateMCCTree(bestTree, annotateRelationshipDetails);

            writeTree(bestTree);
            Log.info("MCC tree written to: " + outputFileName);
        } else {
            Log.err("ERROR: No best tree found");
        }

        Log.info("\nDone!");
    }

    /**
     * Reads SR trees from input file.
     */
    private List<SRTree> readTrees() throws Exception {
        Log.info("Reading trees from: " + inputFileName);

        List<SRTree> trees = new ArrayList<>();

        // Try Nexus format first
        try {
            NexusParser parser = new NexusParser();
            parser.parseFile(new File(inputFileName));

            for (Tree tree : parser.trees) {
                if (tree instanceof SRTree) {
                    trees.add((SRTree) tree);
                } else {
                    // Try to convert to SRTree
                    SRTree srTree = new SRTree();
                    srTree.assignFrom(tree);
                    srTree.orientateTree();
                    trees.add(srTree);
                }
            }
        } catch (Exception e) {
            // If Nexus parsing fails, try reading as Newick
            Log.info("Nexus parsing failed, trying Newick format...");
            trees = readNewickTrees();
        }

        return trees;
    }

    /**
     * Reads trees from Newick format file using TreeParser.
     */
    private List<SRTree> readNewickTrees() throws IOException {
        List<SRTree> trees = new ArrayList<>();
        BufferedReader reader = new BufferedReader(new FileReader(inputFileName));

        String line;
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }

            try {
                TreeParser parser = new TreeParser(line);
                SRTree srTree = new SRTree();
                srTree.assignFrom(parser);
                srTree.orientateTree();
                trees.add(srTree);
            } catch (Exception e) {
                Log.warning("Warning: Failed to parse tree: " + e.getMessage());
            }
        }

        reader.close();
        return trees;
    }

    /**
     * Writes the MCC tree to output file in proper Nexus format with taxa labels.
     */
    private void writeTree(SRTree tree) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(outputFileName));

        // Collect all leaf nodes (taxa)
        List<String> taxaNames = new ArrayList<>();
        collectTaxaNames(tree.getRoot(), taxaNames);

        // Write Nexus header
        writer.println("#NEXUS");
        writer.println();

        // Write Taxa block
        writer.println("Begin taxa;");
        writer.println("\tDimensions ntax=" + taxaNames.size() + ";");
        writer.println("\t\tTaxlabels");
        for (String taxon : taxaNames) {
            writer.println("\t\t\t" + taxon);
        }
        writer.println("\t\t\t;");
        writer.println("End;");
        writer.println();

        // Write Trees block with Translate
        writer.println("Begin trees;");
        writer.println("\tTranslate");
        for (int i = 0; i < taxaNames.size(); i++) {
            String comma = (i < taxaNames.size() - 1) ? "," : "";
            writer.println("\t\t" + (i + 1) + " " + taxaNames.get(i) + comma);
        }
        writer.println(";");

        // Write the tree using toNewick() which includes metadata
        writer.println("tree MCC_tree = " + tree.getRoot().toNewick() + ";");
        writer.println("End;");

        writer.close();
    }

    /**
     * Recursively collects all taxa names from the tree.
     */
    private void collectTaxaNames(Node node, List<String> taxaNames) {
        if (node.isLeaf()) {
            taxaNames.add(node.getID());
        } else {
            if (node.getLeft() != null) {
                collectTaxaNames(node.getLeft(), taxaNames);
            }
            if (node.getRight() != null) {
                collectTaxaNames(node.getRight(), taxaNames);
            }
        }
    }

    /**
     * Writes the relationship summary to a file.
     */
    private void writeSummary(RelationshipSystem relationshipSystem, int totalTreesUsed) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(summaryFileName));
        writer.println("SR Tree Annotator - Relationship Summary");
        writer.println("=========================================");
        writer.println("Total trees analyzed: " + totalTreesUsed);
        writer.println();
        writer.println(relationshipSystem.getSummary());
        writer.close();
    }

    /**
     * Command-line interface via BEAST Application framework.
     */
    public static void main(String[] args) throws Exception {
        new Application(new SRTreeAnnotator(), "SR Tree Annotator", args);
    }
}
