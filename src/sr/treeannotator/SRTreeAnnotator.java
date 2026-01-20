package sr.treeannotator;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.parser.NexusParser;
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
public class SRTreeAnnotator {

    private int burninPercentage = 10;
    private String inputFileName;
    private String outputFileName;
    private String summaryFileName = null;  // Optional summary output file
    private PrintStream progressStream = System.out;
    private boolean useSumCredibility = false;
    private boolean annotateRelationshipDetails = false;  // Whether to add ancestral/descendant taxa annotations

    public SRTreeAnnotator(String inputFileName, String outputFileName) {
        this.inputFileName = inputFileName;
        this.outputFileName = outputFileName;
    }

    public void setBurninPercentage(int burninPercentage) {
        this.burninPercentage = burninPercentage;
    }

    public void setUseSumCredibility(boolean useSumCredibility) {
        this.useSumCredibility = useSumCredibility;
    }

    public void setSummaryFileName(String summaryFileName) {
        this.summaryFileName = summaryFileName;
    }

    public void setAnnotateRelationshipDetails(boolean annotateRelationshipDetails) {
        this.annotateRelationshipDetails = annotateRelationshipDetails;
    }

    /**
     * Main annotation process.
     */
    public void annotate() throws Exception {
        progressStream.println("SR Tree Annotator");
        progressStream.println("Using relationship-based credibility (section 1.1.3)");
        progressStream.println();

        // Read trees from file
        List<SRTree> trees = readTrees();

        if (trees.isEmpty()) {
            throw new IllegalArgumentException("No trees found in input file");
        }

        progressStream.println("Total trees in file: " + trees.size());

        // Apply burnin
        int burninCount = (burninPercentage * trees.size()) / 100;
        List<SRTree> analyzedTrees = trees.subList(burninCount, trees.size());
        int totalTreesUsed = analyzedTrees.size();

        progressStream.println("Burnin: " + burninPercentage + "% (" + burninCount + " trees)");
        progressStream.println("Trees to analyze: " + totalTreesUsed);
        progressStream.println();

        // Phase 1: Collect relationships and compute probabilities
        progressStream.println("Step 1: Collecting relationships and attributes from trees...");
        RelationshipSystem relationshipSystem = new RelationshipSystem();

        int counter = 0;
        for (SRTree tree : analyzedTrees) {
            // Collect relationships WITH heights for annotation
            relationshipSystem.add(tree, true);
            counter++;
            if (counter % 100 == 0) {
                progressStream.print(".");
                if (counter % 1000 == 0) {
                    progressStream.print(" " + counter);
                }
                progressStream.flush();
            }
        }
        progressStream.println();
        progressStream.println("Collected relationships from " + totalTreesUsed + " trees.");

        // Calculate posterior probabilities
        relationshipSystem.calculatePosteriorProbabilities(totalTreesUsed);

        // Write relationship summary to file if specified
        if (summaryFileName != null) {
            writeSummary(relationshipSystem, totalTreesUsed);
            progressStream.println("Relationship summary written to: " + summaryFileName);
        }

        // Step 2: Find MCC tree
        progressStream.println("Step 2: Finding maximum credibility tree...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

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
                progressStream.print("*");
                reported++;
                progressStream.flush();
            }
            counter++;
        }

        progressStream.println();
        progressStream.println();

        if (useSumCredibility) {
            progressStream.println("Highest Sum Relationship Credibility: " + bestScore);
        } else {
            progressStream.println("Highest Log Relationship Credibility: " + bestScore);
        }

        // Step 3: Annotate MCC tree
        if (bestTree != null) {
            progressStream.println("\nStep 3: Annotating MCC tree with relationship probabilities and statistics...");
            relationshipSystem.annotateMCCTree(bestTree, annotateRelationshipDetails);

            writeTree(bestTree);
            progressStream.println("MCC tree written to: " + outputFileName);
        } else {
            progressStream.println("ERROR: No best tree found");
        }

        progressStream.println("\nDone!");
    }

    /**
     * Reads SR trees from input file.
     */
    private List<SRTree> readTrees() throws Exception {
        progressStream.println("Reading trees from: " + inputFileName);

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
            progressStream.println("Nexus parsing failed, trying Newick format...");
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
                progressStream.println("Warning: Failed to parse tree: " + e.getMessage());
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
     * Command-line interface.
     */
    public static void main(String[] args) {
        try {
            if (args.length < 2) {
                printUsage();
                return;
            }

            String inputFile = args[0];
            String outputFile = args[1];
            int burnin = 10;
            boolean useSum = false;
            String summaryFile = null;
            boolean detailed = false;

            // Parse optional arguments
            for (int i = 2; i < args.length; i++) {
                if (args[i].equals("-burnin") && i + 1 < args.length) {
                    burnin = Integer.parseInt(args[i + 1]);
                    i++;
                } else if (args[i].equals("-sum")) {
                    useSum = true;
                } else if (args[i].equals("-summary") && i + 1 < args.length) {
                    summaryFile = args[i + 1];
                    i++;
                } else if (args[i].equals("-detailed")) {
                    detailed = true;
                }
            }

            SRTreeAnnotator annotator = new SRTreeAnnotator(inputFile, outputFile);
            annotator.setBurninPercentage(burnin);
            annotator.setUseSumCredibility(useSum);
            annotator.setAnnotateRelationshipDetails(detailed);
            if (summaryFile != null) {
                annotator.setSummaryFileName(summaryFile);
            }
            annotator.annotate();

        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
            e.printStackTrace();
            printUsage();
        }
    }

    private static void printUsage() {
        System.out.println("SR Tree Annotator - Relationship-based tree summarization for SR trees");
        System.out.println();
        System.out.println("Usage: java sr.treeannotator.SRTreeAnnotator <input_file> <output_file> [options]");
        System.out.println();
        System.out.println("Options:");
        System.out.println("  -burnin <percentage>    Percentage of trees to discard as burnin (default: 10)");
        System.out.println("  -sum                    Use sum of probabilities instead of log (default: log)");
        System.out.println("  -summary <file>         Write relationship summary to specified file");
        System.out.println("  -detailed               Include detailed relationship annotations (ancestral/descendant taxa)");
        System.out.println();
        System.out.println("Example:");
        System.out.println("  java sr.treeannotator.SRTreeAnnotator input.trees output.tree -burnin 20 -summary relationships.txt");
    }
}
