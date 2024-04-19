package sr.util;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.util.Objects;

import javax.swing.BoxLayout;
import javax.swing.GroupLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.WindowConstants;
import javax.swing.border.EtchedBorder;

import beastfx.app.treeannotator.TreeAnnotator;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;

public class SRangesTreeOrientator extends TreeAnnotator {

    private static class SRangesTreeOrientatorOptions extends TreeAnnotator {
        File inFile;
        File outFile = new File("summary.tree");

        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n";

        }
    }

    public SRangesTreeOrientator(SRangesTreeOrientatorOptions options) throws IOException {


        // Display options:
        System.out.println(options + "\n");

        // Initialise reader

        TreeSet treeSet = new FastTreeSet(options.inFile.toString(), 0);
        PrintStream ps = new PrintStream(options.outFile);
        ps.println("#NEXUS");
        ps.println("Begin trees;");
        int i = 0;
        while (treeSet.hasNext()) {
            Tree tree = treeSet.next();
            orientateNodeChildren(tree, tree.getRoot().getNr());
            ps.println("tree STATE_" + i + " = " + tree.getRoot().toNewick() + ";");
            i += 1;
        }
        ps.println("End;");
        System.out.println("\nDone!");


    }

    public static void orientateNodeChildren(Tree tree, int subtreeRootNr) {
        Node subTreeRoot = tree.getNode(subtreeRootNr);
        if (!subTreeRoot.isLeaf()) {
            if(subTreeRoot.getChild(0).metaDataString == null)
                return;
            if ((!subTreeRoot.isFake() && !subTreeRoot.getChild(0).metaDataString.contains("orientation=ancestor"))
                    || (subTreeRoot.isFake() && subTreeRoot.getChild(1).getHeight() != subTreeRoot.getHeight())) {
                Node left = subTreeRoot.getChild(1);
                Node right = subTreeRoot.getChild(0);


                subTreeRoot.removeAllChildren(false);

//                subTreeRoot.setLeft(left);
//                subTreeRoot.setRight(right);
                subTreeRoot.addChild(left);
                subTreeRoot.addChild(right);
            }
            if (!subTreeRoot.isFake()){
                if(subTreeRoot.getChild(0).metaDataString.contains("orientation=descendant"))
                    subTreeRoot.getChild(0).metaDataString = subTreeRoot.getChild(0).metaDataString.replace("orientation=descendant", "orientation=ancestor");
                if(subTreeRoot.getChild(1).metaDataString.contains("orientation=ancestor"))
                    subTreeRoot.getChild(1).metaDataString = subTreeRoot.getChild(1).metaDataString.replace("orientation=ancestor", "orientation=descendant");
            }
            orientateNodeChildren(tree, subTreeRoot.getChild(0).getNr());
            orientateNodeChildren(tree, subTreeRoot.getChild(1).getNr());
        }
    }
//    private void orientateNodeChildren(Tree tree, int nodeNr) {
//        Node subTreeRoot = tree.getNode(nodeNr);
//        if (!subTreeRoot.isLeaf()) {
//            if (!Tools.equalHeightWithPrecision(subTreeRoot.getRight(), subTreeRoot)
//                    && !Tools.equalHeightWithPrecision(subTreeRoot.getLeft(), subTreeRoot)) {
//                Node left = Randomizer.nextBoolean() ? subTreeRoot.getRight() : subTreeRoot.getLeft();
//                Node right = left.getNr() == subTreeRoot.getLeft().getNr() ? subTreeRoot.getRight()
//                        : subTreeRoot.getLeft();
//
//                subTreeRoot.removeAllChildren(false);
//
//                subTreeRoot.addChild(left);
//                subTreeRoot.addChild(right);
//
//                left.metaDataString = "orientation=ancestor";
//                right.metaDataString = "orientation=descendant";
//            } else if (Tools.equalHeightWithPrecision(subTreeRoot.getLeft(), subTreeRoot)) {
//                Node left = subTreeRoot.getRight();
//                Node right = subTreeRoot.getLeft();
//
//                subTreeRoot.removeAllChildren(false);
//
//                subTreeRoot.addChild(left);
//                subTreeRoot.addChild(right);
//
//                left.metaDataString = subTreeRoot.metaDataString;
//                right.metaDataString = subTreeRoot.metaDataString;
//            } else if (Tools.equalHeightWithPrecision(subTreeRoot.getRight(), subTreeRoot)) {
//                subTreeRoot.getLeft().metaDataString = subTreeRoot.metaDataString;
//                subTreeRoot.getRight().metaDataString = subTreeRoot.metaDataString;
//            }
//
//            orientateNodeChildren(tree, subTreeRoot.getLeft().getNr());
//            orientateNodeChildren(tree, subTreeRoot.getRight().getNr());
//        }



    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(SRangesTreeOrientatorOptions options) {

        boolean[] canceled = { false };

        JDialog dialog = new JDialog((JDialog) null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("Isolation with Migration Annotator");

        JLabel logFileLabel = new JLabel("Isolation with migration species tree file:");
        JLabel outFileLabel = new JLabel("Output file:");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");

        Container cp = dialog.getContentPane();
        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
        cp.setLayout(boxLayout);

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);
        layout.setAutoCreateGaps(true);
        layout.setAutoCreateContainerGaps(true);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(outFileLabel))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFileButton)
                        .addComponent(outFileButton)));

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(inFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(inFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton)));

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            dialog.setVisible(false);
        });
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

        JFileChooser inFileChooser = new JFileChooser();
        inFileButton.addActionListener(e -> {
            inFileChooser.setDialogTitle("Select Transmission Trees to orientate");
            if (options.inFile == null)
                inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = inFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.inFile = inFileChooser.getSelectedFile();
                inFilename.setText(inFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
            outFileChooser.setDialogTitle("Select output file name.");
            outFileChooser.setCurrentDirectory(Objects.requireNonNullElseGet(options.inFile, () -> new File(System.getProperty("user.dir"))));

            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
            }
        });

        cp.add(buttonPanel);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        return !canceled[0];
    }

    /**
     * Prepare JFrame to which tntOrientator output streams will be directed.
     */
    private static void setupGUIOutput() {

        JFrame frame = new JFrame();
        frame.setTitle("Transmission Tree Orientator");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JTextArea textArea = new JTextArea(25, 80);
        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
        textArea.setEditable(false);
        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> System.exit(0));
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

        // Redirect streams to output window:
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) {
                SwingUtilities.invokeLater(() -> {
                    if ((char) b == '\r') {
                        int from = textArea.getText().lastIndexOf("\n") + 1;
                        int to = textArea.getText().length();
                        textArea.replaceRange(null, from, to);
                    } else
                        textArea.append(String.valueOf((char) b));
                });
            }
        };

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));

        frame.pack();
        frame.setVisible(true);
    }

    public static String helpMessage = "Transmission Tree Orientator - orientates transmission from TNT package tree based on the metadata.\n"
            + "\n"
            + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
            + "\n"
            + "Option                   Description\n"
            + "--------------------------------------------------------------\n"
            + "-help                    Display usage info.\n"
            + "\n"
            + "If no output file is specified, output is written to a file\n"
            + "named 'summary.tree'.";

    /**
     * Print usage info and exit.
     */
    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    /**
     * Display error, print usage and exit with error.
     */
    public static void printUsageAndError(String errMsg) {
        System.err.println(errMsg);
        System.err.println(helpMessage);
        System.exit(1);
    }

    /**
     * Retrieve ACGAnnotator options from command line.
     *
     * @param args    command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, SRangesTreeOrientatorOptions options) {
        int i = 0;
        while (args[i].startsWith("-")) {
            switch (args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;
                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }

        if (i >= args.length)
            printUsageAndError("No input file specified.");
        else
            options.inFile = new File(args[i]);

        if (i + 1 < args.length)
            options.outFile = new File(args[i + 1]);
    }

    /**
     * Main method for TnTOrientator. Sets up GUI if needed then uses the
     * TnTOrientator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
        SRangesTreeOrientatorOptions options = new SRangesTreeOrientatorOptions();

        if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException
                     | IllegalAccessException e) {
                Log.warning.println("Error setting cross-platform look and feel.");
            }

            try {
                SwingUtilities.invokeAndWait(() -> {
                    if (!getOptionsGUI(options))
                        System.exit(0);

                    setupGUIOutput();
                });
            } catch (InterruptedException | InvocationTargetException e) {
                e.printStackTrace();
            }

        } else {
            getCLIOptions(args, options);
        }

        // Run ACGAnnotator
        try {
            new SRangesTreeOrientator(options);
        } catch (Exception e) {
            if (args.length == 0) {
                JOptionPane.showMessageDialog(null, e.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                System.err.println("Error: " + e.getMessage());
                e.printStackTrace();
                System.err.println();
                System.err.println(helpMessage);
            }

            System.exit(1);
        }
    }

}

