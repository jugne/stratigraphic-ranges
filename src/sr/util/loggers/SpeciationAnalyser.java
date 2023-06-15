package sr.util.loggers;

import beastfx.app.treeannotator.TreeAnnotator;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import sr.util.Tools;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.util.List;
import java.util.*;

import static sr.util.Tools.removeLastSubstring;

/**
 * Analyses transmission trees from TnT transmission tree log file and
 * outputs their transmission histories as a table. Each row represents a tree in the log file.
 * Each cell is a possible transmission event.
 *
 * If direct transmission tree (full sampling or mapped tree), cell value of:
 * 0 - transmission event not found in the tree;
 * 1 - direct transmission event found in the tree;
 * n - indirect transmission event with n-1 intermediate transmissions from unsampled hosts found in the tree.
 *
 * If indirect transmission tree, cell value of:
 * 0 - transmission event not found in the tree;
 * 1 - direct or indirect (do not know which) transmission event found in the tree;
 *
 * Each host sample has to be named following a scheme: hostName_sampleName
 * hostName can be any number/word, except "0". 0 is reserved for unsampled root in the case of incomplete sampling
 * direct transmission tree.
 * sampleName can be any letter or number (and can also include further uses of character "_").
 *
 *  * @author ugne.stolz@protonmail.com
 *  * @date 04.08.22
 *
 */

public class SpeciationAnalyser extends TreeAnnotator {
    private static class SpeciationAnalyserOptions extends TreeAnnotator {
        File inTransmissionTreeFile;
        File outFile;
        int burnIn;

        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input transmission tree file: " + inTransmissionTreeFile + "\n" +
                    "Output file name: " + outFile + "\n" +
                    "Burn in: " + burnIn + "\n";
        }
    }

    public SpeciationAnalyser(SpeciationAnalyser.SpeciationAnalyserOptions options) throws IOException {
        // Display options:
        System.out.println(options + "\n");

        // Initialise reader
        TreeSet treeSet = new MemoryFriendlyTreeSet(options.inTransmissionTreeFile.toString(), options.burnIn);
        PrintStream ps = new PrintStream(options.outFile);

        String outDir = options.inTransmissionTreeFile.getParent();
        String logFile = options.inTransmissionTreeFile.getName();
        int pos = logFile.lastIndexOf(".");
        String justName = pos >= 0 ? logFile.substring(0, pos) : logFile;
        File infectionTimesOut = new File(outDir, justName+"_speciationTimes.txt");

        PrintStream psInfectionTimes = new PrintStream(infectionTimesOut);

        boolean first = true;
        List<String> hostsList = null;
        List<Integer> skipIDs = new ArrayList<>();
        treeSet.reset();
        while (treeSet.hasNext()) {
            Tree tree = treeSet.next();
            Set<String> hosts = new HashSet<String>();
            if (first) {
                hosts.add("unsampled");
                for (int i=0; i<tree.getLeafNodeCount(); i++) {
                    String hostName = removeLastSubstring("_", tree.getNode(i).getID());
                    hosts.add(hostName);
                }
                hostsList = new ArrayList<>(hosts);
                for (int i = 0; i < hostsList.size(); i++) {
                    if (i == hostsList.size() - 1 && !Objects.equals(hostsList.get(i), "unsampled"))
                        psInfectionTimes.print(hostsList.get(i));
                    else if(!Objects.equals(hostsList.get(i), "unsampled"))
                        psInfectionTimes.print(hostsList.get(i) +"\t");

                    for (int j = 0; j < hostsList.size(); j++) {
                        if (i==j || Objects.equals(hostsList.get(j), "unsampled")){
                            skipIDs.add(i * hostsList.size() + j);
                            continue;
                        }
                        if (i == hostsList.size() - 1 && j == hostsList.size() - 2)
                            ps.print(hostsList.get(i) + "_" + hostsList.get(j));
                        else
                            ps.print(hostsList.get(i) + "_" + hostsList.get(j) + "\t");
                    }
                }
                ps.print("\n");
                psInfectionTimes.print("\n");
            }

            first = false;
            Integer[] speciations = new Integer[hostsList.size() * hostsList.size()];
            Arrays.fill(speciations, 0);

            Double[] speciationTimes = new Double[hostsList.size()-1];
            Arrays.fill(speciationTimes, 0.0);

            if(!tree.getRoot().metaDataString.contains("host")){
                addHostMetadata(tree.getRoot());
            }

            List<String> hostsListNoUnsampled = new ArrayList<>(hostsList);
            hostsListNoUnsampled.remove("unsampled");


            for (Node leaf : tree.getExternalNodes()){
                fillTransmissions(leaf, speciations, hostsList);
                if(hostsListNoUnsampled.size()==0)
                    System.out.print("");
                fillInfectionTimes(leaf, speciationTimes, hostsListNoUnsampled);
            }
            for (int i = 0; i < hostsList.size() * hostsList.size(); i++) {
                if (skipIDs.contains(i))
                    continue;
                if (i == hostsList.size() * hostsList.size() - 2) {
                    ps.print(speciations[i]);
                } else {
                    ps.print(speciations[i] + "\t");
                }
            }

            for (int i=0; i<hostsListNoUnsampled.size();i++){
                if(i == hostsListNoUnsampled.size()-1){
                    psInfectionTimes.print(speciationTimes[i]);
                } else {
                    psInfectionTimes.print(speciationTimes[i]+"\t");
                }
            }

            ps.print("\n");
            psInfectionTimes.print("\n");
        }
        System.out.println("\nDone!");


    }

    private void fillTransmissions(Node leaf, Integer[] transmissions,
                                   List<String> hostsList){
        String recipient= leaf.getMetaData("host").toString().split("\\.")[0];
        Node child = leaf;
        Node parent = leaf.getParent();
        int nTransmissions = 0;

        if ( !parent.isFake() && child.metaDataString.contains("right")){
            nTransmissions +=1;
        }

        while(!parent.isRoot() &&
                (parent.getMetaData("host").toString().contains("unsampled") ||
                        Objects.equals(parent.getMetaData("host").toString().split("\\.")[0],
                                leaf.getMetaData("host").toString().split("\\.")[0]))){
            if(parent.getChildCount() == 1 || (!parent.isFake() && child.metaDataString.contains("right"))){
                nTransmissions +=1; // add 1 for every unobserved transmission (obtained by stochastic mapping)
            }
            child = parent;
            parent = parent.getParent();
        }
        String donor = parent.getMetaData("host").toString().split("\\.")[0];

        int idx = hostsList.indexOf(donor);
        nTransmissions = nTransmissions == 0 ? 1 : nTransmissions;
        idx = idx * hostsList.size() + hostsList.indexOf(recipient);
        transmissions[idx] = nTransmissions;
    }
    private void fillInfectionTimes(Node leaf, Double[] infectionTimes,
                                   List<String> hostsListNoUnsampled) {
        String recipient= leaf.getMetaData("host").toString().split("\\.")[0];
        Node child = leaf;
        Node parent = leaf.getParent();
        double infectionTime = Double.NaN;
        while (parent.isFake() || (parent.getChildCount()>1 && child.metaDataString.contains("donor") )){
            child = parent;
            parent = parent.getParent();

            if(parent == null){
                infectionTime=Double.NaN;
                break;
            }
            infectionTime = parent.getHeight();
        }

        if (parent != null){
            infectionTime = parent.getHeight();
        }

        int idx = hostsListNoUnsampled.indexOf(recipient);
        infectionTimes[idx] = infectionTime;
    }

    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(SpeciationAnalyser.SpeciationAnalyserOptions options) {

        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("sRanges speciation analyser");

        JLabel netFileLabel = new JLabel("Tree file:");
        JLabel outFileLabel = new JLabel("Output file:");

        JTextField netFilename = new JTextField(20);
        netFilename.setEditable(false);
        JButton netFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
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
                        .addComponent(netFileLabel)
                        .addComponent(outFileLabel))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(netFilename)
                        .addComponent(outFilename))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(netFileButton)
                        .addComponent(outFileButton))
        );

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(netFileLabel)
                        .addComponent(netFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(netFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton))
        );

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Run");
        runButton.addActionListener((e) -> dialog.setVisible(false));
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);


        JFileChooser netFileChooser = new JFileChooser();
        netFileButton.addActionListener(e -> {
            netFileChooser.setDialogTitle("Select TnT transmission tree file to remap");
            if (options == null)
                netFileChooser.setCurrentDirectory(options.inTransmissionTreeFile);
            else
                netFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = netFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.inTransmissionTreeFile = netFileChooser.getSelectedFile();
                netFilename.setText(netFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });


        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
            outFileChooser.setDialogTitle("Select output XML file name");
            outFileChooser.setCurrentDirectory(Objects.requireNonNull(options.inTransmissionTreeFile));

//            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
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
            public void write(int b) throws IOException {
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

    public static String helpMessage = "Indirect transmission analyser - outputs transmissions that are compatible with the tree topology.\n"
            + "\n"
            + "Usage: appstore IndirectTransmissionAnalyser [-help] |  logFile [outputFile]\n"
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
     * Retrieve ReMapTool options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, SpeciationAnalyser.SpeciationAnalyserOptions options) {
        int i=0;
        while (args.length > i && args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;
                case "-tree":
                    if (args.length<=i+1) {
                        printUsageAndError("-tree must be followed by a transmission tree file path.");
                    }

                    try {
                        options.inTransmissionTreeFile = new File(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing network file.");
                    }

                    i += 1;
                    break;

                case "-burnIn":
                    if (args.length<=i+1) {
                        printUsageAndError("-out must be followed by an output file path.");
                    }

                    try {
                        options.burnIn = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing output file path.");
                    }

                    i += 1;
                    break;

                case "-out":
                    if (args.length<=i+1) {
                        printUsageAndError("-out must be followed by an output file path.");
                    }

                    try {
                        options.outFile = new File(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing output file path.");
                    }

                    i += 1;
                    break;


                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }
    }

    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
        SpeciationAnalyser.SpeciationAnalyserOptions options = new SpeciationAnalyser.SpeciationAnalyserOptions();

        if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
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
            new SpeciationAnalyser(options);

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

    /**
     * Adds metadata on host occupying transmission tree edge.
     * Should be used on trees with full transmission histories.
     * That means either fully sampled trees or trees with stochastically mapped hidden transmission events.
     *
     * Relies on orientation metadata being added before
     * (public parent mathod makes sure of this).
     *
     * @param subTreeRoot the node at which to start
     */
    private void addHostMetadata(Node subTreeRoot){
//        Node subTreeRoot = getNode(subtreeRootNr);
        String metaData = subTreeRoot.metaDataString;

        if (subTreeRoot.isLeaf()) {
            subTreeRoot.setMetaData("host", Tools.removeLastSubstring("_",subTreeRoot.getID()));
            subTreeRoot.metaDataString = String.format("%s,%s=%s",
                    metaData, "host", Tools.removeLastSubstring("_",subTreeRoot.getID()));
            return;
        }
        if (subTreeRoot.isFake()){
            subTreeRoot.setMetaData("host", firstLeafId(subTreeRoot.getDirectAncestorChild()));
            subTreeRoot.metaDataString = String.format("%s,%s=%s",
                    metaData, "host", firstLeafId(subTreeRoot.getDirectAncestorChild()));
        }else {
            Node left = subTreeRoot.getLeft().getMetaData("orientation").equals("left") ?
                    subTreeRoot.getLeft() : subTreeRoot.getRight();
            subTreeRoot.setMetaData("host", firstLeafId(left));
            subTreeRoot.metaDataString = String.format("%s,%s=%s",
                    metaData, "host", firstLeafId(left));
        }
        addHostMetadata(subTreeRoot.getLeft());
        if(subTreeRoot.getChildCount()!=1){
            addHostMetadata(subTreeRoot.getRight());
        }
    }

    private String firstLeafId(Node n) {

        if (n.isLeaf()) {
            return n.getID().split("_")[0];
        } else {
            if (Tools.equalHeightWithPrecision(n, n.getChild(0))) {
                return Tools.removeLastSubstring("_", n.getChild(0).getID());
            }
            if (n.getChildCount() > 1 && Tools.equalHeightWithPrecision(n, n.getChild(1))) {
                return Tools.removeLastSubstring("_", n.getChild(1).getID());
            } else if (n.getChildCount() == 1) {
                return "unsampled";
            } else {
                Node left = n.getChild(0).metaDataString.contains("orientation=left") ? n.getChild(0)
                        : n.getChild(1);
                return firstLeafId(left);
            }

 /*           Node left = n.getChild(0).metaDataString.contains("orientation=donor") ? n.getChild(0)
                    : n.getChild(1);

            if (Tools.equalHeightWithPrecision(n, left)) {
                return left.getID().split("_")[0];
            } else if(n.getChildCount()>1){
                Node right = n.getChild(0).metaDataString.contains("orientation=donor") ? n.getChild(1)
                        : n.getChild(0);
                if(Tools.equalHeightWithPrecision(n, right)){
                    return right.getID().split("_")[0];
                }
            } else if (n.getChildCount()==1){
                return "unsampled";
            }

            return firstLeafId(left);*/
        }
    }
}
