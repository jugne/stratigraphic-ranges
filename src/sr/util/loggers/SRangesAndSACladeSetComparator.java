package sr.util.loggers;


import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.*;
import java.text.DecimalFormat;
import java.util.*;

import javax.imageio.ImageIO;

import beast.base.inference.parameter.RealParameter;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Runnable;
import beast.base.core.Log;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.FrequencySet;
import beast.base.util.Randomizer;

import static sr.util.Tools.*;

@Description("Match clades from two tree sets and print support for both sets so "
        + "they can be plotted in an X-Y plot")
public class SRangesAndSACladeSetComparator extends Runnable {
    final public Input<TreeFile> src1Input = new Input<>("SATree","source SA tree (set or MCC tree) file");
    final public Input<TreeFile> SATreeOffsetInput = new Input<>("SALog","SA log file containing offset for each tree in SATree");
    final public Input<TreeFile> src2Input = new Input<>("sRangesTree","source sRanges tree (set or MCC tree) file");

    public Input<String> sepStringInput = new Input<>("sep",
            "separator string for ranges",
            "_");

    public Input<Boolean> wholeRangeInput = new Input<>("wholeRange",
            "Does the clade have to include the whole range of the species? " +
                    "If true, the clade includes the whole range of the species, " +
                    "if false, the clade may include the whole range or just part of it. " +
                    "Default is false.",
            false);
    final public Input<OutFile> outputInput = new Input<>("out", "output file, or stdout if not specified",
            new OutFile("[[none]]"));
    final public Input<OutFile> svgOutputInput = new Input<>("svg", "svg output file. if not specified, no SVG output is produced.",
            new OutFile("[[none]]"));
    final public Input<OutFile> pngOutputInput = new Input<>("png", "png output file. if not specified, no PNG output is produced.",
            new OutFile("[[none]]"));
    final public Input<Integer> burnInPercentageInput = new Input<>("burnin", "percentage of trees to used as burn-in (and will be ignored)", 10);
    final public Input<Integer> thinningInput = new Input<>("thin", "thin out tree set. When thin=`n` only the first out of every n trees "
            + "is processed, and n-1 trees are skipped. This can be useful for large tree sets", 1);

    final public Input<Boolean> verboseInput = new Input<>("verbose", "print information about clades of interest, and if no output file is specified, all clade information", true);
    final public Input<Double> thresholdInput = new Input<>("threshold", "posterior support level of clades that will be ignored", 0.0);

    final public Input<String> taxonSetSAInput = new Input<>("taxonsetSA",
            "comma separated taxa in SA tree that are equivalent to taxons in sRanges tree. " +
                    "Use only if there are spelling differences or species renaming.",
            "");

    final public Input<String> taxonSetSRangesInput = new Input<>("taxonsetSRanges",
            "comma separated taxa in sRanges tree that are equivalent to taxons in SA tree. " +
                    "Use only if there are spelling differences or species renaming.",
            "");

    final public Input<String> compareTaxonInput = new Input<>("compareTaxon",
            "comma separated taxa in sRanges tree that the figure and log file should be produced for " +
                    "Use only if not all taxa in the sRanges tree should be compared.");

    final public Input<Double> scalingInput = new Input<>("scaling",
            "scaling factor for clade ages " +
                    "If not specified the highest tree in both trees is used as scaling factor.");

    private double n;
    private boolean verbose;
    private double threshold = 0;

    final String header = "<svg version=\"1.2\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" class=\"graph\" aria-labelledby=\"title\" role=\"img\" height=\"1200\">\n" +
            "<g class=\"grid x-grid\" id=\"xGrid\">\n" +
            "  <line x1=\"90\" x2=\"90\" y1=\"10\" y2=\"1010\" style=\"stroke:#000;stroke-width:2\"></line>\n" +
            "</g>\n" +
            "<g class=\"grid y-grid\" id=\"yGrid\">\n" +
            "  <line x1=\"90\" x2=\"1090\" y1=\"1010\" y2=\"1010\" style=\"stroke:#000;stroke-width:2\"></line>\n" +
            "</g>\n" +
            "<line x1=\"90\" x2=\"1090\" y1=\"1010\" y2=\"10\" style=\"stroke:#000;stroke-width:2\"></line>\n" +
            "<line x1=\"90\" x2=\"840\" y1=\"760\" y2=\"10\" style=\"stroke:#00f;stroke-width:1\"></line>\n" +
            "<line x1=\"340\" x2=\"1090\" y1=\"1010\" y2=\"260\" style=\"stroke:#00f;stroke-width:1\"></line>\n" +
            "  <g class=\"labels x-labels\">\n" +
            "  <text x=\"90\" y=\"1030\">0.0</text>\n" +
            "  <text x=\"290\" y=\"1030\">0.2</text>\n" +
            "  <text x=\"490\" y=\"1030\">0.4</text>\n" +
            "  <text x=\"690\" y=\"1030\">0.6</text>\n" +
            "  <text x=\"890\" y=\"1030\">0.8</text>\n" +
            "  <text x=\"1090\" y=\"1030\">1.0</text>\n" +
            "  <text x=\"520\" y=\"1040\" class=\"label-title\">file1</text>\n" +
            "</g>\n" +
            "<g class=\"labels y-labels\">\n" +
            "  <text x=\"60\" y=\"15\">1.0</text>\n" +
            "  <text x=\"60\" y=\"215\">0.8</text>\n" +
            "  <text x=\"60\" y=\"415\">0.6</text>\n" +
            "  <text x=\"60\" y=\"615\">0.4</text>\n" +
            "  <text x=\"60\" y=\"815\">0.2</text>\n" +
            "  <text x=\"60\" y=\"1015\">0.0</text>\n" +
            "  <text x=\"40\" y=\"540\" class=\"label-title\" transform=\"rotate(90,40,540)\">file2</text>\n" +
            "</g>\n" +
            "<g class=\"data\" data-setname=\"Our first data set\">\n";
    //			"  <circle cx=\"90\" cy=\"192\" data-value=\"7.2\" r=\"4\"></circle>\n" +
//			"  <circle cx=\"240\" cy=\"141\" data-value=\"8.1\" r=\"4\"></circle>\n" +
//			"  <circle cx=\"388\" cy=\"179\" data-value=\"7.7\" r=\"4\"></circle>\n" +
//			"  <circle cx=\"531\" cy=\"200\" data-value=\"6.8\" r=\"4\"></circle>\n" +
//			"  <circle cx=\"677\" cy=\"104\" data-value=\"6.7\" r=\"4\"></circle>\n" +
    final String footer =	"</g>\n" +
            "</svg>\n" +
            "\n";


    @Override
    public void initAndValidate() {
    }


    double maxHeight = 0.0;
    int problemCount = 0;
    int interestCount = 0;
    int inconsistentHeightIntervals = 0;
    double meanHeightsDifference = 0;

    HashMap<String, String> SAtoSRangesTaxonMap = new HashMap<>();


    @Override
    public void run() throws Exception {
        if (scalingInput.get()!=null){
            maxHeight = scalingInput.get();
        }
        if ((!taxonSetSAInput.get().isEmpty() && taxonSetSRangesInput.get().isEmpty()) ||
                (taxonSetSAInput.get().isEmpty() && !taxonSetSRangesInput.get().isEmpty())) {
            throw new IllegalArgumentException("taxonSetSA and taxonSetSRanges must be both specified or both empty");
        }

        if (!taxonSetSAInput.get().isEmpty() && !taxonSetSRangesInput.get().isEmpty()){
            List<String> keys = new ArrayList<String>(Arrays.asList(taxonSetSAInput.get().split(",")));
            List<String> values = new ArrayList<String>(Arrays.asList(taxonSetSRangesInput.get().split(",")));
            if (keys.size() != values.size()) {
                throw new IllegalArgumentException("taxonSetSA and taxonSetSRanges must have the same number of elements");
            }
            for (int i = 0; i < keys.size(); i++) {
                SAtoSRangesTaxonMap.put(keys.get(i), values.get(i));
            }
        }

        verbose = verboseInput.get();
        threshold = thresholdInput.get();
        CladeSetWithHeights cladeSet1 = getCladeSet(src1Input.get().getPath(), false);
        double n1 = n;
        CladeSetWithHeights cladeSet2 = getCladeSet(src2Input.get().getPath(), true);
        double n2 = n;
        process(src1Input.get(), src2Input.get(), "", cladeSet1, n1, cladeSet2, n2);
    }

    void process(TreeFile tree1, TreeFile tree2, String suffix,
                 CladeSetWithHeights cladeSet1, double n1,
                 CladeSetWithHeights cladeSet2, double n2
    )  throws Exception {
        PrintStream out = System.out;
        if (outputInput.get() != null && !outputInput.get().getName().equals("[[none]]")) {
            String str = normalise(outputInput.get().getPath(), suffix);
            Log.warning("Writing to file " + str);
            out = new PrintStream(str);
        }
        PrintStream svg = null;
        if (svgOutputInput.get() != null && !svgOutputInput.get().getName().equals("[[none]]")) {
            String str = normalise(svgOutputInput.get().getPath(), suffix);
            Log.warning("Writing to file " + str);
            svg = new PrintStream(str);
            svg.println(header.replaceAll("file1", tree1.getPath()).replaceAll("file2", tree2.getPath()));
        }

        Graphics2D g = null;
        BufferedImage bi = null;
        if (pngOutputInput.get() != null && !pngOutputInput.get().getName().equals("[[none]]")) {
            bi = new BufferedImage(1200, 1200, BufferedImage.TYPE_INT_ARGB);
            g = (Graphics2D) bi.getGraphics();
            initPNG(g, tree1, tree2);
        }

//		CladeSetWithHeights cladeSet1 = getCladeSet(tree1.getPath());
//		double n1 = n;

//		CladeSetWithHeights cladeSet2 = getCladeSet(tree2.getPath());
//		double n2 = n;

        // create map of clades to support values in set1
        Map<String, Double> cladeMap = new LinkedHashMap<>();
        Map<String, Integer> cladeToIndexMap = new LinkedHashMap<>();
        Map<String, Double> cladeHeightMap = new LinkedHashMap<>();
        for (int i = 0; i < cladeSet1.getCladeCount(); i++) {
            String clade = cladeSet1.getClade(i, false, SAtoSRangesTaxonMap);
            BitSet bitset = cladeSet1.get(i);
            int support = cladeSet1.getFrequency(i);
            if (cladeSet1 instanceof SummaryCladeSetWithHeights) {
                cladeMap.put(clade, ((SummaryCladeSetWithHeights)cladeSet1).posteriors.get(bitset));
            } else {
                cladeMap.put(clade, support/ n1);
            }
            cladeHeightMap.put(clade, cladeSet1.getMeanNodeHeight(i));
            cladeToIndexMap.put(clade, i);
        }

        // process clades in set2
        double maxDiff = 0, meanDiff=0, meanDiff2=0;
        int meanDiffCount = 0, meanDiff2Count = 0;
        double [] hist = new double[40];
        problemCount = 0;
        interestCount = 0;
        meanHeightsDifference = 0;
        inconsistentHeightIntervals = 0;

        double clade1ThresholdSupport = getThresholdSupport(cladeSet1, n1);
        double clade2ThresholdSupport = getThresholdSupport(cladeSet2, n2);

        for (int i = 0; i < cladeSet2.getCladeCount(); i++) {
            String clade = cladeSet2.getClade(i, true, SAtoSRangesTaxonMap);
            double support2;
            if (cladeSet2 instanceof SummaryCladeSetWithHeights) {
                BitSet bitset = cladeSet2.get(i);
                support2 = ((SummaryCladeSetWithHeights)cladeSet2).posteriors.get(bitset);
            } else {
                int support = cladeSet2.getFrequency(i);
                support2 = support/n2;
            }
            double h2 = cladeSet2.getMeanNodeHeight(i);
            if (cladeMap.containsKey(clade)) {
                // clade is also in set1
                double h1 = cladeHeightMap.get(clade);
                double [] heights1 = cladeSet1.nodeHeights.get(cladeSet1.get(cladeToIndexMap.get(clade)));
                Arrays.sort(heights1);
                double lo1 = heights1[(int)(heights1.length * 0.025)];
                double hi1 = heights1[(int)(heights1.length * 0.975)];
                double max1 = heights1[heights1.length-1];

                double [] heights2 = cladeSet2.nodeHeights.get(cladeSet2.get(i));
                Arrays.sort(heights2);
                double lo2 = heights2[(int)(heights2.length * 0.025)];
                double hi2 = heights2[(int)(heights2.length * 0.975)];
                double max2 = heights2[heights2.length-1];

                double support1 = cladeMap.get(clade);
                output(out, svg, clade,support1, support2, g, h1, h2,
                        lo1, lo2, hi1, hi2, max1, max2);
                // System.out.println((h1 - h2) + " " + (100 * (h1 - h2) / h1));

                maxDiff = Math.max(maxDiff, Math.abs(cladeMap.get(clade) - support2));
                if (support1 > 0.01 && support2 > 0.01) {
                    meanDiff += Math.abs(cladeMap.get(clade) - support2);
                    meanHeightsDifference += Math.abs(hi1-hi2)/(h1+h2)/2.0;
                    meanDiffCount++;
                    if (lo1 > hi2 || lo2 > hi1) {
                        inconsistentHeightIntervals++;
                        System.err.println("Inconsistent clade height found for clade: "  + clade.replaceAll(" ", ""));
                    }
                }
                if (support2 >= clade2ThresholdSupport || support1 >= clade1ThresholdSupport) {
                    meanDiff2 += (support1 - support2) * (support1 - support2);
                    meanDiff2Count++;
                }
                cladeMap.remove(clade);

                // record difference in 95%HPD (if support > 1% in both clade sets)
                if (support1 > 0.01 && support2 > 0.01) {
                    if ((hi1-lo1) < (hi2-lo2)) {
                        double w = (hist.length/2) * (hi1-lo1)/(hi2-lo2);
                        if (w < 0) {
                            w = 0;
                        }
                        hist[(int)w] += support1 + support2;
                    } else {
                        double w = hist.length - 1 - (hist.length/2) * (hi2-lo2)/(hi1-lo1);
                        if (w >= hist.length) {
                            w = hist.length - 1;
                        }
                        hist[(int)w] += support1 + support2;
                    }
                }
            } else {
                // clade is not in set1
                double [] heights2 = cladeSet2.nodeHeights.get(cladeSet2.get(i));
                Arrays.sort(heights2);
                double lo2 = heights2[(int)(heights2.length * 0.025)];
                double hi2 = heights2[(int)(heights2.length * 0.975)];
                double max2 = heights2[heights2.length-1];
                output(out, svg, clade, 0.0, support2, g, 0, h2, 0, lo2, 0, hi2, 0, max2);
                maxDiff = Math.max(maxDiff, support2);
                if (support2> 0.01) {
                    meanDiff += support2;
                    meanDiffCount++;
                }
                if (support2>clade2ThresholdSupport) {
                    meanDiff2 += support2*support2;
                    meanDiff2Count += 1;
                }
            }
        }

        // process left-overs of clades in set1 that are not in set2
        for (String clade : cladeMap.keySet()) {
            double h1 = cladeHeightMap.get(clade);
            double [] heights1 = cladeSet1.nodeHeights.get(cladeSet1.get(cladeToIndexMap.get(clade)));
            Arrays.sort(heights1);
            double lo1 = heights1[(int)(heights1.length * 0.025)];
            double hi1 = heights1[(int)(heights1.length * 0.975)];
            double max1 = heights1[heights1.length-1];
            output(out, svg, clade, cladeMap.get(clade), 0.0, g, h1, 0.0, lo1, 0, hi1, 0, max1, 0.0);
            double s = cladeMap.get(clade);
            maxDiff = Math.max(maxDiff, s);
            if (s> 0.01) {
                meanDiff += s;
                meanDiffCount++;
            }
            if (s>clade1ThresholdSupport) {
                meanDiff2 += s*s;
                meanDiff2Count += 1;
            }
        }

        double sqrtMeanSumSquared = Math.sqrt(meanDiff2/meanDiff2Count)*100;
        String measure3String = "Sqrt(Mean squared difference) in clade support for most probable N-2 clades: ";

        final DecimalFormat formatter = new DecimalFormat("#.##");
        if (svg != null) {
            svg.println("<text x='110' y='25'>Max difference in clade support: " + formatter.format(maxDiff * 100)+ "%</text>");
            svg.println("<text x='110' y='45'>Mean difference in clade support (when sum over 1%): " + formatter.format(meanDiff/meanDiffCount * 100)+ "%</text>");
            svg.println("<text x='110' y='65'>" + measure3String + formatter.format(sqrtMeanSumSquared)+ "%</text>");
            svg.println("<text x='110' y='85'>" + interestCount + " clades >25% difference "+ problemCount + " problematic</text>");
            svg.println("<text x='110' y='105'>" + inconsistentHeightIntervals + " inconsistent height intervals " + formatter.format(100.0*meanHeightsDifference/meanDiffCount) + " average % mean height diff</text>");
            svg.println(footer);
        }
        if (bi != null) {
            // draw histogram of 95%HPD interval fractions
            double max = 0;
            for (double d : hist) {
                max = Math.max(max, d);
            }
            int width = 10, height = 90;
            g.setComposite(AlphaComposite.SrcOver.derive(1.0f));
            g.drawRect(100, 0, width*hist.length, height);
            for (int i = 0; i < hist.length; i++) {
                g.drawRect(100 + i * width, height-(int)(hist[i] * height / max), width, (int)(hist[i] * height / max));
            }

            g.setFont(new Font("Arial", Font.PLAIN, 18));
            g.setColor(Color.black);
            g.drawString("Max difference in clade support: " + formatter.format(maxDiff * 100)+ "%", 510, 15);
            g.drawString("Mean difference in clade support (when sum over 1%): " + formatter.format(meanDiff/meanDiffCount * 100)+ "%", 510, 35);
            g.drawString(measure3String + formatter.format(sqrtMeanSumSquared)+ "%", 510, 55);
            g.drawString(interestCount + " clades >25% difference "+ problemCount + " problematic", 510, 75);
            g.drawString(inconsistentHeightIntervals + " inconsistent height intervals " + formatter.format(100.0*meanHeightsDifference/meanDiffCount) + " average % mean height diff", 510, 95);

            String str = normalise(pngOutputInput.get().getPath(), suffix);
            Log.warning("Writing to file " + str);
            ImageIO.write(bi, "png", new File(str));
        }
        Log.info("Maximum difference in clade support: " + maxDiff);
        Log.info("Mean difference in clade support (when sum over 1%): " + meanDiff/meanDiffCount);
        Log.info(measure3String + sqrtMeanSumSquared);
        Log.info(interestCount + " clades >25% difference "+ problemCount + " problematic");
        Log.info(inconsistentHeightIntervals + " inconsistent height intervals " + formatter.format(100.0*meanHeightsDifference/meanDiffCount) + " average % mean height diff");
        Log.info.println("Done");
    }

    private double getThresholdSupport(CladeSetWithHeights cladeSet, double totalCount) {
        double[] supportValues = new double[cladeSet.getCladeCount()];
        for (int i = 0; i < cladeSet.getCladeCount(); i++) {
            supportValues[i] = getSupport(cladeSet, i, totalCount);
        }
        Arrays.sort(supportValues);

        return supportValues[supportValues.length-(cladeSet.taxonNames.size()-2)];
    }

    private double getSupport(CladeSetWithHeights cladeSet, int cladeIndex, double totalCount) {
        double support = 0.0;
        if (cladeSet instanceof SummaryCladeSetWithHeights) {
            BitSet bitset = cladeSet.get(cladeIndex);
            support = ((SummaryCladeSetWithHeights)cladeSet).posteriors.get(bitset);
        } else {
            support = cladeSet.getFrequency(cladeIndex) / totalCount;
        }
        return support;
    }

    private String normalise(String str, String suffix) {
        if (str.contains(".")) {
            int k = str.lastIndexOf('.');
            str = str.substring(0, k) + suffix + str.substring(k);
        } else {
            str += suffix;
        }
        return str;
    }

    private void initPNG(Graphics2D g, TreeFile tree1, TreeFile tree2) {
        g.setColor(Color.white);
        g.fillRect(0, 0, 1200, 1200);

        g.setColor(Color.black);
        g.drawRect(100, 100, 1000, 1000);

        // diagonals
        int h = 1200;
        g.drawLine(100, h-100, 1100, h-1100);
        g.setColor(Color.blue);
        g.drawLine(100, h-300,  900, h-1100);
        g.drawLine(300, h-100, 1100, h- 900);

        g.drawString("0.0", 100, 1130);
        g.drawString("0.2", 300, 1130);
        g.drawString("0.4", 500, 1130);
        g.drawString("0.6", 700, 1130);
        g.drawString("0.8", 900, 1130);
        g.drawString("1.0",1100, 1130);

        g.drawString("0.0", 60, h-100);
        g.drawString("0.2", 60, h-300);
        g.drawString("0.4", 60, h-500);
        g.drawString("0.6", 60, h-700);
        g.drawString("0.8", 60, h-900);
        g.drawString("1.0", 60, h-1100);

        g.setColor(Color.black);
        g.setFont(new Font("Arial",Font.PLAIN, 50));
        g.drawString("sampled ancestor tree", 120, 1170);

        AffineTransform orig = g.getTransform();
        g.rotate(-Math.PI/2);
        g.drawString("stratigraphic range tree", -880, 40);
        g.setTransform(orig);

        g.setColor(Color.red);
        g.setComposite(AlphaComposite.SrcOver.derive(0.25f));
    }

    private void output(PrintStream out, PrintStream svg, String clade, Double support1, double support2, Graphics2D g, double h1, double h2,
                        double lo1, double lo2, double hi1, double hi2, double max1, double max2) {
        if (compareTaxonInput.get()!=null){
            List<String> cladeList = Arrays.asList(clade
                    .replaceAll(" ", "")
                    .replaceAll("\\{", "")
                    .replaceAll("}", "")
                    .split(","));
            List<String> compareOnly = Arrays.asList(compareTaxonInput.get().split(","));
            if (!compareOnly.containsAll(cladeList)){
                return;
            }
        }
        if (verbose || System.out != out) {
            out.println(clade.replaceAll(" ", "") + " " + support1 + " " + support2 + " " + h1 +
                    " " + lo1 + " " + hi1 + " " + h2 + " " + lo2 + " " + hi2);
        }
//		if ((support1 < 0.1 && support2 > 0.9) ||
//			(support2 < 0.1 && support1 > 0.9)) {
        if (Math.abs(support1 - support2) > 0.9) {
            if (verbose) {
                Log.warning("Problem clade: " + clade.replaceAll(" ", "") + " " + support1 + " " + support2);
            }
            problemCount++;
        }

        if (Math.abs(support1 - support2) > 0.25) {
            if (verbose) {
                Log.warning("Clade of interest (>25% difference): " + clade.replaceAll(" ", "") + " " + support1 + " " + support2);
            }
            interestCount++;
        }

        if (svg != null && threshold <= support1 + support2) {
            svg.println("  <circle style=\"opacity:0.25;fill:#a00000\" cx=\""+ (90 +1000* support1 + Randomizer.nextInt(10) - 5) +
                    "\" cy=\""+ (10 + 1000 - 1000 * support2 + Randomizer.nextInt(10) - 5) +"\" "
                    + "data-value=\"7.2\" r=\"" + (support1 + support2) * 10 + "\"></circle>");

            if ((support1 + support2) > 0.1) {
                int d1 = 90;
                int d2 = 10;
                int x1 = (int)(d1 + 1000.0 * lo1 / maxHeight);
                int y1 = (int)(1000+d2 - 1000.0 * h2/ maxHeight);
                int x2 = (int)(d1 + 1000.0 * hi1 / maxHeight);
                int y2 = (int)(1000 +d2- 1000.0 * h2/ maxHeight);
                svg.println("<line stroke=\"#0000a040\" x1=\""+x1+"\" y1=\""+y1+"\" x2=\""+x2+"\" y2=\""+y2+"\"/>");
                x1 = (int)(d1 + 1000.0 * h1 / maxHeight);
                y1 = (int)(1000 + d2 - 1000.0 * lo2/ maxHeight);
                x2 = (int)(d1 + 1000.0 * h1 / maxHeight);
                y2 = (int)(1000 + d2 - 1000.0 * hi2/ maxHeight);
                svg.println("<line stroke=\"#0000a040\" x1=\""+x1+"\" y1=\""+y1+"\" x2=\""+x2+"\" y2=\""+y2+"\"/>");
            }
        }

        if (g != null) {
            double x = (100 + 1000 * support1 + Randomizer.nextInt(10) - 5);
            double y = (     1100 - 1000 * support2 + Randomizer.nextInt(10) - 5);
            double r = 1+(support1 + support2) * 10;
            g.setColor(Color.red);
            g.setComposite(AlphaComposite.SrcOver.derive(0.25f));
            g.fillOval((int)(x-r/2), (int)(y-r/2), (int) r, (int) r);

            g.setColor(Color.blue);
            float alpha = (float)(0.1 + ((support1 + support2)/2.0)*0.9);
            g.setComposite(AlphaComposite.SrcOver.derive(alpha));
            x = 100 + 1000.0 * h1 / maxHeight;
            y = 1100 - 1000.0 * h2/ maxHeight;
            r = 3 + Math.max(support1, support2) * 13;
            g.fillOval((int)(x-r/2), (int)(y-r/2), (int) r, (int) r);


            if ((support1 + support2) > 0.1) {
                g.setComposite(AlphaComposite.SrcOver.derive(alpha * alpha));
                int x1 = (int)(100 + 1000.0 * lo1 / maxHeight);
                int y1 = (int)(1100 - 1000.0 * h2/ maxHeight);
                int x2 = (int)(100 + 1000.0 * hi1 / maxHeight);
                int y2 = (int)(1100 - 1000.0 * h2/ maxHeight);
                g.drawLine(x1, y1, x2, y2);
                x1 = (int)(100 + 1000.0 * h1 / maxHeight);
                y1 = (int)(1100 - 1000.0 * lo2/ maxHeight);
                x2 = (int)(100 + 1000.0 * h1 / maxHeight);
                y2 = (int)(1100 - 1000.0 * hi2/ maxHeight);
                g.drawLine(x1, y1, x2, y2);
            }

        }
    }

    private CladeSetWithHeights getCladeSet(String path, boolean rangesTree) throws IOException {
        Log.warning("Processing " + path);
        MemoryFriendlyTreeSet srcTreeSet = new MemoryFriendlyTreeSet(path, burnInPercentageInput.get());
        BufferedReader offsetFile = null;
        int offsetIndex = 0;
        if (!rangesTree && SATreeOffsetInput.get()!=null){
            offsetFile = new BufferedReader(new FileReader(SATreeOffsetInput.get()));
            String line = offsetFile.readLine();
            int n = -1;
            while(line != null && !line.startsWith("#")){
                n+=1;
                line = offsetFile.readLine();
            }
            offsetFile = new BufferedReader(new FileReader(SATreeOffsetInput.get()));
            line = offsetFile.readLine();
            while(line.startsWith("#"))
                line = offsetFile.readLine();
            String[] params = line.split("\t");
            offsetIndex = findIndexOf(params, "offset");
            for (int i=0; i<Math.max(0, (burnInPercentageInput.get() * n)/100); i++){
                line = offsetFile.readLine();
            }
            line = offsetFile.readLine();
            params = line.split("\t");
            offset = Double.parseDouble(params[offsetIndex]);
        }
        else {
            offset =0.;
        }

        srcTreeSet.reset();
        Tree tree = srcTreeSet.next();
        CladeSetWithHeights cladeSet1 = new CladeSetWithHeights(tree, rangesTree);
        n = 1;
        int thin = thinningInput.get();
        if (scalingInput.get()==null) {
            maxHeight = Math.max(maxHeight, tree.getRoot().getHeight());
        }

        while (srcTreeSet.hasNext()) {
            // System.out.println(n);
            if (!rangesTree && SATreeOffsetInput.get()!=null){
                String line = offsetFile.readLine();
                String[] params = line.split("\t");
                offset = Double.parseDouble(params[offsetIndex]);
            } else{
                offset = 0.;
            }
            tree = srcTreeSet.next();
            cladeSet1.add(tree, rangesTree);
            n++;
            if (scalingInput.get()==null) {
                maxHeight = Math.max(maxHeight, tree.getRoot().getHeight());
            }
            int j = 1;
            while (j < thin && srcTreeSet.hasNext()) {
                tree = srcTreeSet.next();
                j++;
            }
        }

        if (n==1) {
            // might be a summary tree
            cladeSet1 = new SummaryCladeSetWithHeights(tree, rangesTree);
        }
        return cladeSet1;
    }

    public class SummaryCladeSetWithHeights extends CladeSetWithHeights {

        public SummaryCladeSetWithHeights(Tree tree, boolean rangesTree) {
            this(tree, tree.getTaxonset(), rangesTree);
        }

        public SummaryCladeSetWithHeights(Tree tree, TaxonSet taxonSet, boolean rangesTree) {
            this.taxonSet = taxonSet;
            this.taxonNames = this.taxonSet.getTaxaNames();
            if (rangesTree)
                this.taxonNames = taxaNamesToSpeciesNames(this.taxonNames, sepStringInput.get());
            add(tree, rangesTree);
        }

        void addClades(Node node, BitSet bits) {

            if (node.isLeaf()) {
                if (taxonSet != null) {

                    int index =  taxonSet.getTaxonIndex(node.getID());
                    bits.set(index);
                } else {
                    bits.set(node.getNr());
                }
            } else {

                BitSet bits2 = new BitSet();
                for (Node child : node.getChildren()) {
                    addClades(child, bits2);
                }

                add(bits2, 1);
                Double [] heightHPD = (Double []) node.getMetaData("height_95%_HPD");
                double lo = node.getHeight(), hi = node.getHeight();
                if (heightHPD != null) {
                    lo = heightHPD[0];
                    hi = heightHPD[1];
                }
                addNodeHeight(bits2, node.getHeight(), lo, hi); // TODO ?= tree.getNodeHeight(node)
                double posterior = 1;
                Double posterior_ = (Double) node.getMetaData("posterior");
                if (posterior_ != null) {
                    posterior = posterior_;
                }
                posteriors.put(bits2, posterior);

                if (bits != null) {
                    bits.or(bits2);
                }
            }
        }

        void addNodeHeight(BitSet bits, double height, double lo, double hi) {
            totalNodeHeight.put(bits, (getTotalNodeHeight(bits) + (height+offset)));
            nodeHeights.put(bits, new double[]{lo+offset, height+offset, hi+offset});
        }

        Map<BitSet, Double> posteriors = new HashMap<>();
    }

    public class CladeSetWithHeights extends FrequencySet<BitSet> {
        //
        // Public stuff
        //

        public CladeSetWithHeights() {}

        /**
         * @param tree
         */
        public CladeSetWithHeights(Tree tree, boolean rangesTree) {
            this(tree, tree.getTaxonset(), rangesTree);
        }

        /**
         * @param taxonSet  a set of taxa used to label the tips
         */
        public CladeSetWithHeights(Tree tree, TaxonSet taxonSet, boolean rangesTree) {
            this.taxonSet = taxonSet;
            this.taxonNames = taxonSet.getTaxaNames();
            if (rangesTree)
                this.taxonNames = taxaNamesToSpeciesNames(this.taxonNames, sepStringInput.get());
            add(tree, rangesTree);
        }

        /** get number of unique clades */
        public int getCladeCount()
        {
            return size();
        }

        /** get clade bit set */
        public String getClade(int index, boolean rangesTree, HashMap<String,String> taxonMap) {
            BitSet bits = get(index);

            StringBuffer buffer = new StringBuffer("{");
            boolean first = true;
            for (String taxonId : getTaxaSet(bits, rangesTree, taxonMap)) {
                if (!first) {
                    buffer.append(", ");
                } else {
                    first = false;
                }
                buffer.append(taxonId);
            }
            buffer.append("}");
            return buffer.toString();
        }

        private SortedSet<String> getTaxaSet(BitSet bits, boolean rangesTree, HashMap<String,String> taxonMap){

            SortedSet<String> taxaSet = new TreeSet<>();

            for (int i = 0; i < bits.length(); i++) {
                if (bits.get(i)) {
                    if (rangesTree)
                        taxaSet.add((String) taxonNames.toArray()[i]); //TODO ?= taxonList.getTaxonId(i)
                    else if (!taxonMap.isEmpty() && taxonMap.containsKey(taxonSet.asStringList().get(i)))
                        taxaSet.add(taxonMap.get(taxonSet.asStringList().get(i)));
                    else
                        taxaSet.add(taxonSet.asStringList().get(i));
                }
            }
            return taxaSet;
        }

        /** get clade frequency */
        int getCladeFrequency(int index)
        {
            return getFrequency(index);
        }


        public void add(Tree tree, Boolean rangesTree) {
            if (taxonSet == null) {
                taxonSet = tree.getTaxonset();
                taxonNames = taxonSet.getTaxaNames();
                if (rangesTree)
                    taxonNames = taxaNamesToSpeciesNames(taxonNames, sepStringInput.get());
            }

            totalTrees += 1;

            // Recurse over the tree and add all the clades (or increment their
            // frequency if already present). The root clade is not added.
            addClades(tree.getRoot(), null, rangesTree);
            addedForThisTree = new HashSet<>();
        }

        void addClades(Node node, BitSet bits, Boolean rangesTree) {

            if (node.isLeaf()) {
                if (taxonSet != null) {
                    String id = node.getID();
                    int index = taxonSet.getTaxonIndex(id);
                    if (wholeRangeInput.get() && id.contains("last"))
                        return;
                    if (rangesTree){
                        id = removeLastSubstring(sepStringInput.get(), id);
                        index = getSpeciesNumber(this.taxonNames, id);
                    }


                    if (index == -1) {
                        throw new RuntimeException("Taxon " + id + " not found in taxonset");
                    }
                    bits.set(index);
                } else {
                    bits.set(node.getNr());
                }
            } else {

                BitSet bits2 = new BitSet();
                for (Node child : node.getChildren()) {
                    addClades(child, bits2, rangesTree);
                }

                if (bits2.cardinality()==1 || this.addedForThisTree.contains(bits2)){ //&& node.getHeight()<this.nodeHeights.get(bits2)[this.nodeHeights.get(bits2).length-1]) {
//                    totalNodeHeight.put(bits, (getTotalNodeHeight(bits) + height));
//                    this.nodeHeights.get(bits2)[this.nodeHeights.get(bits2).length-1] = node.getHeight();
//                    System.out.println("Warning: clade " + bits2 + " already added for this tree");
                } else {
                    add(bits2, 1);
                    addNodeHeight(bits2, node.getHeight());
                    this.addedForThisTree.add(bits2);
                }

                if (bits != null) {
                    bits.or(bits2);
                }
            }
        }

        public double getMeanNodeHeight(int i) {
            BitSet bits = get(i);

            return getTotalNodeHeight(bits) / getFrequency(i);
        }

        public double getMaxNodeHeight(int i) {
            BitSet bits = get(i);
            double[] d = nodeHeights.get(bits);
            Arrays.sort(d);
            return d[d.length-1];
        }

        double getTotalNodeHeight(BitSet bits) {
            Double tnh = totalNodeHeight.get(bits);
            if (tnh == null) return 0.0;
            return tnh;
        }


        private void addNodeHeight(BitSet bits, double height) {
            totalNodeHeight.put(bits, (getTotalNodeHeight(bits) + height));
            if (!nodeHeights.containsKey(bits)) {
                nodeHeights.put(bits, new double[]{height+offset});
            } else {
                double [] heights = nodeHeights.get(bits);
                double [] newHeights = new double[heights.length + 1];
                System.arraycopy(heights, 0, newHeights, 0, heights.length);
                newHeights[heights.length] = height+offset;
                nodeHeights.put(bits, newHeights);
            }
        }

        // Generifying found that this code was buggy. Luckily it is not used anymore.

//	    /** adds all the clades in the CladeSet */
//	    public void add(CladeSet cladeSet)
//	    {
//	        for (int i = 0, n = cladeSet.getCladeCount(); i < n; i++) {
//	            add(cladeSet.getClade(i), cladeSet.getCladeFrequency(i));
//	        }
//	    }

        private BitSet annotate(Tree tree, Node node, String freqAttrName) {
            BitSet b = null;
            if (node.isLeaf()) {
                int index;
                if (taxonSet != null) {
                    index = taxonSet.getTaxonIndex(node.getID());
                } else {
                    index = node.getNr();
                }
                b = new BitSet(tree.getLeafNodeCount());
                b.set(index);

            } else {

                for (Node child : node.getChildren()) {
                    BitSet b1 = annotate(tree, child, freqAttrName);
                    if( child.isRoot() ) {
                        b = b1;
                    } else {
                        b.or(b1);
                    }
                }
                final int total = getFrequency(b);
                if( total >= 0 ) {
                    node.setMetaData(freqAttrName, total / (double)totalTrees );
                }
            }
            return b;
        }

        /**
         * Annotate clades of tree with posterior probability
         * @param tree
         * @param freqAttrName name of attribute to set per node
         * @return sum(log(all clades probability))
         */
        public double annotate(Tree tree, String freqAttrName) {
            annotate(tree, tree.getRoot(), freqAttrName);

            double logClade = 0.0;
            for(Node internalNode : tree.getInternalNodes()) {
                final double f = (Double) internalNode.getMetaData(freqAttrName);
                logClade += Math.log(f);
            }
            return logClade;
        }

        public boolean hasClade(int index, Tree tree) {
            BitSet bits = get(index);

            Node[] mrca = new Node[1];
            findClade(bits, tree.getRoot(), mrca);

            return (mrca[0] != null);
        }

        private int findClade(BitSet bitSet, Node node, Node[] cladeMRCA) {

            if (node.isLeaf()) {

                if (taxonSet != null) {
                    int index = taxonSet.getTaxonIndex(node.getID());
                    if (bitSet.get(index)) return 1;
                } else {
                    if (bitSet.get(node.getNr())) return 1;
                }
                return -1;
            } else {
                int count = 0;
                for (Node child : node.getChildren()) {
                    int childCount = findClade(bitSet, child, cladeMRCA);

                    if (childCount != -1 && count != -1) {
                        count += childCount;
                    } else count = -1;
                }

                if (count == bitSet.cardinality()) cladeMRCA[0] = node;

                return count;
            }
        }

        //
        // Private stuff
        //
        TaxonSet taxonSet = null;
        Set<String> taxonNames = null;
        Map<BitSet, Double> totalNodeHeight = new HashMap<>();
        Map<BitSet, double[]> nodeHeights = new HashMap<>();

        Set<BitSet> addedForThisTree = new HashSet<>();
        int totalTrees = 0;
    }

    double offset = 0.;
    public static void main(String[] args) throws Exception {
        new Application(new SRangesAndSACladeSetComparator(), "Clade Set Comparator", args);

    }

    // Helper method to find the index of a specific value in an array
    public static int findIndexOf(String[] array, String value) {
        for (int i = 0; i < array.length; i++) {
            if (array[i].equals(value)) {
                return i;
            }
        }
        return -1;  // Value not found
    }


}