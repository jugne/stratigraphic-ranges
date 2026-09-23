package evolution.tree;

import beast.base.core.Input;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.distribution.Prior;
import beast.base.inference.distribution.Uniform;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import sa.evolution.operators.SAScaleOperator;
import sr.evolution.operators.LeftRightChildSwap;
import sr.evolution.operators.SRLeafToSampledAncestorJump;
import sr.evolution.operators.SRWilsonBalding;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.RandomSRangeTree;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;
import sr.speciation.SRangesBirthDeathModel;

import javax.xml.parsers.DocumentBuilderFactory;
import java.io.ByteArrayInputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Store/restore of the SRTree state through the state file, i.e. the resume path:
 * {@code State.toXML} -> {@code SRTree.toXML}, then {@code StateNode.copy}, {@code fromXML}
 * and {@code assignFromFragile} as done by {@code State.fromXML}.
 * <p>
 * The states come from a real MCMC run with all sRanges operators, so sampled ancestors,
 * multi-fossil ranges and both child orientations are all exercised. Three scenarios:
 * (1) restore into a fresh tree, (2) restore repeatedly into the same tree, (3) restore
 * stored states in shuffled order after the run. Restored trees must match the live tree in
 * structure, orientation, heights, ranges and likelihood.
 */
public class SRTreeStateRoundTripTest {

    private static final int CHAIN_LENGTH = 30000;
    private static final int LOG_EVERY = 100;

    // shared, read-only setup
    private TaxonSet taxonSet;
    private TraitSet trait;
    private RealParameter origin, birthRate, deathRate, samplingRate, removalProbability, rho;

    @Test
    public void statesSurviveTheStateFileRoundTrip() throws Exception {
        Randomizer.setSeed(127);

        taxonSet = new TaxonSet();
        taxonSet.setID("taxonSet");
        List<Taxon> taxa = new ArrayList<>();
        for (String name : new String[]{"A_first", "A_last", "B_first", "B_last", "C", "D"}) {
            Taxon t = new Taxon();
            t.setID(name);
            taxa.add(t);
        }
        taxonSet.initByName("taxon", taxa);

        trait = new TraitSet();
        trait.initByName("traitname", "date-backward", "taxa", taxonSet,
                "value", "A_first=3.0,A_last=1.5,B_first=2.5,B_last=0.5,C=2.0,D=0.0");
        trait.setID("dateTrait");

        SRTree tree = newSRTree("Tree.t:tree");

        origin = new RealParameter("6.0");
        origin.initByName("lower", "3.0", "upper", "Infinity");
        origin.setID("origin");

        State state = new State();
        state.initByName("stateNode", tree, "stateNode", origin);
        state.setID("state");

        ConstantPopulation populationModel = new ConstantPopulation();
        populationModel.initByName("popSize", "1.0");
        RandomSRangeTree init = new RandomSRangeTree();
        init.initByName("estimate", false, "initial", tree, "nodetype", SRNode.class.getName(),
                "taxonset", taxonSet, "populationModel", populationModel, "stratigraphicRange", newRanges());

        birthRate = fixed("2.0");
        deathRate = fixed("1.0");
        samplingRate = fixed("0.5");
        removalProbability = fixed("0.9");
        rho = fixed("0.5");
        SRangesBirthDeathModel model = newModel(tree);

        Prior originPrior = new Prior();
        Uniform uniform = new Uniform();
        uniform.initByName("lower", "0.", "upper", "20.0");
        originPrior.initByName("x", origin, "distr", uniform);

        CompoundDistribution posterior = new CompoundDistribution();
        List<Distribution> dists = new ArrayList<>();
        dists.add(model);
        dists.add(originPrior);
        posterior.initByName("distribution", dists);

        SRWilsonBalding wilsonBalding = new SRWilsonBalding();
        wilsonBalding.initByName("tree", tree, "weight", "20.0");
        LeftRightChildSwap childSwap = new LeftRightChildSwap();
        childSwap.initByName("tree", tree, "weight", "20.0");
        SRLeafToSampledAncestorJump saJump = new SRLeafToSampledAncestorJump();
        saJump.initByName("tree", tree, "weight", "20.0", "removalProbability", removalProbability);
        ScaleOperator originScaler = new ScaleOperator();
        originScaler.initByName("parameter", origin, "scaleFactor", "0.9", "weight", "3.0");
        SAScaleOperator rootScaler = new SAScaleOperator();
        rootScaler.initByName("rootOnly", "true", "tree", tree, "scaleFactor", "0.9", "weight", "1.0");

        RoundTripLogger roundTrip = new RoundTripLogger();
        roundTrip.state = state;
        roundTrip.liveTree = tree;
        roundTrip.liveModel = model;
        roundTrip.persistentTarget = newSRTree("Tree.t:tree");
        roundTrip.persistentModel = newModel(roundTrip.persistentTarget);
        roundTrip.initByName("logEvery", Integer.toString(LOG_EVERY), "log", tree);

        // for scenario 3 (models must be created while origin is still above the initial root height)
        SRTree reused = newSRTree("Tree.t:tree");
        SRangesBirthDeathModel reusedModel = newModel(reused);

        MCMC mcmc = new MCMC();
        mcmc.initByName("chainLength", Integer.toString(CHAIN_LENGTH), "state", state,
                "distribution", posterior,
                "operator", wilsonBalding, "operator", childSwap, "operator", saJump,
                "operator", originScaler, "operator", rootScaler,
                "logger", roundTrip);
        mcmc.run();

        assertTrue("no failures expected, first was:\n" + roundTrip.firstFailure, roundTrip.firstFailure == null);
        assertTrue("expected many round trips, got " + roundTrip.roundTrips, roundTrip.roundTrips >= CHAIN_LENGTH / LOG_EVERY);
        assertTrue("MCMC should have visited states with sampled ancestors", roundTrip.statesWithSampledAncestors > 0);
        assertTrue("MCMC should have visited states where a range spans internal nodes", roundTrip.statesWithRangeInternalNodes > 0);
        assertTrue("MCMC should have visited both orientations at the root", roundTrip.rootOrientations.size() == 2);

        // scenario 3: restore every stored state, in shuffled order, into one tree and into fresh trees
        List<Integer> order = new ArrayList<>();
        for (int i = 0; i < roundTrip.storedXML.size(); i++) order.add(i);
        Collections.shuffle(order, new java.util.Random(7));
        for (int i : order) {
            // restore the tree and the origin parameter, as resuming from a state file does
            restoreFromStateXML(roundTrip.storedXML.get(i), reused, origin);
            assertEquals("shuffled restore into reused tree, stored state " + i,
                    roundTrip.storedDescriptions.get(i), describe(reused));
            assertEquals("likelihood after shuffled restore, stored state " + i,
                    roundTrip.storedLogP.get(i), reusedModel.calculateLogP(), 1e-9);

            SRTree fresh = newSRTree("Tree.t:tree");
            restoreFromStateXML(roundTrip.storedXML.get(i), fresh);
            assertEquals("shuffled restore into fresh tree, stored state " + i,
                    roundTrip.storedDescriptions.get(i), describe(fresh));
        }
    }

    // ---------------------------------------------------------------- helpers

    private RealParameter fixed(String value) {
        RealParameter p = new RealParameter(value);
        p.initByName("estimate", false, "lower", "0.0");
        return p;
    }

    /** Every tree needs its own StratigraphicRange objects: SRTree.initSRanges mutates them. */
    private List<StratigraphicRange> newRanges() {
        List<StratigraphicRange> ranges = new ArrayList<>();
        for (String name : new String[]{"A", "B"}) {
            StratigraphicRange r = new StratigraphicRange();
            r.initByName("firstOccurrence", taxonSet.taxonsetInput.get().get(taxonSet.getTaxonIndex(name + "_first")),
                    "lastOccurrence", taxonSet.taxonsetInput.get().get(taxonSet.getTaxonIndex(name + "_last")));
            ranges.add(r);
        }
        return ranges;
    }

    private SRTree newSRTree(String id) {
        SRTree t = new SRTree();
        t.initByName("trait", trait, "taxonset", taxonSet, "nodetype", SRNode.class.getName(),
                "stratigraphicRange", newRanges());
        t.setID(id);
        return t;
    }

    private SRangesBirthDeathModel newModel(SRTree tree) {
        SRangesBirthDeathModel m = new SRangesBirthDeathModel();
        m.initByName("origin", origin, "tree", tree, "birthRate", birthRate, "deathRate", deathRate,
                "samplingRate", samplingRate, "removalProbability", removalProbability, "rho", rho,
                "conditionOnSampling", true);
        return m;
    }

    /** Exactly what State.fromXML does, for the given state nodes only. */
    static void restoreFromStateXML(String stateXML, StateNode... targets) throws Exception {
        Document doc = DocumentBuilderFactory.newInstance().newDocumentBuilder()
                .parse(new ByteArrayInputStream(stateXML.getBytes()));
        doc.normalize();
        NodeList children = doc.getElementsByTagName("*").item(0).getChildNodes();
        for (StateNode target : targets) {
            org.w3c.dom.Node element = null;
            for (int i = 0; i < children.getLength(); i++) {
                org.w3c.dom.Node child = children.item(i);
                if (child.getNodeType() == org.w3c.dom.Node.ELEMENT_NODE
                        && child.getAttributes().getNamedItem("id").getNodeValue().equals(target.getID())) {
                    element = child;
                }
            }
            assertTrue("state XML has no element for " + target.getID(), element != null);
            StateNode scratch = target.copy();
            scratch.fromXML(element);
            target.assignFromFragile(scratch);
        }
    }

    /** Canonical description: oriented structure, sampled ancestors, heights, ranges and node-to-range map. */
    static String describe(SRTree tree) {
        StringBuilder sb = new StringBuilder();
        describeNode(tree.getRoot(), sb);
        sb.append("\nranges:");
        List<String> ranges = new ArrayList<>();
        for (StratigraphicRange r : tree.getSRanges()) {
            StringBuilder rs = new StringBuilder();
            rs.append(r.getFirstOccurrenceID()).append("->").append(r.getLastOccurrenceID())
                    .append(r.isSingleFossilRange() ? " single " : " multi ").append("nodes=").append(r.getNodeNrs());
            ranges.add(rs.toString());
        }
        Collections.sort(ranges);
        sb.append(ranges);
        sb.append("\nrangeOfNode:");
        for (int i = 0; i < tree.getNodeCount(); i++) {
            StratigraphicRange r = tree.getRangeOfNode(tree.getNode(i));
            sb.append(i).append('=').append(r == null ? "-" : r.getFirstOccurrenceID()).append(' ');
        }
        sb.append("\ninternalRangeNodes:");
        List<Integer> internal = new ArrayList<>(tree.getSRangesInternalNodeNrs());
        Collections.sort(internal);
        sb.append(internal);
        return sb.toString();
    }

    private static void describeNode(Node node, StringBuilder sb) {
        if (node.isLeaf()) {
            sb.append(node.getID()).append('#').append(node.getNr());
            if (node.isDirectAncestor()) sb.append('*');
        } else {
            sb.append('(');
            describeNode(node.getLeft(), sb);
            sb.append(',');
            describeNode(node.getRight(), sb);
            sb.append(')').append('#').append(node.getNr());
        }
        sb.append('@').append(String.format(Locale.ROOT, "%.9f", node.getHeight()));
    }

    /** Logger that performs the round trip at every logged state. */
    public static class RoundTripLogger extends Logger {
        State state;
        SRTree liveTree;
        SRangesBirthDeathModel liveModel;
        SRTree persistentTarget;
        SRangesBirthDeathModel persistentModel;

        int roundTrips = 0;
        int statesWithSampledAncestors = 0;
        int statesWithRangeInternalNodes = 0;
        java.util.Set<String> rootOrientations = new java.util.HashSet<>();
        String firstFailure = null;
        List<String> storedXML = new ArrayList<>();
        List<String> storedDescriptions = new ArrayList<>();
        List<Double> storedLogP = new ArrayList<>();

        @Override
        public void init() {
        }

        @Override
        public void close() {
        }

        @Override
        public void log(long sample) {
            if (sample % LOG_EVERY != 0) return;
            try {
                SRTree live = (SRTree) liveTree.getCurrent();
                String expected = describe(live);
                double expectedLogP = liveModel.calculateLogP();
                String xml = state.toXML(sample);

                storedXML.add(xml);
                storedDescriptions.add(expected);
                storedLogP.add(expectedLogP);
                if (live.getDirectAncestorNodeCount() > 0) statesWithSampledAncestors++;
                if (!live.getSRangesInternalNodeNrs().isEmpty()) statesWithRangeInternalNodes++;
                rootOrientations.add(live.getRoot().getLeft().isLeaf() ? "L" : "I");

                // scenario 1: fresh tree
                SRTree fresh = new SRTree();
                fresh.initByName("trait", persistentTarget.m_traitList.get(), "taxonset", persistentTarget.getTaxonset(),
                        "nodetype", SRNode.class.getName(), "stratigraphicRange", copyRangeInputs(persistentTarget));
                fresh.setID(live.getID());
                restoreFromStateXML(xml, fresh);
                check("fresh restore at sample " + sample, expected, describe(fresh));

                // scenario 2: the same target, restored again and again
                restoreFromStateXML(xml, persistentTarget);
                check("repeated restore at sample " + sample, expected, describe(persistentTarget));
                double restoredLogP = persistentModel.calculateLogP();
                if (Math.abs(restoredLogP - expectedLogP) > 1e-9) {
                    fail("likelihood after restore at sample " + sample + ": " + restoredLogP + " vs " + expectedLogP);
                }
                roundTrips++;
            } catch (Exception e) {
                fail("exception at sample " + sample + ": " + e);
            }
        }

        private static List<StratigraphicRange> copyRangeInputs(SRTree template) {
            List<StratigraphicRange> ranges = new ArrayList<>();
            for (StratigraphicRange r : template.stratigraphicRangeInput.get()) {
                if (r.taxonFirstOccurrenceInput.get() == null) continue; // synthetic single ranges are rebuilt
                StratigraphicRange c = new StratigraphicRange();
                c.initByName("firstOccurrence", r.taxonFirstOccurrenceInput.get(),
                        "lastOccurrence", r.taxonLastOccurrenceInput.get());
                ranges.add(c);
            }
            return ranges;
        }

        private void check(String what, String expected, String actual) {
            if (!expected.equals(actual)) {
                fail(what + "\nexpected:\n" + expected + "\nactual:\n" + actual);
            }
        }

        private void fail(String message) {
            if (firstFailure == null) firstFailure = message;
        }
    }
}