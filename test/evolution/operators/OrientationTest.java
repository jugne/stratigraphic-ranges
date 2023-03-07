package evolution.operators;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.inference.*;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.inference.distribution.Prior;
import beast.base.inference.distribution.Uniform;
import beast.base.util.Randomizer;
import org.junit.Assert;
import org.junit.Test;
import sa.evolution.operators.SAScaleOperator;
import sr.evolution.operators.*;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.RandomSRangeTree;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;
import sr.speciation.SRangesBirthDeathModel;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class OrientationTest {
	double tolerance = 0.005;
	double toleranceOriented = 0.025;
	double orientedFrequency = 0.50;

	Integer chainLength = 5000000;
	Integer logEvery = 50;
	int statesLogged = chainLength / logEvery;

	@Test
	public void topologyDistribution_1() throws Exception {
		
//		Randomizer.setSeed(Randomizer.nextInt());
		Randomizer.setSeed(42);

		// state

		TaxonSet taxonSet = new TaxonSet();
		taxonSet.setID("taxonSet");

		Taxon taxon1 = new Taxon();
		taxon1.setID("1");

		Taxon taxon2 = new Taxon();
		taxon2.setID("2");

		Taxon taxon3 = new Taxon();
		taxon3.setID("3");
		taxonSet.initByName("taxon", taxon1, "taxon", taxon2,"taxon", taxon3);


		SRTree tree = new SRTree();

		StratigraphicRange range1 = new StratigraphicRange();
		range1.initByName("firstOccurrence", taxon2, "lastOccurrence", taxon3);
		StratigraphicRange range2 = new StratigraphicRange();
		range2.initByName("firstOccurrence", taxon1, "lastOccurrence", taxon1);

		TraitSet trait = new TraitSet();
		trait.initByName("traitname", "date-backward", "taxa", taxonSet, "value", "1=2,2=1,3=0");
		trait.setID("dateTrait.t:tree");

		tree.initByName("trait", trait, "taxonset", taxonSet,"nodetype", SRNode.class.getName(), "stratigraphicRange", range1, "stratigraphicRange", range2);
		tree.setID("Tree.t:tree");

		RealParameter origin = new RealParameter("4.");
		origin.initByName("lower", "2.0", "upper", "Infinity");
		origin.setID("origin");

		// Set up state:
		State state = new State();
		state.initByName("stateNode", tree, "stateNode", origin);
		state.setID("state");

		ConstantPopulation populationModel = new ConstantPopulation();
		populationModel.setID("ConstantPopulation0.t");
		populationModel.initByName("popSize", "1.0");

		RandomSRangeTree init = new RandomSRangeTree();
		init.initByName("estimate", false, "initial", tree, "nodetype", SRNode.class.getName(),
				"taxonset", taxonSet, "populationModel", populationModel,
				"stratigraphicRange", range1, "stratigraphicRange", range2);
		init.setID("Randomtree");

		// Set up distributions:

		// dist 1
		CompoundDistribution prior = new CompoundDistribution();
		SRangesBirthDeathModel model = new SRangesBirthDeathModel();
		Prior priorDist = new Prior();
		Uniform uniform3 = new Uniform();

		RealParameter birthRate = new RealParameter("2.0");
		RealParameter deathRate = new RealParameter("1.0");
		RealParameter samplingRate = new RealParameter("0.5");
		RealParameter removalProbability = new RealParameter("0.9");
		RealParameter rho = new RealParameter("0.5");

		birthRate.initByName("estimate", false, "lower", "0.0");
		deathRate.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		samplingRate.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		removalProbability.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		rho.initByName("estimate", false, "lower", "0.0", "upper", "1.0");

		model.initByName("origin", origin, "tree", tree, "birthRate", birthRate, "deathRate", deathRate,
				"samplingRate", samplingRate, "removalProbability", removalProbability, "rho", rho, "conditionOnSampling", true);
		uniform3.initByName("upper", "10.0", "lower", "0.");
		priorDist.initByName("x", origin, "distr", uniform3);
		
		List<Distribution> dist21 = new ArrayList<>();
		dist21.add(model);
		dist21.add(priorDist);
		prior.initByName("distribution", dist21);
		prior.setID("prior");

		SiteModel siteModel = new SiteModel();
		RealParameter mutationRate = new RealParameter("1.0");
		RealParameter shape = new RealParameter("1.0");
		RealParameter proportionInvariant = new RealParameter("0.0");
		JukesCantor substModel = new JukesCantor();
		mutationRate.initByName("estimate", false);
		shape.initByName("estimate", false);
		proportionInvariant.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		siteModel.initByName("mutationRate", mutationRate, "shape", shape, "proportionInvariant", proportionInvariant,
				"substModel", substModel);
		
		StrictClockModel branchRateModel = new StrictClockModel();
		RealParameter clockRate = new RealParameter("1.0");
		clockRate.initByName("estimate", false, "lower", "0.0");
		branchRateModel.initByName("clock.rate", clockRate);

		CompoundDistribution posterior = new CompoundDistribution();
		List<Distribution> distributions = new ArrayList<>();
		distributions.add(prior);
		posterior.initByName("distribution", distributions);

		// SRWilsonBalding
		SRWilsonBalding srWilsonBalding = new SRWilsonBalding();
		srWilsonBalding.initByName("tree", tree, "weight", "20.0");

		// LeftRightChildSwap
		LeftRightChildSwap leftRightChildSwap = new LeftRightChildSwap();
		leftRightChildSwap.initByName("tree", tree, "weight", "20.0");

		//LeafToSampledAncestorJump
		SRLeafToSampledAncestorJump LeafToSampledAncestorJump = new SRLeafToSampledAncestorJump();
		LeafToSampledAncestorJump.initByName("tree", tree, "weight", "20.0", "removalProbability", removalProbability);

		// Origin Scaler
		ScaleOperator originScaler = new ScaleOperator();
		originScaler.initByName("parameter", origin, "scaleFactor", "0.9", "weight","3.0");

		// root Scaler
		SAScaleOperator rootScaler = new SAScaleOperator();
		rootScaler.initByName("rootOnly", "true", "tree", tree, "scaleFactor", "0.9", "weight","1.0");



		// Set up logger:
		OrientedThreeSampleTreeLogger treeReport = new OrientedThreeSampleTreeLogger();
		treeReport.initByName("logEvery", logEvery.toString(),
				"burnin", "0",
				"tree", tree,
				"log", tree,
				"silent", true);


		// Set up MCMC:
		MCMC mcmc = new MCMC();
		mcmc.initByName("chainLength", chainLength.toString(),
						"state", state, 
//				"init", init,
				"distribution", posterior,
				"operator", srWilsonBalding,
				"operator", leftRightChildSwap,
				"operator", LeafToSampledAncestorJump,
				"operator", originScaler,
				"operator", rootScaler,
				"logger", treeReport);

		// Run MCMC:
		mcmc.run();

		int[][] frequencies = treeReport.getAnalysis();

		/*
		 * There are two possible non-oriented topologies for three samples,
		 * two stratigraphic ranges (one with the oldest sample and the other with remaining two):
		 *  ((3)2,1) (((3)2)1)
		 */

		// probabilities calculated by separate python script: https://github.com/jugne/sRanges-material/blob/main/integrate_topologies.py
		double[] probs = new double[] {0.6358071491891069, 0.36419285081089314};


		double tolerance = 0.005;
		double toleranceOriented = 0.025;
		double orientedFrequency = 0.50;

//		int s = 0;
//		double[] ss= new double[8];
//		for (int nonOrientedTopologyNr = 0; nonOrientedTopologyNr < 8; nonOrientedTopologyNr++) {
//			s += Arrays.stream(frequencies[nonOrientedTopologyNr]).sum();
//			ss[nonOrientedTopologyNr] = (double) Arrays.stream(frequencies[nonOrientedTopologyNr]).sum() / (double) (statesLogged+1);
//		}

		for (int nonOrientedTopologyNr = 0; nonOrientedTopologyNr < 2; nonOrientedTopologyNr++) {
			int sumTopology = Arrays.stream(frequencies[6+nonOrientedTopologyNr]).sum();
			double probTopology = (double) sumTopology / (double) statesLogged;
			System.out.println("_____________________");
			System.out.println(probs[nonOrientedTopologyNr]);
			System.out.println("_____________________");
			System.out.println(probTopology);
			Assert.assertEquals(probs[nonOrientedTopologyNr], probTopology, tolerance);

			if (6+nonOrientedTopologyNr == 6) {
				for (int j = 0; j < 4; j += 2) {
					double frequency = (double) (frequencies[6+nonOrientedTopologyNr][j]
							+ frequencies[6+nonOrientedTopologyNr][j + 1]) / (double) sumTopology;
					Assert.assertEquals(frequency, orientedFrequency, toleranceOriented);

				}
			}
		}
	}

	@Test
	public void topologyDistribution_2() throws Exception{

		Randomizer.setSeed(Randomizer.nextInt());
//		Randomizer.setSeed(42);

		// state

		TaxonSet taxonSet = new TaxonSet();
		taxonSet.setID("taxonSet");

		Taxon taxon1 = new Taxon();
		taxon1.setID("1");

		Taxon taxon2 = new Taxon();
		taxon2.setID("2");

		Taxon taxon3 = new Taxon();
		taxon3.setID("3");
		taxonSet.initByName("taxon", taxon1, "taxon", taxon2,"taxon", taxon3);


		SRTree tree = new SRTree();

		StratigraphicRange range1 = new StratigraphicRange();
		range1.initByName("firstOccurrence", taxon1, "lastOccurrence", taxon1);
		StratigraphicRange range2 = new StratigraphicRange();
		range2.initByName("firstOccurrence", taxon2, "lastOccurrence", taxon2);
		StratigraphicRange range3 = new StratigraphicRange();
		range3.initByName("firstOccurrence", taxon3, "lastOccurrence", taxon3);


		TraitSet trait = new TraitSet();
		trait.initByName("traitname", "date-backward", "taxa", taxonSet, "value", "1=2,2=1,3=0");
		trait.setID("dateTrait.t:tree");

		tree.initByName("trait", trait, "taxonset", taxonSet,"nodetype", SRNode.class.getName(),
				"stratigraphicRange", range1, "stratigraphicRange", range2, "stratigraphicRange", range3);
		tree.setID("Tree.t:tree");

		RealParameter origin = new RealParameter("4.");
		origin.initByName("lower", "2.0", "upper", "Infinity");
		origin.setID("origin");

		// Set up state:
		State state = new State();
		state.initByName("stateNode", tree, "stateNode", origin);
		state.setID("state");

		ConstantPopulation populationModel = new ConstantPopulation();
		populationModel.setID("ConstantPopulation0.t");
		populationModel.initByName("popSize", "1.0");

		RandomSRangeTree init = new RandomSRangeTree();
		init.initByName("estimate", false, "initial", tree, "nodetype", SRNode.class.getName(),
				"taxonset", taxonSet, "populationModel", populationModel,
				"stratigraphicRange", range1, "stratigraphicRange", range2, "stratigraphicRange", range3);
		init.setID("Randomtree");

		// Set up distributions:

		// dist 1
		CompoundDistribution prior = new CompoundDistribution();
		SRangesBirthDeathModel model = new SRangesBirthDeathModel();
		Prior priorDist = new Prior();
		Uniform uniform3 = new Uniform();


		RealParameter birthRate = new RealParameter("2.0");
		RealParameter deathRate = new RealParameter("1.0");
		RealParameter samplingRate = new RealParameter("0.5");
		RealParameter removalProbability = new RealParameter("0.9");
		RealParameter rho = new RealParameter("0.5");

		birthRate.initByName("estimate", false, "lower", "0.0");
		deathRate.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		samplingRate.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		removalProbability.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		rho.initByName("estimate", false, "lower", "0.0", "upper", "1.0");

		model.initByName("origin", origin, "tree", tree, "birthRate", birthRate, "deathRate", deathRate,
				"samplingRate", samplingRate, "removalProbability", removalProbability, "rho", rho, "conditionOnSampling", true);
		uniform3.initByName("upper", "10.0", "lower", "0.");
		priorDist.initByName("x", origin, "distr", uniform3);

		List<Distribution> dist21 = new ArrayList<>();
		dist21.add(model);
		dist21.add(priorDist);
		prior.initByName("distribution", dist21);
		prior.setID("prior");


		SiteModel siteModel = new SiteModel();
		RealParameter mutationRate = new RealParameter("1.0");
		RealParameter shape = new RealParameter("1.0");
		RealParameter proportionInvariant = new RealParameter("0.0");
		JukesCantor substModel = new JukesCantor();
		mutationRate.initByName("estimate", false);
		shape.initByName("estimate", false);
		proportionInvariant.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		siteModel.initByName("mutationRate", mutationRate, "shape", shape, "proportionInvariant", proportionInvariant,
				"substModel", substModel);

		StrictClockModel branchRateModel = new StrictClockModel();
		RealParameter clockRate = new RealParameter("1.0");
		clockRate.initByName("estimate", false, "lower", "0.0");
		branchRateModel.initByName("clock.rate", clockRate);

		CompoundDistribution posterior = new CompoundDistribution();
		List<Distribution> distributions = new ArrayList<>();
		distributions.add(prior);
		posterior.initByName("distribution", distributions);

		// SRWilsonBalding
		SRWilsonBalding srWilsonBalding = new SRWilsonBalding();
		srWilsonBalding.initByName("tree", tree, "weight", "20.0");

		// LeftRightChildSwap
		LeftRightChildSwap leftRightChildSwap = new LeftRightChildSwap();
		leftRightChildSwap.initByName("tree", tree, "weight", "20.0");

		//LeafToSampledAncestorJump
		SRLeafToSampledAncestorJump LeafToSampledAncestorJump = new SRLeafToSampledAncestorJump();
		LeafToSampledAncestorJump.initByName("tree", tree, "weight", "20.0", "removalProbability", removalProbability);

		// Origin Scaler
		ScaleOperator originScaler = new ScaleOperator();
		originScaler.initByName("parameter", origin, "scaleFactor", "0.9", "weight","3.0");

		// root Scaler
		SAScaleOperator rootScaler = new SAScaleOperator();
		rootScaler.initByName("rootOnly", "true", "tree", tree, "scaleFactor", "0.9", "weight","1.0");

		Integer chainLength = 5000000;
		Integer logEvery = 10;
		int statesLogged = chainLength / logEvery;

		// Set up logger:
		OrientedThreeSampleTreeLogger treeReport = new OrientedThreeSampleTreeLogger();
		treeReport.initByName("logEvery", logEvery.toString(),
				"burnin", "0",
				"tree", tree,
				"log", tree,
				"silent", true);


		// Set up MCMC:
		MCMC mcmc = new MCMC();
		mcmc.initByName("chainLength", chainLength.toString(),
				"state", state,
				"distribution", posterior,
				"operator", srWilsonBalding,
				"operator", leftRightChildSwap,
				"operator", LeafToSampledAncestorJump,
				"operator", originScaler,
				"operator", rootScaler,
				"logger", treeReport);

		// Run MCMC:
		mcmc.run();

		int[][] frequencies = treeReport.getAnalysis();

		/*
		 * There are eight possible non-oriented topologies for three samples,
		 * each being a stratigraphic range:
		 *  1-((3,2),1); 2-((3,1),2); 3-(3,(2,1)); 4-((3,2)1); 5-(3,(2)1); 6-(2,(3)1); 7-(1,(3)2); 8-(((3)2)1)
		 */

		// probabilities calculated by separate python script: https://github.com/jugne/sRanges-material/blob/main/integrate_topologies.py
		double[] probs = new double[] {0.4139885939652612, 0.019549249188852594, 0.019549249188852594, 0.2851624323374122,
				0.016791078271745332, 0.02559417041533165, 0.13956962629990535, 0.07979560033263923};




//		int s = 0;
//		double[] ss= new double[8];
//		for (int nonOrientedTopologyNr = 0; nonOrientedTopologyNr < 8; nonOrientedTopologyNr++) {
//			s += Arrays.stream(frequencies[nonOrientedTopologyNr]).sum();
//			ss[nonOrientedTopologyNr] = (double) Arrays.stream(frequencies[nonOrientedTopologyNr]).sum() / (double) (statesLogged+1);
//		}

		for (int nonOrientedTopologyNr = 0; nonOrientedTopologyNr < 8; nonOrientedTopologyNr++) {
			int sumTopology = Arrays.stream(frequencies[nonOrientedTopologyNr]).sum();
			double probTopology = (double) sumTopology / (double) statesLogged;
			System.out.println("_____________________");
			System.out.println(probs[nonOrientedTopologyNr]);
			System.out.println("_____________________");
			System.out.println(probTopology);
			Assert.assertEquals(probs[nonOrientedTopologyNr], probTopology, tolerance);

//			if (6+nonOrientedTopologyNr == 6) {
//				for (int j = 0; j < 4; j += 2) {
//					double frequency = (double) (frequencies[6+nonOrientedTopologyNr][j]
//							+ frequencies[6+nonOrientedTopologyNr][j + 1]) / (double) sumTopology;
//					Assert.assertEquals(frequency, orientedFrequency, toleranceOriented);
//
//				}
			}
	}

	@Test
	public void topologyDistribution_3() throws Exception{

		Randomizer.setSeed(Randomizer.nextInt());
//		Randomizer.setSeed(42);

		// state

		TaxonSet taxonSet = new TaxonSet();
		taxonSet.setID("taxonSet");

		Taxon taxon1 = new Taxon();
		taxon1.setID("1");

		Taxon taxon2 = new Taxon();
		taxon2.setID("2");

		Taxon taxon3 = new Taxon();
		taxon3.setID("3");
		taxonSet.initByName("taxon", taxon1, "taxon", taxon2,"taxon", taxon3);


		SRTree tree = new SRTree();

		StratigraphicRange range1 = new StratigraphicRange();
		range1.initByName("firstOccurrence", taxon1, "lastOccurrence", taxon2);
		StratigraphicRange range2 = new StratigraphicRange();
		range2.initByName("firstOccurrence", taxon3, "lastOccurrence", taxon3);

		TraitSet trait = new TraitSet();
		trait.initByName("traitname", "date-backward", "taxa", taxonSet, "value", "1=2,2=1,3=0");
		trait.setID("dateTrait.t:tree");

		tree.initByName("trait", trait, "taxonset", taxonSet,"nodetype", SRNode.class.getName(), "stratigraphicRange", range1, "stratigraphicRange", range2);
		tree.setID("Tree.t:tree");

		RealParameter origin = new RealParameter("4.");
		origin.initByName("lower", "2.0", "upper", "Infinity");
		origin.setID("origin");

		// Set up state:
		State state = new State();
		state.initByName("stateNode", tree, "stateNode", origin);
		state.setID("state");

		ConstantPopulation populationModel = new ConstantPopulation();
		populationModel.setID("ConstantPopulation0.t");
		populationModel.initByName("popSize", "1.0");

		RandomSRangeTree init = new RandomSRangeTree();
		init.initByName("estimate", false, "initial", tree, "nodetype", SRNode.class.getName(),
				"taxonset", taxonSet, "populationModel", populationModel,
				"stratigraphicRange", range1, "stratigraphicRange", range2);
		init.setID("Randomtree");

		// Set up distributions:

		// dist 1
		CompoundDistribution prior = new CompoundDistribution();
		SRangesBirthDeathModel model = new SRangesBirthDeathModel();
		Prior priorDist = new Prior();
		Uniform uniform3 = new Uniform();


		RealParameter birthRate = new RealParameter("2.0");
		RealParameter deathRate = new RealParameter("1.0");
		RealParameter samplingRate = new RealParameter("0.5");
		RealParameter removalProbability = new RealParameter("0.9");
		RealParameter rho = new RealParameter("0.5");

		birthRate.initByName("estimate", false, "lower", "0.0");
		deathRate.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		samplingRate.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		removalProbability.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		rho.initByName("estimate", false, "lower", "0.0", "upper", "1.0");

		model.initByName("origin", origin, "tree", tree, "birthRate", birthRate, "deathRate", deathRate,
				"samplingRate", samplingRate, "removalProbability", removalProbability, "rho", rho, "conditionOnSampling", true);
		uniform3.initByName("upper", "1000.0", "lower", "0.");
		priorDist.initByName("x", origin, "distr", uniform3);

		List<Distribution> dist21 = new ArrayList<>();
		dist21.add(model);
		dist21.add(priorDist);
		prior.initByName("distribution", dist21);
		prior.setID("prior");


		SiteModel siteModel = new SiteModel();
		RealParameter mutationRate = new RealParameter("1.0");
		RealParameter shape = new RealParameter("1.0");
		RealParameter proportionInvariant = new RealParameter("0.0");
		JukesCantor substModel = new JukesCantor();
		mutationRate.initByName("estimate", false);
		shape.initByName("estimate", false);
		proportionInvariant.initByName("estimate", false, "lower", "0.0", "upper", "1.0");
		siteModel.initByName("mutationRate", mutationRate, "shape", shape, "proportionInvariant", proportionInvariant,
				"substModel", substModel);

		StrictClockModel branchRateModel = new StrictClockModel();
		RealParameter clockRate = new RealParameter("1.0");
		clockRate.initByName("estimate", false, "lower", "0.0");
		branchRateModel.initByName("clock.rate", clockRate);

		CompoundDistribution posterior = new CompoundDistribution();
		List<Distribution> distributions = new ArrayList<>();
		distributions.add(prior);
		posterior.initByName("distribution", distributions);

		// SRWilsonBalding
		SRWilsonBalding srWilsonBalding = new SRWilsonBalding();
		srWilsonBalding.initByName("tree", tree, "weight", "20.0");

		// LeftRightChildSwap
		LeftRightChildSwap leftRightChildSwap = new LeftRightChildSwap();
		leftRightChildSwap.initByName("tree", tree, "weight", "20.0");

		//LeafToSampledAncestorJump
		SRLeafToSampledAncestorJump LeafToSampledAncestorJump = new SRLeafToSampledAncestorJump();
		LeafToSampledAncestorJump.initByName("tree", tree, "weight", "20.0", "removalProbability", removalProbability);

		// Origin Scaler
		ScaleOperator originScaler = new ScaleOperator();
		originScaler.initByName("parameter", origin, "scaleFactor", "0.9", "weight","3.0");

		// root Scaler
		SAScaleOperator rootScaler = new SAScaleOperator();
		rootScaler.initByName("rootOnly", "true", "tree", tree, "scaleFactor", "0.9", "weight","1.0");

		Integer chainLength = 5000000;
		Integer logEvery = 50;
		int statesLogged = chainLength / logEvery;

		// Set up logger:
		OrientedThreeSampleTreeLogger treeReport = new OrientedThreeSampleTreeLogger();
		treeReport.initByName("logEvery", logEvery.toString(),
				"burnin", "0",
				"tree", tree,
				"log", tree,
				"silent", true);


		// Set up MCMC:
		MCMC mcmc = new MCMC();
		mcmc.initByName("chainLength", chainLength.toString(),
				"state", state,
				"distribution", posterior,
				"operator", srWilsonBalding,
				"operator", leftRightChildSwap,
				"operator", LeafToSampledAncestorJump,
				"operator", originScaler,
				"operator", rootScaler,
				"logger", treeReport);

		// Run MCMC:
		mcmc.run();

		int[][] frequencies = treeReport.getAnalysis();

		/*
		 * There are three possible non-oriented topologies for three samples and one range between top two samples:
		 *  4-((3,2)1); 5-(3,(2)1); 8-(((3)2)1)
		 */

		// probabilities calculated by separate python script: https://github.com/jugne/sRanges-material/blob/main/integrate_topologies.py
		double[] probs = new double[] {0, 0, 0, 0.5390864435975525, 0.08012736035537182, 0, 0, 0.38078619604707575};

		for (int nonOrientedTopologyNr = 0; nonOrientedTopologyNr < 8; nonOrientedTopologyNr++) {
			int sumTopology = Arrays.stream(frequencies[nonOrientedTopologyNr]).sum();
			double probTopology = (double) sumTopology / (double) statesLogged;
			System.out.println("_____________________");
			System.out.println(probs[nonOrientedTopologyNr]);
			System.out.println("_____________________");
			System.out.println(probTopology);
			Assert.assertEquals(probs[nonOrientedTopologyNr], probTopology, tolerance);

//			if (6+nonOrientedTopologyNr == 6) {
//				for (int j = 0; j < 4; j += 2) {
//					double frequency = (double) (frequencies[6+nonOrientedTopologyNr][j]
//							+ frequencies[6+nonOrientedTopologyNr][j + 1]) / (double) sumTopology;
//					Assert.assertEquals(frequency, orientedFrequency, toleranceOriented);
//
//				}
//			}
		}
	}
	public static class OrientedThreeSampleTreeLogger extends Logger {
		public Input<Integer> burninInput = new Input<>("burnin",
				"Number of samples to skip (burn in)", Input.Validate.REQUIRED);

		public Input<Boolean> silentInput = new Input<>("silent",
				"Don't display final report.", false);

		final public Input<SRTree> treeInput = new Input<>("tree",
				"The species tree to be logged.", Validate.REQUIRED);


		private DecimalFormat df;


		SRTree SRtree;
		int m_nEvery = 1;
		int burnin;
		boolean silent = false;
		int[][] duplicate = new int[8][4];

		// There are 2 possible non-oriented topologies in our case. Without any restrictions there are 8.
		// Each of those can be produced in Beast2 in 4 different orientation if no restrictions are applied.
		// With stratigraphic ranges model restrictions the first topology can have 2 orientations
		// and second can have only one orientation. We still track all four orientations in case something with
		// restrictions goes wrong. This would then be caught up at the testing step.
		int[][] freq = new int[8][4];

		String number = ":(\\d+(\\.\\d+)?)(E-\\d+)?";

		// Beast2 will log all topologies in one of these three oriented groups.
		// The reason for this, is that a sampled ancestor node is added as a fake child
		// with edge of length 0;
		HashMap<Integer, Pattern> rx_1 = new HashMap<>() {
			{ // Group 1
				// (1,(2,3))
				put(1, Pattern.compile("\\(1" + number + ",\\(2" + number + ",3" + number + "\\)(.*)\\)(.*)"));
				// (1,(3,2))
				put(2, Pattern.compile("\\(1" + number + ",\\(3" + number + ",2" + number + "\\)(.*)\\)(.*)"));
				// ((2,3),1)
				put(3, Pattern.compile("\\(\\(2" + number + ",3" + number + "\\)(.*),1" + number + "\\)(.*)"));
				// ((3,2),1)
				put(4, Pattern.compile("\\(\\(3" + number + ",2" + number + "\\)(.*),1" + number + "\\)(.*)"));
			}
		};

		HashMap<Integer, Pattern> rx_2 = new HashMap<>() {
			{ // Group 2
				// (2, (1,3))
				put(1, Pattern.compile("\\(2" + number + ",\\(1" + number + ",3" + number + "\\)(.*)\\)(.*)"));
				// (2,(3,1))
				put(2, Pattern.compile("\\(2" + number + ",\\(3" + number + ",1" + number + "\\)(.*)\\)(.*)"));
				// ((1,3),2)
				put(3, Pattern.compile("\\(\\(1" + number + ",3" + number + "\\)(.*),2" + number + "\\)(.*)"));
				// ((3,1),2)
				put(4, Pattern.compile("\\(\\(3" + number + ",1" + number + "\\)(.*),2" + number + "\\)(.*)"));
			}
		};

		HashMap<Integer, Pattern> rx_3 = new HashMap<>() {
			{ // Group 3
				// (3,(1,2))
				put(1, Pattern.compile("\\(3" + number + ",\\(1" + number + ",2" + number + "\\)(.*)\\)(.*)"));
				// (3,(2,1))
				put(2, Pattern.compile("\\(3" + number + ",\\(2" + number + ",1" + number + "\\)(.*)\\)(.*)"));
				// ((2,1),3)
				put(3, Pattern.compile("\\(\\(2" + number + ",1" + number + "\\)(.*),3" + number + "\\)(.*)"));
				// ((1,2),3)
				put(4, Pattern.compile("\\(\\(1" + number + ",2" + number + "\\)(.*),3" + number + "\\)(.*)"));
			}
		};

		HashMap<Integer, HashMap<Integer, Pattern>> rx = new HashMap<>(){{
			put(1, rx_1);
			put(2, rx_2);
			put(3, rx_3);
		}};

		public OrientedThreeSampleTreeLogger() {
		}

		@Override
		public void initAndValidate() {

			List<BEASTObject> loggers = loggersInput.get();
			final int nLoggers = loggers.size();
			if (nLoggers == 0) {
				throw new IllegalArgumentException("Logger with nothing to log specified");
			}

			if (everyInput.get() != null)
				m_nEvery = everyInput.get();

			burnin = burninInput.get();

			if (silentInput.get() != null)
				silent = silentInput.get();

			SRtree = treeInput.get();

		}

		@Override
		public void init() {
		}

		@Override
		public void log(long nSample) {

			if ((nSample % m_nEvery > 0) || nSample < burnin)
				return;

			SRTree tree = (SRTree) SRtree.getCurrent();
			String newick = toNewick(tree.getRoot());

			boolean firstAncestor = false;
			boolean secondAncestor = false;
			Node first = null;
			List<Node> leaves = tree.getRoot().getAllLeafNodes();
			for (Node l : leaves) {
				if (l.getID().equals("1")) {
					first = l; 
					firstAncestor = l.isDirectAncestor();
				}
				if (l.getID().equals("2")) {
					secondAncestor = l.isDirectAncestor();
				}
			}
			

			for (int i = 0; i < freq.length; i++) {
				duplicate[i] = Arrays.copyOf(freq[i], freq[i].length);
			}


			if (!secondAncestor) {
				if (!firstAncestor) {
					// Non-oriented tree ((3,2),1) number 0, prob 77.8327
					// Non-oriented tree ((3,1),2) number 1, prob 4.3189
					// Non-oriented tree (3,(2,1)) number 2, prob 4.3189
					for (int nonOrientedTopologyNr=0; nonOrientedTopologyNr<3;nonOrientedTopologyNr++) {
						for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
							Matcher m = rx.get(nonOrientedTopologyNr + 1).get(i).matcher(newick);
							if (m.matches()) {
								freq[nonOrientedTopologyNr][i - 1] += 1;
							}
						}
					}
				} else {
					// Non-oriented tree ((3,2)1) number 3:
					// node 1 is an ancestral node of both leaf nodes
					if (first.getParent().getNonDirectAncestorChild().getAllChildNodesAndSelf().size() == 3) {
						int nonOrientedTopologyNr = 3;
						int orientedNr = 1;
						for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
							Matcher m = rx.get(orientedNr).get(i).matcher(newick);
							if (m.matches()) {
								freq[nonOrientedTopologyNr][i - 1] += 1;
							}
						}
					}// Non-oriented tree (3,(2)1) number 4:
					// node 1 is an ancestral node of leaf node 2
					else if (first.getParent().getNonDirectAncestorChild().getID() != null &&
							first.getParent().getNonDirectAncestorChild().getID().equals("2")) {
						int nonOrientedTopologyNr = 4;
						int orientedNr = 3;
						for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
							Matcher m = rx.get(orientedNr).get(i).matcher(newick);
							if (m.matches()) {
								freq[nonOrientedTopologyNr][i - 1] += 1;
							}
						}
					} // Non-oriented tree (2,(3)1) number 5:
					// node 1 is an ancestral node of leaf node 3
					else if (first.getParent().getNonDirectAncestorChild().getID() != null &&
							first.getParent().getNonDirectAncestorChild().getID().equals("3")) {
						int nonOrientedTopologyNr = 5;
						int orientedNr = 2;
						for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
							Matcher m = rx.get(orientedNr).get(i).matcher(newick);
							if (m.matches()) {
								freq[nonOrientedTopologyNr][i - 1] += 1;
							}
						}
					}
				}
			}// Non-oriented tree (1,(3)2) number 6:
			// node 2 is an ancestral node of leaf node 3 
			else if (!firstAncestor) {
				int nonOrientedTopologyNr = 6;
				int orientedNr = 1;
				for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
					Matcher m = rx.get(orientedNr).get(i).matcher(newick);
					if (m.matches()) {
						freq[nonOrientedTopologyNr][i - 1] += 1;
					}
				}
			} // Non-oriented tree (((3)2)1) number 7:
				// 1 is ancestral of 2 and 3, 2 is ancestral of 3
			else {
				int nonOrientedTopologyNr = 7;
				int orientedNr = 1;
				for (int i = 1; i <= freq[nonOrientedTopologyNr].length; i++) {
					Matcher m = rx.get(orientedNr).get(i).matcher(newick);
					if (m.matches()) {
						freq[nonOrientedTopologyNr][i - 1] += 1;
					}
				}
			}

			if (Arrays.deepEquals(duplicate, freq)) {
				throw new IllegalArgumentException("Not Assigned: " + newick);
			}

			int sumBefore = 0;
			int sumAfter = 0;

			for (int nonOrientedTopologyNr = 0; nonOrientedTopologyNr < 8; nonOrientedTopologyNr++) {
				sumAfter += Arrays.stream(freq[nonOrientedTopologyNr]).sum();
				sumBefore += Arrays.stream(duplicate[nonOrientedTopologyNr]).sum();
			}

			if (sumAfter != sumBefore + 1) {
				throw new IllegalArgumentException("Tree assigned to more than one topology");
			}

		}

		@Override
		public void close() {
		}

		/**
		 * Obtain all frequencies
		 *
		 * @return topology frequencies.
		 */
		public int[][] getAnalysis() {
			return freq;
		}

		// Get Newick string, without sorting the nodes before
		String toNewick(Node node) {
			StringBuffer buf = new StringBuffer();
			if (node.getLeft() != null) {
				buf.append("(");
				buf.append(toNewick(node.getLeft()));
				if (node.getRight() != null) {
					buf.append(',');
					buf.append(toNewick(node.getRight()));
				}
				buf.append(")");
			}
			if (node.getID() == null) {
				buf.append(node.getNr() + 1);
			} else {
				buf.append(node.getID());
			}
			buf.append(":");

			double nodeLength;
			if (node.isRoot()) {
				nodeLength = SRtree.getRoot().getHeight() - node.getHeight();
			} else {
				nodeLength = node.getLength();
			}
			appendDouble(buf, nodeLength);

			return buf.toString();
		}

		/**
		 * Appends a double to the given StringBuffer, formatting it using the private
		 * DecimalFormat instance, if the input 'dp' has been given a non-negative
		 * integer, otherwise just uses default formatting.
		 * 
		 * @param buf
		 * @param d
		 */
		private void appendDouble(StringBuffer buf, double d) {
			if (df == null) {
				buf.append(d);
			} else {
				buf.append(df.format(d));
			}
		}
	}
}

