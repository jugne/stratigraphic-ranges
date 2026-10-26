package sr.speciation;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;

public class MixedSRangesBirthDeathModel extends SRangesBirthDeathModel {
	public Input<RealParameter> anagenesisRateInput =
			new Input<RealParameter>("anagenesisRate", "Anagenesis rate, default 0", new RealParameter("0.0"));
	public Input<RealParameter> symProportionInput =
			new Input<RealParameter>("symProportion", "Proportion of symmetric birth, default 0.0", new RealParameter("0.0"));

	private RealParameter anagenesisRate, symProportion;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		if(!((SRTree) combinedTree.getTree()).useMixedSpeciationInput.get()) 
			throw new IllegalArgumentException("Trying to use mixed speciation likelihood with non-mixed tree, please set mixed=true in the tree");

		anagenesisRate = anagenesisRateInput.get();
		if(anagenesisRate.getValue() < 0) throw new IllegalArgumentException("Anagenesis rate must be positive");
		anagenesisRate.setLower(0.0);

		symProportion = symProportionInput.get();
		if(symProportion.getValue() < 0.0 || symProportion.getValue() > 1.0) 
			throw new IllegalArgumentException("Asymmetric proportion must be between 0 and 1");
		symProportion.setBounds(0.0, 1.0);
	}

	@Override
	protected double log_q_tilde(double t) {
		double logq = (1 - symProportion.getValue()) * super.log_q_tilde(t); 
		logq -= t * (anagenesisRate.getValue() + symProportion.getValue() * (lambda + mu + psi));
		return logq;
	}
	
	@Override
	protected double q_tilde(double t) {
		double qt = Math.pow(super.q_tilde(t), 1 - symProportion.getValue()); 
		qt *= Math.exp(- t * (anagenesisRate.getValue() + symProportion.getValue() * (lambda + mu + psi)));
		return qt;
	}
	
	@Override
	protected double birthNodeContribution(Node node) {		
		if(((SRNode) node).isBudding()) return Math.log(lambda) + Math.log(1 - symProportion.getValue());
		return Math.log(lambda) + Math.log(symProportion.getValue());
	}
}
