package sr.speciation;

import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import sa.evolution.speciation.SABDParameterization;

/**
 * Transmission parameterization. For now r input is a bit redundant.
 * Created by Ugne Stolz on 7/09/15.
 * ugne.stolz@protonmail.com
 */
public class TransmissionParameterization extends SABDParameterization {

    public Input<RealParameter> R0Input = new Input<>("R0",
            "Basic reproduction number skyline.", Input.Validate.REQUIRED);

    public Input<RealParameter> becomeUninfectiousRateInput = new Input<>("becomeUninfectiousRate",
            "Become uninfectious rate skyline.", Input.Validate.REQUIRED);

    public Input<RealParameter> samplingProportionInput = new Input<>("samplingProportion",
            "Sampling proportion skyline.", Input.Validate.REQUIRED);

    public Input<RealParameter> removalProbInput = new Input<>("removalProb",
            "Removal prob skyline.", Input.Validate.REQUIRED);

    public Input<RealParameter> originInput =
            new Input<RealParameter>("origin", "The time when the process started", Input.Validate.REQUIRED);

    public double mu() {

        return b()*s()/(1.0-(1-r())*s());
    }

    public double lambda() {

        return r0()*b();
    }

    public double psi() {

        return s()*b()/(1.0-(1-r())*s());
    }

    public double origin() {
        return originInput.get().getValue();
    }
    public double r0() {return R0Input.get().getValue(); }

    public double b() {
        return becomeUninfectiousRateInput.get().getValue();
    }
    public double s() {
        return samplingProportionInput.get().getValue();
    }
    public double r() { return removalProbInput.get().getValue();}

    @Override
    public void initAndValidate() {}
}

