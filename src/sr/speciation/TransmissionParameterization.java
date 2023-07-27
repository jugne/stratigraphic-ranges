package sr.speciation;

import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;

/**
 * Transmission parameterization. For now r input is a bit redundant.
 * Created by Ugne Stolz on 7/09/15.
 * ugne.stolz@protonmail.com
 */
public class TransmissionParameterization extends CalculationNode {

    public Input<RealParameter> R0Input = new Input<>("R0",
            "Basic reproduction number skyline.", Input.Validate.REQUIRED);

    public Input<RealParameter> becomeUninfectiousRateInput = new Input<>("becomeUninfectiousRate",
            "Become uninfectious rate skyline.", Input.Validate.REQUIRED);

    public Input<RealParameter> samplingProportionInput = new Input<>("samplingProportion",
            "Sampling proportion skyline.", Input.Validate.REQUIRED);

    public Input<RealParameter> samplingChangeTimeInput = new Input<>("samplingChangeTime",
            "Age (backwards in time) at which sampling started.", Input.Validate.REQUIRED);

    public Input<RealParameter> removalProbInput = new Input<>("removalProb",
            "Removal prob skyline.", Input.Validate.REQUIRED);

    public Input<RealParameter> originInput =
            new Input<RealParameter>("origin", "The time when the process started", Input.Validate.REQUIRED);

    public double mu(int i) {

        return b(i)*(1-s(i))/(1.0-(1-r(i))*s(i));
    }

    public double lambda(int i) {

        return r0(i)*b(i);
    }

    public double psi(int i) {

        return s(i)*b(i)/(1.0-(1-r(i))*s(i));
    }

    public double origin() {
        return originInput.get().getValue();
    }
    public double r0(int i) {
        if (R0Input.get().getDimension()>1)
            return R0Input.get().getArrayValue(i);
        else
            return R0Input.get().getValue(); }

    public double b(int i) {
        if (becomeUninfectiousRateInput.get().getDimension()>1)
            return becomeUninfectiousRateInput.get().getArrayValue(i);
        else
            return becomeUninfectiousRateInput.get().getValue();
    }
    public double s(int i) {
        if (samplingProportionInput.get().getDimension()>1)
            return samplingProportionInput.get().getArrayValue(i);
        else
            return samplingProportionInput.get().getValue();
    }
    public double r(int i) {
        if (removalProbInput.get().getDimension()>1)
            return removalProbInput.get().getArrayValue(i);
        else
            return removalProbInput.get().getValue();
    }

    @Override
    public void initAndValidate() {

    }
}

