package sr.parameterization;

import beast.base.core.Input;

public class CanonicalParameterization extends Parameterization {

    public Input<SkylineVectorParameter> birthRateInput = new Input<>("birthRate",
            "Birth rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> deathRateInput = new Input<>("deathRate",
            "Death rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> samplingRateInput = new Input<>("samplingRate",
            "Sampling rate skyline.", Input.Validate.REQUIRED);

    public Input<TimedParameter> rhoSamplingInput = new Input<>("rhoSampling",
            "Contemporaneous sampling times and probabilities.");

    public Input<SkylineVectorParameter> removalProbInput = new Input<>("removalProb",
            "Removal prob skyline.", Input.Validate.REQUIRED);


    @Override
    public double[] getBirthRateChangeTimes() {
        return birthRateInput.get().getChangeTimes();
    }

    @Override
    public double[] getDeathRateChangeTimes() {
        return deathRateInput.get().getChangeTimes();
    }

    @Override
    public double[] getSamplingRateChangeTimes() {
        return samplingRateInput.get().getChangeTimes();
    }

    @Override
    public double[] getRemovalProbChangeTimes() {
        return removalProbInput.get().getChangeTimes();
    }

    @Override
    public double[] getRhoSamplingTimes() {
        if (rhoSamplingInput.get() != null)
            return rhoSamplingInput.get().getTimes();
        else
            return EMPTY_TIME_ARRAY;
    }

    @Override
    protected double getBirthRateValues(double time) {
        return birthRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double getDeathRateValues(double time) {
        return deathRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double getSamplingRateValues(double time) {
        return samplingRateInput.get().getValuesAtTime(time);
    }

    @Override
    protected double getRemovalProbValues(double time) {
        return removalProbInput.get().getValuesAtTime(time);
    }

    @Override
    protected double getRhoValues(double time) {
        if (rhoSamplingInput.get() != null)
            return rhoSamplingInput.get().getValuesAtTime(time);
        else
            return 0.;
    }
}
