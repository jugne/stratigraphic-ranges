package sr.parameterization;

import beast.base.core.Input;

public class EpiParameterization extends Parameterization {

    public Input<SkylineVectorParameter> R0Input = new Input<>("R0",
            "Basic reproduction number skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> becomeUninfectiousRateInput = new Input<>("becomeUninfectiousRate",
            "Become uninfectious rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> samplingProportionInput = new Input<>("samplingProportion",
            "Sampling proportion skyline.", Input.Validate.REQUIRED);

    public Input<TimedParameter> rhoSamplingInput = new Input<>("rhoSampling",
            "Contemporaneous sampling times and probabilities.");

    public Input<SkylineVectorParameter> removalProbInput = new Input<>("removalProb",
            "Removal prob skyline.", Input.Validate.REQUIRED);


    private double[] birthRateChangeTimes;

    @Override
    public double[] getBirthRateChangeTimes() {
        birthRateChangeTimes = combineAndSortTimes(birthRateChangeTimes,
                R0Input.get().getChangeTimes(),
                becomeUninfectiousRateInput.get().getChangeTimes());

        return birthRateChangeTimes;
    }

    private double[] deathRateChangeTimes;

    @Override
    public double[] getDeathRateChangeTimes() {
        deathRateChangeTimes = combineAndSortTimes(deathRateChangeTimes,
                becomeUninfectiousRateInput.get().getChangeTimes(),
                samplingProportionInput.get().getChangeTimes(),
                removalProbInput.get().getChangeTimes());

        return deathRateChangeTimes;
    }

    private double[] samplingRateChangeTimes;

    @Override
    public double[] getSamplingRateChangeTimes() {
        samplingRateChangeTimes = combineAndSortTimes(samplingRateChangeTimes,
                samplingProportionInput.get().getChangeTimes(),
                becomeUninfectiousRateInput.get().getChangeTimes());

        return samplingRateChangeTimes;
    }

    @Override
    public double[] getRemovalProbChangeTimes() {
        return removalProbInput.get().getChangeTimes();
    }

    @Override
    public double[] getRhoSamplingTimes() {
        if (rhoSamplingInput.get() == null)
            return EMPTY_TIME_ARRAY;

        return rhoSamplingInput.get().getTimes();
    }


    @Override
    protected double getBirthRateValues(double time) {
        double res = R0Input.get().getValuesAtTime(time);
        double buVals = becomeUninfectiousRateInput.get().getValuesAtTime(time);

        res *= buVals;

        return res;
    }

    @Override
    protected double getDeathRateValues(double time) {
        double res = becomeUninfectiousRateInput.get().getValuesAtTime(time);
        double samplingProp = samplingProportionInput.get().getValuesAtTime(time);
        double removalProb = removalProbInput.get().getValuesAtTime(time);

        res *= (1 - samplingProp)
                    / (1.0 - (1.0-removalProb)*samplingProp);

        return res;
    }

    @Override
    protected double getSamplingRateValues(double time) {
        double res = samplingProportionInput.get().getValuesAtTime(time);
        double buRate  = becomeUninfectiousRateInput.get().getValuesAtTime(time);
        double removalProb  = removalProbInput.get().getValuesAtTime(time);

        res = res*buRate/(1 - (1-removalProb)*res);

        return res;
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
