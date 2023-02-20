package sr.parameterization;

import beast.base.core.Input;

public class FBDParameterization extends Parameterization {

    public Input<SkylineVectorParameter> diversificationRateInput = new Input<>("diversificationRate",
            "Diversification rate skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> turnoverInput = new Input<>("turnover",
            "Turnover skyline.", Input.Validate.REQUIRED);

    public Input<SkylineVectorParameter> samplingProportionInput = new Input<>("samplingProportion",
            "Sampling proportion skyline.", Input.Validate.REQUIRED);

    public Input<TimedParameter> rhoSamplingInput = new Input<>("rhoSampling",
            "Contemporaneous sampling times and probabilities.");


    private double[] birthRateChangeTimes;

    @Override
    public double[] getBirthRateChangeTimes() {
        birthRateChangeTimes = combineAndSortTimes(birthRateChangeTimes,
                    turnoverInput.get().getChangeTimes());

        return birthRateChangeTimes;
    }




    private double[] deathRateChangeTimes;

    @Override
    public double[] getDeathRateChangeTimes() {
        deathRateChangeTimes = combineAndSortTimes(deathRateChangeTimes,
                turnoverInput.get().getChangeTimes(),
                diversificationRateInput.get().getChangeTimes());

        return deathRateChangeTimes;
    }

    private double[] samplingRateChangeTimes;

    @Override
    public double[] getSamplingRateChangeTimes() {
        samplingRateChangeTimes = combineAndSortTimes(samplingRateChangeTimes,
                samplingProportionInput.get().getChangeTimes(),
                diversificationRateInput.get().getChangeTimes(),
                turnoverInput.get().getChangeTimes());

        return samplingRateChangeTimes;
    }

    @Override
    public double[] getRemovalProbChangeTimes() {
        return EMPTY_TIME_ARRAY;
    }

    @Override
    public double[] getRhoSamplingTimes() {
        if (rhoSamplingInput.get() == null)
            return EMPTY_TIME_ARRAY;

        return rhoSamplingInput.get().getTimes();
    }

    @Override
    protected double getBirthRateValues(double time) {
        double res = diversificationRateInput.get().getValuesAtTime(time);
        double toVals = turnoverInput.get().getValuesAtTime(time);

        res /= 1.0 - toVals;

        return res;
    }

    @Override
    protected double getDeathRateValues(double time) {
        double res = diversificationRateInput.get().getValuesAtTime(time);
        double toVals = turnoverInput.get().getValuesAtTime(time);

        res *= toVals/(1.0-toVals);

        return res;
    }

    @Override
    protected double getSamplingRateValues(double time) {
        double res = diversificationRateInput.get().getValuesAtTime(time);
        double sVals = samplingProportionInput.get().getValuesAtTime(time);
        double toVals  = turnoverInput.get().getValuesAtTime(time);

        res *= sVals/(1.0-sVals)
                    * toVals/(1.0-toVals);

        return res;
    }

    @Override
    protected double getRemovalProbValues(double time) {
        return 0.;
    }

    @Override
    protected double getRhoValues(double time) {
        if (rhoSamplingInput.get() != null)
            return rhoSamplingInput.get().getValuesAtTime(time);
        else
            return 0.;
    }
}
