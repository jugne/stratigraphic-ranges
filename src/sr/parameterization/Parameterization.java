package sr.parameterization;

import beast.base.inference.CalculationNode;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import sr.util.Utils;

import java.util.*;

/**
 * Full parameterization for a multi-type birth-death skyline model with sampled ancestors.
 *
 * Objects of this class expose a variety of methods for interrogating
 * the canonical (lambda, mu, M, psi, rho, r, t_origin) parameters at different
 * times.  By "time" we mean a number that increases into the future from
 * some point in the past defined by the start of the birth-death process.
 * (We us "ages" to refer to numbers that increase into the past.)
 *
 * In accordance with birth-death skyline convention, the time period spanned by
 * the birth-death process is broken up into intervals at boundaries t_i.
 * Interval i includes the times (t_i-1,t_i], i.e. it does NOT include the time
 * at the earlier boundary.
 */
public abstract class Parameterization extends CalculationNode {

    public Input<Function> originInput = new Input<>("origin",
            "Time between start of process and the end.");

    public Input<Tree> treeInput = new Input<>("tree",
            "If specified, condition on root time rather than origin time.",
            Input.Validate.XOR, originInput);

    private boolean dirty;

    private SortedSet<Double> intervalEndTimesSet = new TreeSet<>(Utils::precisionLimitedComparator);

    private double[] intervalEndTimes, storedIntervalEndTimes;

    private double[] birthRates, deathRates, samplingRates, removalProbs, rhoValues;

    private double[] storedBirthRates, storedDeathRates, storedSamplingRates,
            storedRemovalProbs, storedRhoValues;

    final static double[] EMPTY_TIME_ARRAY = new double[0];


    @Override
    public void initAndValidate() {

        dirty = true;
        update();
    }

    public abstract double[] getBirthRateChangeTimes();
    public abstract double[] getDeathRateChangeTimes();
    public abstract double[] getSamplingRateChangeTimes();
    public abstract double[] getRemovalProbChangeTimes();
    public abstract double[] getRhoSamplingTimes();

    protected abstract double getBirthRateValues(double time);
    protected abstract double getDeathRateValues(double time);
    protected abstract double getSamplingRateValues(double time);
    protected abstract double getRemovalProbValues(double time);
    protected abstract double getRhoValues(double time);


    public double getTotalProcessLength() {
        if (originInput.get() != null)
            return originInput.get().getArrayValue();
        else
            return treeInput.get().getRoot().getHeight();
    }

    public boolean conditionedOnRoot() {
        return originInput.get() == null;
    }

    private void update() {
        if (!dirty)
            return;

        updateModelEventTimes();

        if (birthRates == null) {

            birthRates = new double[intervalEndTimes.length];
            deathRates = new double[intervalEndTimes.length];
            samplingRates = new double[intervalEndTimes.length];
            removalProbs = new double[intervalEndTimes.length];
            rhoValues = new double[intervalEndTimes.length];

            storedBirthRates = new double[intervalEndTimes.length];
            storedDeathRates = new double[intervalEndTimes.length];
            storedSamplingRates = new double[intervalEndTimes.length];
            storedRemovalProbs = new double[intervalEndTimes.length];
            storedRhoValues = new double[intervalEndTimes.length];
        }

        updateValues();

        dirty = false;
    }

    private void addTimes(double[] times) {
        if (times == null)
            return;

        for (double time : times)
            intervalEndTimesSet.add(time);
    }

    private void updateModelEventTimes() {

        intervalEndTimesSet.clear();

        addTimes(getBirthRateChangeTimes());
        addTimes(getDeathRateChangeTimes());
        addTimes(getSamplingRateChangeTimes());
        addTimes(getRemovalProbChangeTimes());
        addTimes(getRhoSamplingTimes());

        intervalEndTimesSet.add(getTotalProcessLength()); // End time of final interval

        if (intervalEndTimes == null) {
            intervalEndTimes = new double[intervalEndTimesSet.size()];
            storedIntervalEndTimes = new double[intervalEndTimesSet.size()];
        }

        List<Double> timeList = new ArrayList<>(intervalEndTimesSet);
        for (int i = 0; i< intervalEndTimesSet.size(); i++)
            intervalEndTimes[i] = timeList.get(i);
    }

    public double[] getIntervalEndTimes() {
        update();

        return intervalEndTimes;
    }

    public int getTotalIntervalCount() {
        update();

        return intervalEndTimes.length;
    }

    /**
     * Finds the index of the time interval t lies in.  Note that if t
     * lies on a boundary between intervals, the interval returned will be
     * the _earlier_ of these two intervals.
     *
     * @param t time for which to identify interval
     * @return index identifying interval.
     */
    public int getIntervalIndex(double t) {
        update();

        int index = Arrays.binarySearch(intervalEndTimes, t);

        if (index < 0)
            index = -index - 1;

        // return at most the index of the last interval (m-1)
        return Math.max(0, Math.min(index, intervalEndTimes.length-1));
    }

    void updateValues() {

        for (int interval = 0; interval < intervalEndTimes.length; interval++) {

            double t = intervalEndTimes[interval];
            System.arraycopy(getBirthRateValues(t), 0, birthRates[interval], 0, 1);
            System.arraycopy(getDeathRateValues(t), 0, deathRates[interval], 0, 1);
            System.arraycopy(getSamplingRateValues(t), 0, samplingRates[interval], 0, 1);
            System.arraycopy(getRemovalProbValues(t), 0, removalProbs[interval], 0, 1);
            System.arraycopy(getRhoValues(t), 0, rhoValues[interval], 0, 1);
        }
    }

    public double[] getBirthRates() {
        update();

        return birthRates;
    }

    public double[] getDeathRates() {
        update();

        return deathRates;
    }

    public double[] getSamplingRates() {
        update();

        return samplingRates;
    }

    public double[] getRemovalProbs() {
        update();

        return removalProbs;
    }

    public double[] getRhoValues() {
        update();

        return rhoValues;
    }


    /**
     * Return time of node, i.e. T - node_age.
     *
     * @param node node whose time to query.
     * @return time of node.
     */
    public double getNodeTime(Node node, double finalSampleOffset) {
        return getTotalProcessLength() - node.getHeight() - finalSampleOffset;
    }

    /**
     * Return the age corresponding to the given time relative to the most recent sample.
     *
     * @param time time to query age for
     * @return age corresponding to time
     */
    public double getAge(double time, double finalSampleOffset) {
        return getTotalProcessLength() - time - finalSampleOffset;
    }

    /**
     * Retrieve index of interval containing node.
     *
     * @param node node whose interval to query.
     * @return index of interval.
     */
    public int getNodeIntervalIndex(Node node, double finalSampleOffset) {
        return getIntervalIndex(getNodeTime(node, finalSampleOffset));
    }


    protected SortedSet<Double> changeTimeSet = new TreeSet<>(Utils::precisionLimitedComparator);

    /**
     * Combine times from individual time arrays, removing duplicates.
     *
     * @param changeTimeArrays One or more arrays to combine.
     * @return combined time array
     */
    protected double[] combineAndSortTimes(double[] destArray, double[] ... changeTimeArrays) {
        changeTimeSet.clear();

        for (double[] changeTimeArray : changeTimeArrays) {
            for (double t : changeTimeArray)
                changeTimeSet.add(t);
        }

        if (destArray == null || destArray.length != changeTimeSet.size())
            destArray = new double[changeTimeSet.size()];

        int i=0;
        for (double changeTime : changeTimeSet)
            destArray[i++] = changeTime;

        return destArray;
    }

    @Override
    protected boolean requiresRecalculation() {
        dirty = true;
        return true;
    }

    @Override
    protected void store() {
        System.arraycopy(intervalEndTimes, 0, storedIntervalEndTimes, 0, intervalEndTimes.length);

        for (int interval=0; interval<intervalEndTimes.length; interval++) {

            System.arraycopy(birthRates[interval], 0, storedBirthRates[interval], 0, 1);
            System.arraycopy(deathRates[interval], 0, storedDeathRates[interval], 0, 1);
            System.arraycopy(samplingRates[interval], 0, storedSamplingRates[interval], 0, 1);
            System.arraycopy(removalProbs[interval], 0, storedRemovalProbs[interval], 0, 1);
            System.arraycopy(rhoValues[interval], 0, storedRhoValues[interval], 0, 1);

        }

        super.store();
    }

    @Override
    protected void restore() {

        double[] scalarTmp;
        double[] vectorTmp;

        scalarTmp = intervalEndTimes;
        intervalEndTimes = storedIntervalEndTimes;
        storedIntervalEndTimes = scalarTmp;

        vectorTmp = birthRates;
        birthRates = storedBirthRates;
        storedBirthRates = vectorTmp;

        vectorTmp = deathRates;
        deathRates = storedDeathRates;
        storedDeathRates = vectorTmp;

        vectorTmp = samplingRates;
        samplingRates = storedSamplingRates;
        storedSamplingRates = vectorTmp;

        vectorTmp = removalProbs;
        removalProbs = storedRemovalProbs;
        storedRemovalProbs = vectorTmp;

        vectorTmp = rhoValues;
        rhoValues = storedRhoValues;
        storedRhoValues = vectorTmp;

        super.restore();
    }
}