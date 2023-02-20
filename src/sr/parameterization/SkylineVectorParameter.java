package sr.parameterization;


import beast.base.inference.parameter.RealParameter;
import sr.util.Utils;

import java.io.PrintStream;

public class SkylineVectorParameter extends SkylineParameter {

    Double[] values;
    Double[] storedValues;
    double valuesAtTime;

    boolean inputIsScalar;

    public SkylineVectorParameter() { }

    public SkylineVectorParameter(RealParameter changeTimesParam,
                                  RealParameter skylineValuesParam) {
        super(changeTimesParam, skylineValuesParam);
    }

    public SkylineVectorParameter(RealParameter changeTimesParam,
                                  RealParameter skylineValuesParam,
                                  int nTypes) {
        super(changeTimesParam, skylineValuesParam, nTypes, null);
    }

    public SkylineVectorParameter(RealParameter changeTimesParam,
                                  RealParameter skylineValuesParam,
                                  int nTypes, RealParameter origin) {
        super(changeTimesParam, skylineValuesParam, nTypes, origin);
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (skylineValuesInput.get().getDimension() % nIntervals != 0)
            throw new IllegalArgumentException("Value parameter dimension must " +
                    "be a multiple of the number of intervals.");

        int valsPerInterval = skylineValuesInput.get().getDimension()/nIntervals;
        inputIsScalar = valsPerInterval==1;


            if (!inputIsScalar)
                throw new IllegalArgumentException("SkylineVector has an incorrect " +
                        "number of elements.");


        values = new Double[nIntervals];
        storedValues = new Double[nIntervals];

        valuesAtTime = 0.;
    }

    @Override
    protected void updateValues() {

        for (int interval=0; interval<nIntervals; interval++) {
                if (inputIsScalar)
                    values[interval] = skylineValuesInput.get().getArrayValue(interval);
                else
                    values[interval] = skylineValuesInput.get().getArrayValue(interval );

        }

        if (timesAreAges)
            Utils.reverseArray(values);
    }

    /**
     * Retrieve value of vector at a chosen time (not age).
     *
     * @param time when to evaluate the skyline parameter.
     * @return value of the vector at the chosen time.
     */
    protected double getValuesAtTime(double time) {
        update();

        int intervalIdx = getIntervalIdx(time);

        System.arraycopy(values[intervalIdx], 0, valuesAtTime, 0, 1);

        return valuesAtTime;
    }


    @Override
    protected void store() {
        super.store();

        for (int interval=0; interval<nIntervals; interval++)
            System.arraycopy(values[interval], 0, storedValues[interval], 0, 1);
    }

    @Override
    protected void restore() {
        super.restore();

        Double[] tmp;
        tmp = values;
        values = storedValues;
        storedValues = tmp;
    }

    /*
     * Loggable implementation
     */

    @Override
    public void init(PrintStream out) {

        for (int interval=0; interval<nIntervals; interval++) {

            if (interval < nIntervals-1)
                out.print(getID() + "i" + interval + "_endtime\t");

            if (inputIsScalar) {
                out.print(getID());

                if (nIntervals > 1)
                    out.print("i" + interval);

                out.print("\t");

            } else {
                    out.print(getID());

                    if (nIntervals > 1)
                        out.print("i" + interval + "_");

                    out.print("\t");
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {

        for (int interval=0; interval<nIntervals; interval++) {

            if (interval<nIntervals-1)
                out.print(times[interval] + "\t");
            out.print(values[interval] + "\t");

        }
    }

    @Override
    public void close(PrintStream out) {
    }
}
