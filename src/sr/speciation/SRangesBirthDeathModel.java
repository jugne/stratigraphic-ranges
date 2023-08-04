package sr.speciation;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.core.Citation;
import beast.base.core.Description;


import beast.base.evolution.tree.TreeDistribution;
import beast.base.inference.parameter.RealParameter;
import sr.evolution.tree.SRTree;
import sr.evolution.sranges.StratigraphicRange;

import static beast.base.core.Log.warning;
import static sr.util.Tools.getEndSubrangeLength;
import static sr.util.Tools.getStartSubrangeLength;

/**
 * @author Alexandra Gavryushkina
 * @author Ugne Stolz
 */


@Description("A variant of the fossilized birth-death model under budding (asymmetric) speciation with stratigraphic ranges")
@Citation("Stadler T, Gavryushkina A, Warnock RCM, Drummond AJ, Heath TA (2017) \n" +
        "The fossilized birth-death model under different types of speciation")
@Citation("Gavryushkina A, Warnock RCM, Drummond AJ, Heath TA, Stadler T (2017) \n" +
        "Bayesian total-evidence dating under the fossilized birth-death model with stratigraphic ranges.")
public class SRangesBirthDeathModel extends TreeDistribution {

    public Input<TransmissionParameterization> parameterizationInput = new Input<>("parameterization", "The parameterization to use.");
    final public Input<Boolean> integrateOverRangesInput = new Input<Boolean>("integrateOverRanges",
            "If true, there should be no samples in between range bounds and these samples are integrated over." +
                    "Otherwise, no integration happens. You should not provide in between samples, " +
                    "when the data they provide is minimal and sampling rate inconsistent with the rest of the tree.", true);

    //     Next two inputs only relevant for transmission applications when we want to include intermediate samples in
//     the range and also set a period before or after sampling where we integrate over the sampling
    final public Input<Boolean> startPeriodInput = new Input<Boolean>("startPeriod",
            "If true, integrate over sampling between empty start sample and first non-empty sample of the range.",
            false);

    final public Input<Boolean> endPeriodInput = new Input<Boolean>("endPeriod",
            "If true, integrate over sampling between last non-empty sample of the range and empty end sample at the end.",
            false);

    public Input<RealParameter> rhoProbability =
            new Input<RealParameter>("rho", "Probability of an individual to be sampled at present", (RealParameter)null);

    // if the tree likelihood is condition on sampling at least one individual then set to true one of the inputs:
    public Input<Boolean> conditionOnSamplingInput = new Input<Boolean>("conditionOnSampling", "the tree " +
            "likelihood is conditioned on sampling at least one individual if condition on origin or at least one individual on both sides of the root if condition on root", false);
//    public Input<Boolean> conditionOnRhoSamplingInput = new Input<Boolean>("conditionOnRhoSampling", "the tree " +
//            "likelihood is conditioned on sampling at least one individual in present if condition on origin or at lease one extant individual on both sides of the root if condition on root", false);

    public Input<Boolean> conditionOnRootInput = new Input<Boolean>("conditionOnRoot", "the tree " +
            "likelihood is conditioned on the root height otherwise on the time of origin", false);

    protected double r[] = new double[2];
    protected double lambda[] = new double[2];
    protected double mu[] = new double[2];
    protected double psi[] = new double[2];
    protected double A[] = new double[2];
    protected double B[] = new double[2];
    protected double origin;
    protected double rho[] = new double[2];
    protected double intervalEndTimes[] = new double[2];
    protected TransmissionParameterization parameterization;


    private double q_i(double t,double t_i, int i) {
        return 4 * Math.exp(-A[i]*(t-t_i)) / Math.pow((1+B[i]) + (1-B[i])*Math.exp(-A[i]*(t-t_i)), 2.0);
    }

    private double log_q_i(double t, double t_i, int i) {
        return Math.log(q_i(t, t_i, i));
    }

    private double q_i_tilde(double t, double t_i, int i) {
        return Math.sqrt(Math.exp(-(lambda[i] + mu[i] + psi[i])*(t-t_i))*q_i(t, t_i, i));
    }

    private double log_q_i_tilde(double t, double t_i, int i) {
        return 0.5*(-(lambda[i] + mu[i] + psi[i])*(t-t_i) + log_q_i(t, t_i, i));
    }

    private double p_i(double t, double t_i, int i) {
        return (lambda[i] + mu[i] + psi[i] - A[i]*((1+B[i])-(1-B[i])*Math.exp(-A[i]*(t-t_i)))/((1+B[i])+(1-B[i])*Math.exp(-A[i]*(t-t_i))))/(2*lambda[i]);
    }

    private double get_p_i(double lambda, double mu, double psi, double A, double B, double t, double t_i) {
        double v = (1 + B) * Math.exp(A * (t - t_i));
        return (lambda + mu + psi - A*((1 + B) -(1-B) * Math.exp(-A * (t - t_i)))/((1 + B)  +(1-B)* Math.exp(-A * (t - t_i))))/(2*lambda);
    }

    private double log_p_i(double t, double t_i, int i) {
        return Math.log(p_i(t, t_i, i));
    }

    private Node findAncestralRangeLastNode(Node node) {
        Node parent = node.getParent();
        if (node.isDirectAncestor()){
            parent = parent.getParent();
            node = node.getParent();
        }
        if (parent == null) {
            return parent;
        } else {
            if (parent.isFake()) {
                return parent;
            } else if (parent.getChild(0) == node) {
                return findAncestralRangeLastNode(parent);
            } else {
                return null;
            }
        }
    }

    protected boolean updateParameters() {
        if (parameterization.samplingChangeTimeInput.get() !=null ){
            intervalEndTimes[0] = parameterization.samplingChangeTimeInput.get().getValue();
        }

            for (int i = dimension-1; i >= 0; i--){
                double p_i_next;
                if (i + 1 < dimension) {
                    p_i_next = get_p_i(parameterization.lambda(i+1),
                            parameterization.mu(i+1),
                            parameterization.psi(i+1),
                            A[i+1], B[i+1], intervalEndTimes[i],intervalEndTimes[i+1]);
                } else {
                    p_i_next = 1.0;
                }

                lambda[i] = parameterization.lambda(i);
                mu[i] = parameterization.mu(i);
                psi[i] = parameterization.psi(i);
                r[i] = parameterization.r(i);

                if (lambdaExceedsMu && lambda[i] <= mu[i]) {
                    return false;
                }

                if (lambda[i] < 0 || mu[i] < 0 || psi[i] < 0) {
                    return false;
                }

                if (i==1 && rhoProbability.get() != null )
                    rho[i] = rhoProbability.get().getValue();
                A[i] = Math.sqrt((lambda[i] - mu[i] - psi[i]) * (lambda[i] - mu[i] - psi[i]) + 4 * lambda[i] * psi[i]);
                B[i] = ((1-2*(1-rho[i])*p_i_next)*lambda[i] +mu[i] + psi[i])/A[i];
            }

            if (!conditionOnRootInput.get()){
                origin = parameterization.origin();
            }  else {
                origin = Double.POSITIVE_INFINITY;
            }




            return true;
        }

    int dimension = 1;
    protected boolean lambdaExceedsMu = false;
    boolean integrateOverRanges;
    boolean useStartPeriod;
    boolean useEndPeriod;
    @Override
    public void initAndValidate() {
        parameterization = parameterizationInput.get();
        if (!integrateOverRangesInput.get() && (startPeriodInput.get() || endPeriodInput.get())) {
            throw new IllegalArgumentException("You can only use startPeriod and endPeriod when " +
                    "integrateOverRanges is true.");
        }
        super.initAndValidate();
        integrateOverRanges = integrateOverRangesInput.get();
        useStartPeriod = startPeriodInput.get();
        useEndPeriod = endPeriodInput.get();
        if (parameterization.samplingChangeTimeInput.get() !=null ){
            intervalEndTimes = new double[2];
            intervalEndTimes[0] = parameterization.samplingChangeTimeInput.get().getValue();
            intervalEndTimes[1] = 0;
            dimension = 2;
        }
    }

    private int getIntervalNumber(double time){
        for (int i=0; i<dimension; i++){
            if (time >= intervalEndTimes[i]){
                return i;
            }
        }
        throw new RuntimeException("Time outside of interval");
    }

    @Override
    public double calculateLogP()
    {
        SRTree tree = (SRTree) treeInput.get();
        int nodeCount = tree.getNodeCount();
        if (!updateParameters())
            return Double.NEGATIVE_INFINITY;


        double x0 = origin;
        double x1=tree.getRoot().getHeight();
        int originInt = getIntervalNumber(x0);

        if (x0 < x1 ) {
            return Double.NEGATIVE_INFINITY;
        }

        if (!conditionOnRootInput.get()){
            logP = log_q_i(x0, intervalEndTimes[0], 0); // -0.10423623551548937
        } else {
            if (tree.getRoot().isFake()){   //when conditioning on the root we assume the process
                //starts at the time of the first branching event and
                //that means that the root can not be a sampled ancestor
                return Double.NEGATIVE_INFINITY;
            } else {
                int i = getIntervalNumber(x1);
                logP = log_q_i(x1, intervalEndTimes[i], i);
            }
        }

        if (conditionOnSamplingInput.get()) {
            int i = getIntervalNumber(x0);
            double p_i = p_i(x0, intervalEndTimes[i], i);
            if (p_i == 1)
                return Double.NEGATIVE_INFINITY; // Following BDSKY's behaviour
            logP -= Math.log(1-p_i); // 0.2981273085767977
        }

        for (int i = 0; i < nodeCount; i++) {
            Node node = tree.getNode(i);
            int j = getIntervalNumber(node.getHeight());
            if (node.isLeaf()) {
                if (dimension>1 && node.getHeight() >= intervalEndTimes[0]){
//                    warning.println("Warning: sampling times before the sampling change time (looking from the root) are not supported yet in this special case implementation.");
                    return Double.NEGATIVE_INFINITY;
                }
                if  (!node.isDirectAncestor())  {
                    Node fossilParent = node.getParent();
                    if (node.getHeight() > intervalEndTimes[j] + 0.000000000005 || rho[j] == 0.) {

                        if ((tree).belongToSameSRange(i, fossilParent.getNr())) {
                            logP += Math.log(psi[j]) - log_q_i_tilde(node.getHeight(), intervalEndTimes[j], j) + log_p_i(node.getHeight(), intervalEndTimes[j], j); // -3.3539504971651484 -7.0136348676169185 -10.67549143975073 -14.33843393100682 -18.023070003492332 -21.71312230023691 -25.411834805439813 -29.110547310642716 -32.82223722234542
                        } else {
                            logP += Math.log(psi[j]) - log_q_i(node.getHeight(), intervalEndTimes[j], j) + log_p_i(node.getHeight(), intervalEndTimes[j], j);
                        }
                    } else {
                        logP += Math.log(rho[j]);
                    }
                }
            } else {
                if (node.isFake()) {
                    if (r[j] == 1) {
                        System.out.println("r = 1 but there are sampled ancestors in the tree");
                        System.exit(0);
                    }
                    logP += Math.log(psi[j]) + Math.log(1 - r[j]); //-64.06813215404485
                    Node parent = node.getParent();

                    Node child = node.getNonDirectAncestorChild();
                    Node DAchild = node.getDirectAncestorChild();
                    if (parent==null && j!=originInt)
                        logP+= log_q_i(intervalEndTimes[originInt], intervalEndTimes[j], j);
                    if (parent != null && tree.belongToSameSRange(parent.getNr(),DAchild.getNr())) {
                        logP += - log_q_i_tilde(node.getHeight(), intervalEndTimes[j], j) + log_q_i(node.getHeight(), intervalEndTimes[j], j); // -84.03106379870141
                    }
                    if (child != null && tree.belongToSameSRange(i,child.getNr())) {
                        logP += - log_q_i(node.getHeight(), intervalEndTimes[j], j) +  log_q_i_tilde(node.getHeight(), intervalEndTimes[j], j); // -66.36368285840686
                    }
                } else {
                    if (node.getParent()==null && j!=originInt)
                        logP+= log_q_i(intervalEndTimes[originInt], intervalEndTimes[j], j);
                    if (node.getParent() != null && getIntervalNumber(node.getParent().getHeight()) != j) {
                        int k = getIntervalNumber(node.getParent().getHeight());
                        logP += log_q_i(intervalEndTimes[k], intervalEndTimes[j], j);
                    }
                    logP += Math.log(lambda[j]) + log_q_i(node.getHeight(), intervalEndTimes[j], j); // -208.80332901142688
                }
            }
        }

        // integrate over fossils in the range. This seems to suggest that we take out the psi in the previous equations
        for (StratigraphicRange range:(tree).getSRanges()) {
            Node first =  tree.getNode(range.getNodeNrs().get(0));
            if (integrateOverRanges && !range.isSingleFossilRange()) {
                double tFirst = first.getHeight();
                int i = getIntervalNumber(tFirst); // this is enough since we assume that ranges only happen in the second period, looking
                // from the origin. This is NOT a proper skyline.
                int rangeSize =  range.getNodeNrs().size();
                if (rangeSize > 1){
                    if (useStartPeriod)
                        logP += (psi[i]*(1-r[i]))*(getStartSubrangeLength(range, tree));
                    if (useEndPeriod)
                        logP += (psi[i]*(1-r[i]))*(getEndSubrangeLength(range, tree));
                } else {
                    double tLast = tree.getNode(range.getNodeNrs().get(range.getNodeNrs().size()-1)).getHeight();
                    logP += (psi[i]*(1-r[i]))*(tFirst - tLast); // -199.3578764324858
                }

            }
            Node ancestralLast = findAncestralRangeLastNode(first);
            if (ancestralLast != null) {
                double tOld = ancestralLast.getHeight();
                double tYoung = first.getHeight();
                int i = getIntervalNumber(tOld); // again, not proper skyline. Only works because we assume that ranges only happen in the second period, looking
                // from the origin.
                logP += Math.log(1-q_i(tYoung, intervalEndTimes[i], i)/q_i_tilde(tYoung, intervalEndTimes[i], i)*q_i_tilde(tOld, intervalEndTimes[i], i)/q_i(tOld, intervalEndTimes[i], i));
            }
        }
        if (logP==Double.POSITIVE_INFINITY)
            System.out.println();
        return logP;
    }
    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

}

