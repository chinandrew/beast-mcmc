package dr.inference.operators.hmc;

import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;

public class StateNoUTurnOperator extends SimpleMCMCOperator implements GibbsOperator, Loggable {

    class Options {
        private double logProbErrorTol = 100.0;
        private int findMax = 100;
        private int maxHeight = 10;
    }

    private final Options options = new Options();

    public StateNoUTurnOperator(ReversibleHMCProvider hmcProvider,
                                boolean adaptiveStepsize,
                                int adaptiveDelay,
                                double weight) {

        this.hmcProvider = hmcProvider;
        this.adaptiveStepsize = adaptiveStepsize;
        this.adaptiveDelay = adaptiveDelay;
        if (hmcProvider instanceof SplitHamiltonianMonteCarloOperator) {
            this.splitHMCmultiplier = ((SplitHamiltonianMonteCarloOperator) hmcProvider).travelTimeMultipler;
            this.splitHMCinner = ((SplitHamiltonianMonteCarloOperator) hmcProvider).inner;
            this.splitHMCouter = ((SplitHamiltonianMonteCarloOperator) hmcProvider).outer;
        }
        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return "HBPS No-U-Turn Operator";
    }

    private WrappedVector drawInertia() {
        double[] inertia_arr = new double[1];
        inertia_arr[0] = MathUtils.nextExponential(1);
        return new WrappedVector.Raw(inertia_arr);
    }

    @Override
    public double doOperation() {

        if (hmcProvider instanceof SplitHamiltonianMonteCarloOperator) {
            updateRS();
            if (splitHMCmultiplier.shouldGetMultiplier(getCount())) {
                ((SplitHamiltonianMonteCarloOperator) hmcProvider).relativeScale = splitHMCmultiplier.getMultiplier();
            }
        }

        final WrappedVector initialPosition = new WrappedVector.Raw(hmcProvider.getInitialPosition());
        final WrappedVector initialMomentum = hmcProvider.drawMomentum();
        final WrappedVector initialInertia = drawInertia();
        final WrappedVector gradient = new WrappedVector.Raw(hmcProvider.getGradientProvider().getGradientLogDensity());

        final ParticleState particleState = new ParticleState(initialPosition, initialMomentum, initialInertia, gradient);
        if (updatePreconditioning) { //todo: should preconditioning, use a schedular
            hmcProvider.providerUpdatePreconditioning();
        }

        if (stepSizeInformation == null) {
            stepSizeInformation = findReasonableStepSize(initialPosition.getBuffer(),
                    hmcProvider.getGradientProvider().getGradientLogDensity(), hmcProvider.getStepSize());
        }
        initializeNumEvents();
        double[] position = takeOneStep(getCount() + 1, particleState);

        hmcProvider.setParameter(position);
        return 0;
    }

    private double[] takeOneStep(long m, ParticleState particleState) {

        double[] endPosition = Arrays.copyOf(particleState.position.getBuffer(), particleState.position.getBuffer().length);

        final double initialJointDensity = hmcProvider.getJointProbability(particleState.momentum);
        double logSliceU = Math.log(getUniform()) + initialJointDensity;

        TreeState trajectoryTree = new TreeState(particleState, 1, true);

        int height = 0;

        while (trajectoryTree.flagContinue) {
            double[] tmp = updateTrajectoryTree(trajectoryTree, height, logSliceU, initialJointDensity);
            if (tmp != null) {
                endPosition = tmp;
            }

            height++;

            if (height > options.maxHeight) {
                trajectoryTree.flagContinue = false;
            }
        }
        if (adaptiveStepsize && getCount() > adaptiveDelay) {
            stepSizeInformation.update(m, trajectoryTree.cumAcceptProb, trajectoryTree.numAcceptProbStates);
            if (printStepsize) System.err.println("step size is " + stepSizeInformation.getStepSize());
        }
        return endPosition;
    }

    private double[] updateTrajectoryTree(TreeState trajectoryTree, int depth, double logSliceU,
                                          double initialJointDensity) {

        double[] endPosition = null;

        final double uniform1 = getUniform();
        int direction = (uniform1 < 0.5) ? -1 : 1;
        TreeState nextTrajectoryTree = buildTree(trajectoryTree.getState(direction),
                direction, logSliceU, depth, stepSizeInformation.getStepSize(), initialJointDensity);

        if (nextTrajectoryTree.flagContinue) {

            final double uniform = getUniform();
            final double acceptProb = (double) nextTrajectoryTree.numNodes / (double) trajectoryTree.numNodes;
            if (uniform < acceptProb) {
                endPosition = nextTrajectoryTree.getSample().position.getBuffer();
            }
        }

        trajectoryTree.mergeNextTree(nextTrajectoryTree, direction);

        return endPosition;
    }

    private TreeState buildTree(ParticleState particleState, int direction,
                                double logSliceU, int height, double stepSize, double initialJointDensity) {

        if (height == 0) {
            return buildBaseCase(particleState, direction, logSliceU, stepSize, initialJointDensity);
        } else {
            return buildRecursiveCase(particleState, direction, logSliceU, height, stepSize,
                    initialJointDensity);
        }
    }


    private TreeState buildBaseCase(ParticleState inParticleState, int direction,
                                    double logSliceU, double stepSize, double initialJointDensity) {
        recordOneBaseCall();
        // Make deep copy of position and momentum
        WrappedVector position = new WrappedVector.Raw(Arrays.copyOf(inParticleState.position.getBuffer(), inParticleState.position.getBuffer().length));
        WrappedVector momentum = new WrappedVector.Raw(Arrays.copyOf(inParticleState.momentum.getBuffer(), inParticleState.momentum.getBuffer().length));
        WrappedVector gradient = new WrappedVector.Raw(Arrays.copyOf(inParticleState.gradient.getBuffer(), inParticleState.gradient.getBuffer().length));
        WrappedVector inertia = new WrappedVector.Raw(Arrays.copyOf(inParticleState.inertia.getBuffer(), inParticleState.inertia.getBuffer().length));
        ParticleState particleState = new ParticleState(position, momentum, inertia, gradient);
        hmcProvider.setParameter(position.getBuffer());

        // "one reversibleHMC integral
        hmcProvider.reversiblePositionMomentumUpdate(position, momentum, inertia, gradient, direction, stepSize);

        recordEvents();

        double logJointProbAfter = hmcProvider.getJointProbability(momentum);

        final int numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);

        final boolean flagContinue = (logSliceU < options.logProbErrorTol + logJointProbAfter);

        // Values for dual-averaging
        final double acceptProb = Math.min(1.0, Math.exp(logJointProbAfter - initialJointDensity));
        final int numAcceptProbStates = 1;

        hmcProvider.setParameter(inParticleState.position.getBuffer());

        return new TreeState(inParticleState, numNodes, flagContinue, acceptProb, numAcceptProbStates);
    }

    private TreeState buildRecursiveCase(ParticleState particleState, int direction,
                                         double logSliceU, int height, double stepSize, double initialJointDensity) {

        TreeState subtree = buildTree(particleState, direction, logSliceU,
                height - 1, // Recursion
                stepSize, initialJointDensity);

        if (subtree.flagContinue) {

            TreeState nextSubtree = buildTree(subtree.getState(direction), direction,
                    logSliceU, height - 1, stepSizeInformation.getStepSize(), initialJointDensity);

            subtree.mergeNextTree(nextSubtree, direction);

        }
        return subtree;
    }

    private static boolean computeStopCriterion(boolean flagContinue, TreeState state) {
        return computeStopCriterion(flagContinue,
                state.getPosition(1).getBuffer(), state.getPosition(-1).getBuffer(),
                state.getMomentum(1).getBuffer(), state.getMomentum(-1).getBuffer());
    }

    private StepSize findReasonableStepSize(double[] initialPosition, double[] initialGradient,
                                            double forcedInitialStepSize) {

        if (forcedInitialStepSize != 0) {
            return new StepSize(forcedInitialStepSize);
        } else {
            double stepSize = 0.1;

            WrappedVector momentum = hmcProvider.drawMomentum();
            WrappedVector inertia = drawInertia();
            int count = 1;
            int dim = initialPosition.length;
            WrappedVector position = new WrappedVector.Raw(Arrays.copyOf(initialPosition, dim));
            WrappedVector gradient = new WrappedVector.Raw(Arrays.copyOf(initialGradient, dim));

            double probBefore = hmcProvider.getJointProbability(momentum);

            hmcProvider.reversiblePositionMomentumUpdate(position, momentum, inertia, gradient, 1, stepSize);

            double probAfter = hmcProvider.getJointProbability(momentum);

            double a = ((probAfter - probBefore) > Math.log(0.5) ? 1 : -1);

            double probRatio = Math.exp(probAfter - probBefore);

            while (Math.pow(probRatio, a) > Math.pow(2, -a)) {

                probBefore = probAfter;
                hmcProvider.reversiblePositionMomentumUpdate(position, momentum, inertia, gradient, 1, stepSize);

                probAfter = hmcProvider.getJointProbability(momentum);
                probRatio = Math.exp(probAfter - probBefore);

                stepSize = Math.pow(2, a) * stepSize;
                count++;

                if (count > options.findMax) {
                    throw new RuntimeException("Cannot find a reasonable step-size in " + options.findMax + " " +
                            "iterations");
                }
            }
            hmcProvider.setParameter(initialPosition);
            return new StepSize(stepSize);
        }
    }


    private static boolean computeStopCriterion(boolean flagContinue,
                                                double[] positionPlus, double[] positionMinus,
                                                double[] momentumPlus, double[] momentumMinus) {

        double[] positionDifference = subtractArray(positionPlus, positionMinus);

        return flagContinue &&
                getDotProduct(positionDifference, momentumMinus) >= 0 &&
                getDotProduct(positionDifference, momentumPlus) >= 0;
    }

    private static double getDotProduct(double[] x, double[] y) {

        assert (x.length == y.length);
        final int dim = x.length;

        double total = 0.0;
        for (int i = 0; i < dim; i++) {
            total += x[i] * y[i];
        }
        return total;
    }

    private static double[] subtractArray(double[] a, double[] b) {

        assert (a.length == b.length);
        final int dim = a.length;

        double[] result = new double[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = a[i] - b[i];
        }

        return result;
    }


    private double getUniform() {
        double tmp;
        if (randomFlg) {
            tmp = MathUtils.nextDouble();
        } else {
            if (count % 10 == 0) {
                ++count;
            }
            tmp = count % 10 / 10.;
            System.err.println(tmp);
            ++count;
        }
        return tmp;
    }

    private class TreeState {

        private TreeState(ParticleState particleState,
                          int numNodes, boolean flagContinue) {
            this(particleState, numNodes, flagContinue, 0.0, 0);
        }

        private TreeState(ParticleState particleState,
                          int numNodes, boolean flagContinue,
                          double cumAcceptProb, int numAcceptProbStates) {
            this.particleState = new ParticleState[3];

            for (int i = 0; i < 3; ++i) {
                this.particleState[i] = particleState;
            }

            // Recursion variables
            this.numNodes = numNodes;
            this.flagContinue = flagContinue;

            // Dual-averaging variables
            this.cumAcceptProb = cumAcceptProb;
            this.numAcceptProbStates = numAcceptProbStates;
        }

        private WrappedVector getPosition(int direction) {
            return particleState[getIndex(direction)].position;
        }

        private WrappedVector getMomentum(int direction) {
            return particleState[getIndex(direction)].momentum;
        }

        private WrappedVector getGradient(int direction) {
            return particleState[getIndex(direction)].gradient;
        }

        private WrappedVector getInertia(int direction) {
            return particleState[getIndex(direction)].inertia;
        }

        private ParticleState getSample() {
            /*
            Returns a state chosen uniformly from the acceptable states along a hamiltonian dynamics trajectory tree.
            The sample is updated recursively while building trees.
            */
            return particleState[getIndex(0)];
        }

        private ParticleState getState(int direction) {
            return particleState[getIndex(direction)];
        }

        private void setState(int direction, WrappedVector position, WrappedVector momentum, WrappedVector inertia, WrappedVector gradient) {
            this.particleState[getIndex(direction)].position = position;
            this.particleState[getIndex(direction)].momentum = momentum;
            this.particleState[getIndex(direction)].inertia = inertia;
            this.particleState[getIndex(direction)].gradient = gradient;
        }


        private void setSample(WrappedVector position, WrappedVector momentum, WrappedVector inertia, WrappedVector gradient) {
            setState(0, position, momentum, inertia, gradient);
        }

        private void setSample(ParticleState particleState) {
            setState(0, particleState.position, particleState.momentum, particleState.inertia, particleState.gradient);
        }

        private int getIndex(int direction) { // valid directions: -1, 0, +1
            assert (direction >= -1 && direction <= 1);
            return direction + 1;
        }

        private void mergeNextTree(TreeState nextTree, int direction) {
            setState(direction,
                    nextTree.getPosition(direction),
                    nextTree.getMomentum(direction),
                    nextTree.getGradient(direction),
                    nextTree.getInertia(direction));


            updateSample(nextTree);

            numNodes += nextTree.numNodes;
            flagContinue = computeStopCriterion(nextTree.flagContinue, this);

            cumAcceptProb += nextTree.cumAcceptProb;
            numAcceptProbStates += nextTree.numAcceptProbStates;
        }

        private void updateSample(TreeState nextTree) {
            double uniform = getUniform();
            if (nextTree.numNodes > 0
                    && uniform < ((double) nextTree.numNodes / (double) (numNodes + nextTree.numNodes))) {
                setSample(nextTree.getSample());
            }
        }

        final private ParticleState[] particleState;
        private int numNodes;
        private boolean flagContinue;

        private double cumAcceptProb;
        private int numAcceptProbStates;
    }

    private void initializeNumEvents() {
        numBaseCalls = 0;
        numBoundaryEvents = 0;
        numGradientEvents = 0;
    }

    private void recordOneBaseCall() {
        numBaseCalls++;
    }

    private void recordEvents() {
        numGradientEvents += hmcProvider.getNumGradientEvent();
        numBoundaryEvents += hmcProvider.getNumBoundaryEvent();
    }

    @Override
    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[4];
        columns[0] = new NumberColumn("base calls") {
            @Override
            public double getDoubleValue() {
                return numBaseCalls;
            }
        };
        columns[1] = new NumberColumn("step size") {
            @Override

            public double getDoubleValue() {
                if (stepSizeInformation != null) return stepSizeInformation.getStepSize();
                else return 0;
            }
        };
        columns[2] = new NumberColumn("gradient events") {
            @Override
            public double getDoubleValue() {
                return numGradientEvents;
            }
        };
        columns[3] = new NumberColumn("boundary events") {
            @Override
            public double getDoubleValue() {
                return numBoundaryEvents;
            }
        };
        return columns;
    }

    private void updateRS() {
        if (splitHMCmultiplier != null && splitHMCmultiplier.shouldUpdateSCM(getCount())) {
            splitHMCmultiplier.updateSCM(splitHMCmultiplier.getInnerCov(), splitHMCinner.getInitialPosition(), getCount());
            splitHMCmultiplier.updateSCM(splitHMCmultiplier.getOuterCov(), splitHMCouter.getInitialPosition(), getCount());
        }
    }

    private ReversibleHMCProvider hmcProvider;
    private StepSize stepSizeInformation;
    private boolean adaptiveStepsize;
    private int adaptiveDelay;
    private int numBaseCalls;
    private int numBoundaryEvents;
    private int numGradientEvents;

    private SplitHMCtravelTimeMultiplier splitHMCmultiplier = null;
    private ReversibleHMCProvider splitHMCinner = null;
    private ReversibleHMCProvider splitHMCouter = null;

    private final boolean updatePreconditioning = false;
    private final boolean printStepsize = false;


    final private boolean randomFlg = true;
    private int count;
}

