/*
 * NewHamiltonianMonteCarloOperator.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.TaskPool;
import dr.util.Transform;
import dr.xml.Reportable;

import java.util.function.BinaryOperator;

import static dr.math.matrixAlgebra.ReadableVector.Utils.innerProduct;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */


/**
 * Used for NUTS, which requires a final inertia update.
 * The HamiltonianBouncyParticleOperator is used for standard HBPS
 * without nuts and omits this last step.
 */
public class ReversibleHamiltonianBouncyParticleOperator extends AbstractParticleOperator implements Reportable, ReversibleHMCProvider {

    public ReversibleHamiltonianBouncyParticleOperator(GradientWrtParameterProvider gradientProvider,
                                             PrecisionMatrixVectorProductProvider multiplicationProvider,
                                             PrecisionColumnProvider columnProvider,
                                             double weight, Options runtimeOptions, NativeCodeOptions nativeOptions,
                                             boolean refreshVelocity, Parameter mask,
                                             MassPreconditioner massPreconditioner,
                                             MassPreconditionScheduler.Type preconditionSchedulerType) {
        super(gradientProvider, multiplicationProvider, columnProvider, weight, runtimeOptions, nativeOptions,
                refreshVelocity, mask, null, massPreconditioner, preconditionSchedulerType);
    }


    @SuppressWarnings("all")
    private double getBounceTime(double v_phi_v, double v_phi_x, double[] inertia) {
        double a = v_phi_v;
        double b = 2.0 * v_phi_x;
        double c = - 2.0 * inertia[0];
        return (-b + Math.sqrt(b * b - 4 * a * c)) / 2 / a;
    }

    private double drawInitialInertia() {
        return MathUtils.nextExponential(1);
    }

    private WrappedVector drawInitialVelocity() {

        if (storedVelocity != null) {
            return storedVelocity;
        } else {
            ReadableVector mass = preconditioning.mass;
            double[] velocity = new double[mass.getDim()];

            for (int i = 0, len = velocity.length; i < len; i++) {
                velocity[i] = MathUtils.nextGaussian() / Math.sqrt(mass.get(i));
            }

            if (mask != null) {
                applyMask(velocity);
            }

            return new WrappedVector.Raw(velocity);
        }
    }

    private MinimumTravelInformation getTimeToBoundary(ReadableVector position, ReadableVector velocity) {

        assert (position.getDim() == velocity.getDim());

        int index = -1;
        double minTime = Double.MAX_VALUE;

        for (int i = 0, len = position.getDim(); i < len; ++i) {

            // TODO Here is where we check that x_j > x_i for categorical dimensions

            // TODO I believe we can simply the condition below (for fixed boundaries) with:
            // double travelTime = -position.get(i) / velocity.get(i); // This is only true for boundaries at 0
            // if (travelTime > 0.0 && missingDataMask[positionIndex] == 0.0)

            double travelTime = Math.abs(position.get(i) / velocity.get(i));
            if (travelTime > 0.0 && headingTowardsBinaryBoundary(velocity.get(i), i)) {

                if (travelTime < minTime) {
                    index = i;
                    minTime = travelTime;
                }
            }
        }

        return new MinimumTravelInformation(minTime, index);
    }


    private BounceState doBounce(double remainingTime, double bounceTime,
                                 MinimumTravelInformation boundaryInfo,
                                 WrappedVector position, WrappedVector velocity,
                                 WrappedVector gradient, WrappedVector action,
                                 double[] inertia, double v_Phi_x, double v_Phi_v) {

        double timeToBoundary = boundaryInfo.time;
        int boundaryIndex = boundaryInfo.index[0];
        final BounceState finalBounceState;
        final Type eventType;
        int eventIndex;
        if (remainingTime < Math.min(timeToBoundary, bounceTime)) { // No event during remaining time

            updatePosition(position, velocity, remainingTime);
//            updateGradient(gradient, remainingTime, action);
            updateInertia(inertia, remainingTime, v_Phi_x, v_Phi_v);
            finalBounceState = new BounceState(Type.NONE, -1, 0.0);
        } else {
            if (timeToBoundary < bounceTime) { // Reflect against the boundary
                eventType = Type.BINARY_BOUNDARY;
                eventIndex = boundaryIndex;

                updatePosition(position, velocity, timeToBoundary);
                updateGradient(gradient, timeToBoundary, action);
                position.set(boundaryIndex, 0.0);
                velocity.set(boundaryIndex, -1 * velocity.get(boundaryIndex));
                updateInertia(inertia, timeToBoundary, v_Phi_x, v_Phi_v);


                remainingTime -= timeToBoundary;

            } else { // Bounce caused by the gradient
                eventType = Type.GRADIENT;
                eventIndex = -1;
                updatePosition(position, velocity, bounceTime);

                updateGradient(gradient, bounceTime, action);
                updateVelocity(velocity, gradient, preconditioning.mass);

                zeroInertia(inertia);
                remainingTime -= bounceTime;
            }
            finalBounceState = new BounceState(eventType, eventIndex, remainingTime);
        }
        return finalBounceState;
    }

    private void zeroInertia(double[] inertia){
        inertia[0] = 0.0;
    }

    private void updateInertia(double[] inertia, double time, double v_Phi_x, double  v_Phi_v) {
        inertia[0] = inertia[0] - (time*time/2*v_Phi_v+time*v_Phi_x);
    }

    private static void updateVelocity(WrappedVector velocity, WrappedVector gradient, ReadableVector mass) {

        ReadableVector gDivM = new ReadableVector.Quotient(gradient, mass);

        double vg = innerProduct(velocity, gradient);
        double ggDivM = innerProduct(gradient, gDivM);

        for (int i = 0, len = velocity.getDim(); i < len; ++i) {
            velocity.set(i, velocity.get(i) - 2 * vg / ggDivM * gDivM.get(i));
        }
    }

    @Override
    public String getOperatorName() {
        return "Hamiltonian bouncy particle operator";
    }

    @Override
    double integrateTrajectory(WrappedVector position, WrappedVector momentum) {return 0.0;}

    double integrateTrajectory(WrappedVector position, WrappedVector velocity, WrappedVector inertia) {
        WrappedVector gradient = getInitialGradient();
        WrappedVector action = new WrappedVector.Raw(new double[] {0.0}); // initialize to empty array since it gets computed on line 219
        BounceState bounceState = new BounceState(drawTotalTravelTime());

        initializeNumEvent();

        while (bounceState.remainingTime > 0) {

            if (bounceState.type == Type.BINARY_BOUNDARY) {
                updateAction(action, velocity, bounceState.index);
            } else {
                action = getPrecisionProduct(velocity);
            }

            double v_Phi_x = -innerProduct(velocity, gradient);
            double v_Phi_v = innerProduct(velocity, action);

            double bounceTime = getBounceTime(v_Phi_v, v_Phi_x, inertia.getBuffer());
            MinimumTravelInformation travelInfo = getTimeToBoundary(position, velocity);

            bounceState = doBounce(
                    bounceState.remainingTime, bounceTime, travelInfo,
                    position, velocity, gradient, action, inertia.getBuffer(), v_Phi_x, v_Phi_v
            );

            recordOneMoreEvent();
        }
        // storeVelocity(velocity);
        return 0.0;
    }


    @Override
    public void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum, WrappedVector inertia, WrappedVector gradient,
                                                 int direction, double time) {

        preconditioning.totalTravelTime = time;
        if (direction == -1) {
            // negate momentum
            negateVector(momentum);
        }
        // integrate
        integrateTrajectory(position, momentum, inertia);
        if (direction == -1) {
            //negate momentum again
            negateVector(momentum);
        }
        ReadableVector.Utils.setParameter(position, parameter);
    }

    @Override
    public void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum, WrappedVector gradient,
                                                 int direction, double time){}
    @Override
    public void providerUpdatePreconditioning() {
        updatePreconditioning(new WrappedVector.Raw(this.getInitialPosition()));
    }

    @Override
    public double[] getInitialPosition() {
        return parameter.getParameterValues();
    }

    @Override
    public double getParameterLogJacobian() { // transform is not allowed yet.
        return 0;
    }

    @Override
    public int getNumGradientEvent() {
        return numGradientEvents;
    }

    @Override
    public int getNumBoundaryEvent() {
        return numBoundaryEvents;
    }

    @Override
    public double[] getMask() {
        return maskVector;
    }

    @Override
    public Transform getTransform() {
        return null;
    }

    @Override
    public GradientWrtParameterProvider getGradientProvider() {
        return gradientProvider;
    }

    @Override
    public void setParameter(double[] position) {
        ReadableVector.Utils.setParameter(position, parameter);
    }

    @Override
    public WrappedVector drawInitialMomentum() {
        return drawInitialVelocity();
    }

    @Override
    public WrappedVector drawMomentum() {
        return drawInitialVelocity();
    }


    @Override
    public double getJointProbability(WrappedVector momentum) {
        return gradientProvider.getLikelihood().getLogLikelihood() - getKineticEnergy(momentum) - getParameterLogJacobian();
    }

    @Override
    public double getJointProbability(WrappedVector momentum, WrappedVector inertia) {
        return gradientProvider.getLikelihood().getLogLikelihood() - getKineticEnergy(momentum) - getParameterLogJacobian() - inertia.getBuffer()[0];
    }

    @Override
    public double getLogLikelihood() {
        return gradientProvider.getLikelihood().getLogLikelihood();
    }

    @Override
    public double getKineticEnergy(ReadableVector momentum) {

        final int dim = momentum.getDim();

        double energy = 0.0;
        for (int i = 0; i < dim; i++) {
            energy += momentum.get(i) * momentum.get(i);
        }
        return energy / 2.0;
    }
    @Override
    public double getStepSize() {
        return preconditioning.totalTravelTime;
    }

    private void negateVector(WrappedVector vector) {
        for (int i = 0; i < vector.getDim(); i++) {
            vector.set(i, -vector.get(i));
        }
    }
}
