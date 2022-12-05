/*
 * NewHamiltonianMonteCarloOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

import static dr.math.matrixAlgebra.ReadableVector.Utils.innerProduct;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class BouncyParticleOperator extends AbstractParticleOperator implements Loggable {

    public BouncyParticleOperator(GradientWrtParameterProvider gradientProvider,
                                             PrecisionMatrixVectorProductProvider multiplicationProvider,
                                             PrecisionColumnProvider columnProvider,
                                             double weight, Options runtimeOptions, NativeCodeOptions nativeOptions,
                                             boolean refreshVelocity, Parameter mask,
                                             MassPreconditioner massPreconditioner,
                                             MassPreconditionScheduler.Type preconditionSchedulerType) {
        super(gradientProvider, multiplicationProvider, columnProvider, weight, runtimeOptions, nativeOptions,
                refreshVelocity, mask, null, massPreconditioner, preconditionSchedulerType);
    }

    @Override
    public String getOperatorName() {
        return "Hamiltonian bouncy particle operator";
    }

    @Override
    double integrateTrajectory(WrappedVector position, WrappedVector momentum) {

        WrappedVector velocity = drawInitialVelocity();
        WrappedVector gradient = getInitialGradient();
        WrappedVector action = getPrecisionProduct(velocity);
        double[] inertia = new double[1]; // array so retains state
        inertia[0] = drawInitialInertia();
        BounceState bounceState = new BounceState(drawTotalTravelTime());

        initializeNumEvent();


        double x_Phi_x_old;
        double x_Phi_x_new;

        while (bounceState.remainingTime > 0) {

            if (bounceState.type == Type.BINARY_BOUNDARY) {
                updateAction(action, velocity, bounceState.index);
            } else {
                action = getPrecisionProduct(velocity);
            }
            x_Phi_x_old =  -innerProduct(position, gradient);
            double v_Phi_x = -innerProduct(velocity, gradient);
            double v_Phi_v = innerProduct(velocity, action);
            System.out.println( (-v_Phi_x + Math.sqrt(v_Phi_x * v_Phi_x + v_Phi_v * 2 * inertia[0])) / v_Phi_v);
            // (-xprecv + (xprecv ** 2 + vprecv * 2 * inertia)**0.5) / vprecv
            double bounceTime = getBounceTime(v_Phi_v, v_Phi_x, inertia);

//            for (int i = 0, len = 85; i < len; ++i) {
//                updatePosition(position, velocity, 0.0001);
//                updateGradient(gradient, 0.0001, action);
//                x_Phi_x_new =  -innerProduct(position, gradient);
//                System.out.println(inertia[0] - x_Phi_x_new / 2.0 + x_Phi_x_old / 2.0);
//                System.out.println((i+1)*0.0001);
//            }
            double actualBounceTime = 0.0085190844;
            updatePosition(position, velocity, actualBounceTime);
            updateGradient(gradient, actualBounceTime, action);
            x_Phi_x_new =  -innerProduct(position, gradient);

            System.out.println("U diff");
            System.out.println(v_Phi_x*actualBounceTime + v_Phi_v*actualBounceTime*actualBounceTime/2);
            System.out.println(- x_Phi_x_new / 2.0 + x_Phi_x_old / 2.0);
            System.out.println("inertia");
            System.out.println(inertia[0]);
            System.out.println("computed bouncetime ##################################");
            System.out.println(bounceTime);
            System.out.println("should be near zero");
            System.out.println(inertia[0] - x_Phi_x_new / 2.0 + x_Phi_x_old / 2.0);


            MinimumTravelInformation travelInfo = getTimeToBoundary(position, velocity);

            if (printEventLocations) System.err.println(position);
            bounceState = doBounce(
                    bounceState.remainingTime, bounceTime, travelInfo,
                    position, velocity, gradient, action, inertia
            );

            recordOneMoreEvent();
        }
        storeVelocity(velocity);
        return 0.0;
    }

    private BounceState doBounce(double remainingTime, double bounceTime,
                                 MinimumTravelInformation boundaryInfo,
                                 WrappedVector position, WrappedVector velocity,
                                 WrappedVector gradient, WrappedVector action,
                                 double[] inertia) {

        double x_Phi_x_old;
        double x_Phi_x_new;
        double timeToBoundary = boundaryInfo.time;
        int boundaryIndex = boundaryInfo.index[0];
        final BounceState finalBounceState;
        final Type eventType;
        int eventIndex;
        if (remainingTime < Math.min(timeToBoundary, bounceTime)) { // No event during remaining time

            updatePosition(position, velocity, remainingTime);
            finalBounceState = new BounceState(Type.NONE, -1, 0.0);
        } else {
            if (timeToBoundary < bounceTime) { // Reflect against the boundary
                eventType = Type.BINARY_BOUNDARY;
                eventIndex = boundaryIndex;
                x_Phi_x_old =  -innerProduct(position, gradient);

                updatePosition(position, velocity, timeToBoundary);
                updateGradient(gradient, timeToBoundary, action);
                position.set(boundaryIndex, 0.0);
                velocity.set(boundaryIndex, -1 * velocity.get(boundaryIndex));
                x_Phi_x_new =  -innerProduct(position, gradient);
                updateInertia(inertia, x_Phi_x_old, x_Phi_x_new);


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

    private void updateInertia(double[] inertia, double x_Phi_x_old, double  x_Phi_x_new) {
        inertia[0] = inertia[0] - x_Phi_x_new / 2.0 + x_Phi_x_old / 2.0;
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

    private double drawInitialInertia() {
        return MathUtils.nextExponential(1);
    }

    @SuppressWarnings("all")
    private double getBounceTime(double v_phi_v, double v_phi_x, double[] inertia) {
        double a = v_phi_v;
        double b = 2.0 * v_phi_x;
        double c = - 2.0 * inertia[0];
//        System.out.println( (-v_phi_x + Math.sqrt(v_phi_x *v_phi_x + v_phi_v * 2 * inertia[0])) / v_phi_v);
        return (-b + Math.sqrt(b * b - 4 * a * c)) / 2 / a;
    }

    private static void updateVelocity(WrappedVector velocity, WrappedVector gradient, ReadableVector mass) {

        ReadableVector gDivM = new ReadableVector.Quotient(gradient, mass);

        double vg = innerProduct(velocity, gradient);
        double ggDivM = innerProduct(gradient, gDivM);

        for (int i = 0, len = velocity.getDim(); i < len; ++i) {
            velocity.set(i, velocity.get(i) - 2 * vg / ggDivM * gDivM.get(i));
        }
    }

    private WrappedVector storedVelocity;

    @Override
    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[preconditioning.mass.getDim()];
        for (int i = 0; i < preconditioning.mass.getDim(); ++i) {
            final int index = i;
            columns[i] = new NumberColumn("v" + index) {
                @Override
                public double getDoubleValue() {
                    if (storedVelocity != null) {
                        return storedVelocity.get(index);
                    } else {
                        return 0.0;
                    }
                }
            };
        }

        return columns;
    }

    private final static boolean printEventLocations = false;
}
