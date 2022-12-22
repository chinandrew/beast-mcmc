package dr.inference.operators.hmc;
import dr.math.matrixAlgebra.WrappedVector;

public class ParticleState {
    public WrappedVector position;
    public WrappedVector momentum;
    public WrappedVector inertia;
    public WrappedVector gradient;
    public ParticleState(WrappedVector position,
                         WrappedVector momentum,
                         WrappedVector inertia,
                         WrappedVector gradient){
        this.position = position;
        this.momentum = momentum;
        this.inertia = inertia;
        this.gradient = gradient;
    }
}
