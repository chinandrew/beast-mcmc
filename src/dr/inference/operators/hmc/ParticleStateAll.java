package dr.inference.operators.hmc;
import dr.math.matrixAlgebra.WrappedVector;

public class ParticleStateAll {
    public WrappedVector q;
    public WrappedVector p;
    public WrappedVector inertia;
    public double logp;
    public WrappedVector phi_q;
    public WrappedVector phi_p;
    public double q_phi_p;
    public double p_phi_p;
    public ParticleStateAll(WrappedVector q,
                            WrappedVector p,
                            WrappedVector inertia,
                            double logp,
                            WrappedVector phi_q,
                            WrappedVector phi_p,
                            double q_phi_p,
                            double p_phi_p){
        this.q = q;
        this.p = p;
        this.inertia = inertia;
        this.logp = logp;
        this.phi_q = phi_q;
        this.phi_p = phi_p;
        this.q_phi_p = q_phi_p;
        this.p_phi_p = p_phi_p;
    }
}
