package dr.inference.operators;

import dr.inference.model.GeneralParameterBounds;
import dr.inference.model.Parameter;
import dr.inference.model.TransformedParameter;

public class TransformedParameterOperator extends AbstractAdaptableOperator {
    private boolean isAdaptable;
    private SimpleMCMCOperator subOperator;
    private TransformedParameter parameter;
    private boolean checkValid;
    private GeneralParameterBounds generalBounds;

    public TransformedParameterOperator(SimpleMCMCOperator operator, GeneralParameterBounds generalBounds) {

        this.subOperator = operator;
        setWeight(operator.getWeight());
        this.isAdaptable = operator instanceof AbstractAdaptableOperator;
        this.parameter = (TransformedParameter) operator.getParameter();

        this.generalBounds = generalBounds;
        this.checkValid = generalBounds != null;
    }


    @Override
    protected void setAdaptableParameterValue(double value) {
        if (isAdaptable) {
            ((AbstractAdaptableOperator) subOperator).setAdaptableParameterValue(value);
        }
    }

    @Override
    protected double getAdaptableParameterValue() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getAdaptableParameterValue();
        }
        return 0;
    }

    @Override
    public double getRawParameter() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getRawParameter();
        }
        throw new RuntimeException("not actually adaptable parameter");
    }

    @Override
    public String getAdaptableParameterName() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getAdaptableParameterName();
        }
        throw new RuntimeException("not actually adaptable parameter");
    }

    @Override
    public String getOperatorName() {
        return "transformedParameterOperator." + subOperator.getOperatorName();
    }

    @Override
    public double doOperation() {
        double[] oldValues = parameter.getParameterUntransformedValues();
        double ratio = subOperator.doOperation();
        double[] newValues = parameter.getParameterUntransformedValues();


        if (checkValid) { // GH: below is sloppy, but best I could do without refactoring how Parameter handles bounds
            if (generalBounds == null && !parameter.isWithinBounds()) {
                return Double.NEGATIVE_INFINITY;
            } else if (!generalBounds.satisfiesBounds(parameter)) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // Compute Jacobians
        ratio += parameter.diffLogJacobian(oldValues, newValues);

        return ratio;
    }

    @Override
    public Parameter getParameter() {
        return subOperator.getParameter();
    }
}
