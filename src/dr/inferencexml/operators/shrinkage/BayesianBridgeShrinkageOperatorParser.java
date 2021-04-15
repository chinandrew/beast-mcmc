package dr.inferencexml.operators.shrinkage;

import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.shrinkage.BayesianBridgeShrinkageOperator;
import dr.math.distributions.GammaDistribution;
import dr.xml.*;

import static dr.evoxml.MaskedPatternsParser.MASK;
import static dr.inference.operators.MCMCOperator.WEIGHT;

public class BayesianBridgeShrinkageOperatorParser extends AbstractXMLObjectParser {

    public final static String BAYESIAN_BRIDGE_PARSER = "bayesianBridgeGibbsOperator";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(WEIGHT);

        BayesianBridgeStatisticsProvider bayesianBridge =
                (BayesianBridgeStatisticsProvider) xo.getChild(BayesianBridgeStatisticsProvider.class);

        if (bayesianBridge == null) {
            bayesianBridge = parseAutoCorrelatedRates(xo);
        }

        GammaDistribution globalScalePrior = null;

        DistributionLikelihood prior = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
        if (prior != null) {
            if (prior.getDistribution() instanceof GammaDistribution) {
                globalScalePrior = (GammaDistribution) prior.getDistribution();
            } else {
                throw new XMLParseException("Gibbs sampler only implemented for a gamma distributed global scale");
            }
        }

        Parameter mask = null;
        if (xo.hasChildNamed(MASK)) {
            mask = (Parameter) xo.getElementFirstChild(MASK);
        }

        return new BayesianBridgeShrinkageOperator(bayesianBridge, globalScalePrior, mask, weight);
    }

    private BayesianBridgeStatisticsProvider parseAutoCorrelatedRates(XMLObject xo) throws XMLParseException {

        final AutoCorrelatedBranchRatesDistribution rates =
                (AutoCorrelatedBranchRatesDistribution) xo.getChild(AutoCorrelatedBranchRatesDistribution.class);

        if (!(rates.getPrior() instanceof BayesianBridgeDistributionModel)) {
            throw new XMLParseException("Gibbs sampler only implemented for a Bayesian Bridge prior");
        }

        final BayesianBridgeDistributionModel prior = (BayesianBridgeDistributionModel) rates.getPrior();

        return new BayesianBridgeStatisticsProvider() {

            @Override
            public double getCoefficient(int i) { return rates.getIncrement(i); }

            @Override
            public Parameter getGlobalScale() {
                return prior.getGlobalScale();
            }

            @Override
            public Parameter getLocalScale() {
                return prior.getLocalScale();
            }

            @Override
            public Parameter getExponent() {
                return prior.getExponent();
            }

            @Override
            public int getDimension() {
                return rates.getDimension();
            }
        };
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            new XORRule(
                    new ElementRule(BayesianBridgeStatisticsProvider.class),
                    new ElementRule(AutoCorrelatedBranchRatesDistribution.class)
            ),
            new ElementRule(DistributionLikelihood.class, true),
            new ElementRule(MASK, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),

            }, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return BayesianBridgeShrinkageOperator.class;
    }

    @Override
    public String getParserName() {
        return BAYESIAN_BRIDGE_PARSER;
    }
}
