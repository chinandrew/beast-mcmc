# Code for Section 5.3 of "MCMC using bouncy Hamiltonian dynamics: A unifying framework for Hamiltonian Monte Carlo and piecewise deterministic Markov process samplers" by Chin et al. 2024

## Setting up BEAST
This code is a fork of the [BEAST](https://github.com/beast-dev/beast-mcmc) software, modified with the Hamiltonian Bouncy Particle Sampler of Chin et al. 

1. Clone the `hbps_develop` branch: `git clone -b hbps_develop git@github.com:chinandrew/beast-mcmc.git`
2. Enter the cloned directory and build the software with Apache ant: `cd beast-mcmc; ant`

A `beast.jar` file should now be located in `build/dist`, which is run in the next step.

## Running BEAST
BEAST can be run with the following command:
```
java -jar ${path_to_beast_jar} -seed ${seed} ${path_to_input_xml}
```
where `path_to_beast_jar` and `path_to_input_xml` specify paths to the `beast.jar` and XML file and `seed` specifies the random seed.
XML files can be found in `simulation/xml`, with the three XML files `24t_HBPS.xml`, `24t_HBPS_split.xml`, and `24t_BPS.xml`.

Running BEAST will output a file `corr_24t_HNUTS.log` containing the samples, of which the `precisionMatrix` columns are the ones of interest.

