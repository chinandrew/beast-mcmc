module load conda
source activate env
source activate env
cd hbps/hiv_sim2/
mkdir $output
cd $output
time java -jar /users/achin/beast-mcmc/build/dist/beast.jar -seed $seed -overwrite /users/achin/beast-mcmc/simulation/xml/$xml

