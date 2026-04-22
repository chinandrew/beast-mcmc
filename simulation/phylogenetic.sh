#!/bin/bash
#SBATCH --mail-user=achin23@jhu.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --time=720:00:00

module load conda
source activate env
source activate env
cd /fastscratch/myscratch/achin/hiv_sim
mkdir "${output}_${SLURM_ARRAY_TASK_ID}"
cd "${output}_${SLURM_ARRAY_TASK_ID}"
time java -jar /users/achin/beast-mcmc/build/dist/beast.jar -seed $SLURM_ARRAY_TASK_ID -overwrite /users/achin/beast-mcmc/simulation/xml/$xml

