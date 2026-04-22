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

# sbatch --array=666-670 --nodelist=compute-124 --export=output=hbps_0-5_5M,xml=24t_HBPS_0-5.xml --job-name=hiv0-5 --mem=15G phylogenetic.sh
# sbatch --array=666-670 --nodelist=compute-124 --export=output=hbps_1_2-5M,xml=24t_HBPS_1.xml --job-name=hiv1 --mem=15G phylogenetic.sh
# sbatch --array=666-670 --nodelist=compute-124 --export=output=bps_0-5_5M,xml=24t_BPS.xml --job-name=hivbps --mem=15G phylogenetic.sh
# sbatch --array=666-670 --nodelist=compute-124 --export=output=split_50K,xml=24t_HNUTS_split.xml --job-name=hivsplit --mem=15G phylogenetic.sh