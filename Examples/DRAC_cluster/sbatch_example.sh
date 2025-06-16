#Example of a batch job on DRAC clusters
#To be submitted via 'sbatch sbatch_example' 
#

#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --account=rpp-nahee
#SBATCH --mem=4000M
#SBATCH --output=logs/gen_%a.log
#SBATCH --array=1-10
module --force purge
module load StdEnv/2020 gcc/11.3.0 apptainer scipy-stack/2023b

apptainer exec /cvmfs/software.pacific-neutrino.org/containers/icetray_v1.10 icetray_job.sh $(whoami) ${SLURM_ARRAY_TASK_ID}
