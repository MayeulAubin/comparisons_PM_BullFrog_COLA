#!/bin/bash
#SBATCH -J Compare_PM_COLA_BF        # A single job name for the array
#SBATCH --array=0-10                 # Size of the array
#SBATCH -p                           # The partition to use
#SBATCH -N 1                         # The number of nodes to use (min-max)
#SBATCH -n 64                        # The number of tasks (i.e. cores)
##SBATCH --gres=gpu:1                # Number of GPU, here commented out
#SBATCH -t 0-01:00                   # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=64G                    # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ~/slurm_%x_%a_%A.out      # File to which STDOUT will be written, %j inserts jobid, %x the job name, %a the array index, %A array's master job allocation number.
#SBATCH -e ~/slurm_%x_%a_%A.err      # File to which STDERR will be written, %j inserts jobid, %x the job name, %a the array index, %A array's master job allocation number.

conda activate SBMY
python script_compare_PM_COLA_BullFrog ${SLURM_ARRAY_TASK_ID} -pm 1 2 5 10 20 50 100 200 -co 1 2 5 10 20 50 100 -bf 1 2 5 10 20 50 100 -L 2000 -Np 384 -N 384 -Npm 768 -V --RedshiftLPT 49.0
exit
