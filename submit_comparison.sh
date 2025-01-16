#!/bin/sh
#SBATCH --job-name=N768DKD
#SBATCH -o /home/aubin/data/BullFrog/comparison_article/logs/N768DKD.out
#SBATCH -e /home/aubin/data/BullFrog/comparison_article/logs/N768DKD.err
#SBATCH -N 1                                    # Number of nodes (value or min-max)
#SBATCH -n 64                                   # The number of tasks (i.e. cores) per node
#SBATCH --time=0-12:00:00                         # Runtime in D-HH:MM, min 10 minutes

##SBATCH --exclusive                            # Exclusive node usage
##SBATCH --mem=64G                              # Memory pool for all cores (see also --mem-per-cpu)
##SBATCH --gres=gpu:1                           # Number of GPU
##SBATCH --constraint=i1                        # Constraint e.g. specific node

##SBATCH --array=0-10                           # Size of the array
##SBATCH -o slurm_logs/slurm_%x_%a_%A.out       # File to which STDOUT will be written, %j inserts jobid, %x the job name, %a the array index, %A array's master job allocation number.
##SBATCH -e slurm_logs/slurm_%x_%a_%A.err       # File to which STDERR will be written, %j inserts jobid, %x the job name, %a the array index, %A array's master job allocation number.



export OMP_NUM_THREADS=64
python3 "/home/aubin/BullFrog/comparisons_PM_BullFrog_COLA/compare_PM_COLA_BullFrog_parser.py" \
    --workdir /home/aubin/data/BullFrog/comparison_article/params/ \
    --simdir_root /home/aubin/data/BullFrog/comparison_article/sims/ \
    --name N768DKD \
    -L 2000 \
    --size 384 \
    --Np 384 \
    --Npm 768 \
    --RedshiftLPT 1.0e12 \
    --RedshiftFC 0.0 \
    --verbosity 1 \
    -ts 3 \
    -cola 1 4 10 50 \
    -bf 1 4 10 50 \
    --DKD

exit 0
