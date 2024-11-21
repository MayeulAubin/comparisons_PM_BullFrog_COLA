#!/bin/sh
#SBATCH --job-name=N1536z24
#SBATCH -o /home/aubin/data/BullFrog/comparison_article/logs/N1536z24.out
#SBATCH -e /home/aubin/data/BullFrog/comparison_article/logs/N1536z24.err
#SBATCH -N 1                        # Number of nodes (value or min-max)
#SBATCH -n 128                      # The number of tasks (i.e. cores) per node
#SBATCH --time=08:00:00

##SBATCH --exclusive
##SBATCH --mem=64G                  # Memory pool for all cores (see also --mem-per-cpu)
##SBATCH --array=0-10               # Size of the array
##SBATCH --gres=gpu:1               # Number of GPU
##SBATCH --constraint=i1            # Constraint e.g. specific node

# conda activate simbel2

export OMP_NUM_THREADS=128
python3 "/home/aubin/BullFrog/comparisons_PM_BullFrog_COLA/compare_PM_COLA_BullFrog_parser.py" \
    --workdir /home/aubin/data/BullFrog/comparison_article/params/ \
    --simdir_root /home/aubin/data/BullFrog/comparison_article/sims/ \
    --name N1536z24 \
    -pm 200  \
    -cola 1 4 10 50 \
    -bf 1 4 10 50 \
    -L 2000 \
    --size 768 \
    --Np 768 \
    --Npm 1536 \
    --RedshiftLPT 24.0 \
    --RedshiftFC 0.0 \
    --verbosity 1 \
    --force False

exit 0
