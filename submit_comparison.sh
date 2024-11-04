#!/bin/sh
#SBATCH --job-name=run2
#SBATCH -o /data101/hoellinger/BullFrog/comparisons_PM_BullFrog_COLA/tests_parser/logs/run2.log
#SBATCH -e /data101/hoellinger/BullFrog/comparisons_PM_BullFrog_COLA/tests_parser/logs/run2.log
#SBATCH -N 1                        # Number of nodes (value or min-max)
#SBATCH -n 128                      # The number of tasks (i.e. cores) per node
#SBATCH --time=24:00:00

##SBATCH --exclusive
##SBATCH --mem=64G                  # Memory pool for all cores (see also --mem-per-cpu)
##SBATCH --array=0-10               # Size of the array
##SBATCH --gres=gpu:1               # Number of GPU
##SBATCH --constraint=i1            # Constraint e.g. specific node

conda activate simbel2

export OMP_NUM_THREADS=64
python3 "/data101/hoellinger/BullFrog/comparisons_PM_BullFrog_COLA/compare_PM_COLA_BullFrog_parser.py" \
    --workdir /data101/hoellinger/BullFrog/comparisons_PM_BullFrog_COLA/tests_parser/params/ \
    --simdir_root /data101/hoellinger/BullFrog/comparisons_PM_BullFrog_COLA/tests_parser/sims/ \
    --name run2 \
    -pm 50 \
    -cola 1 2 5 10 20 \
    -bf 1 2 5 10 20 \
    -L 2000 \
    --size 512 \
    --Np 512 \
    --Npm 1024 \
    --RedshiftLPT 24.0 \
    --RedshiftFC 0.0 \
    --verbosity 1 \
    --force False

# export OMP_NUM_THREADS=64
# python3 "/data101/hoellinger/BullFrog/comparisons_PM_BullFrog_COLA/compare_PM_COLA_BullFrog_parser.py" \
#     --workdir /data101/hoellinger/BullFrog/comparisons_PM_BullFrog_COLA/tests_parser/params/ \
#     --simdir_root /data101/hoellinger/BullFrog/comparisons_PM_BullFrog_COLA/tests_parser/sims/ \
#     --name run1 \
#     -pm 50 100 \
#     -cola 1 2 5 10 20 50 \
#     -bf 1 2 5 10 20 50 \
#     -L 2000 \
#     --size 384 \
#     --Np 384 \
#     --Npm 384 \
#     --RedshiftLPT 24.0 \
#     --RedshiftFC 0.0 \
#     --verbosity 1 \
#     --force False

exit 0
