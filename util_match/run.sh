#!/bin/bash
#SBATCH --partition=fast        # select partition (normal or fast)
#SBATCH --time=02:00:00         # set time limit in HH:MM:SS
#SBATCH --nodes=1               # number of nodes
#SBATCH --ntasks=1              # number of processes (for MPI)
#SBATCH --cpus-per-task=1       # OMP_NUM_THREADS (for openMP)
#SBATCH --job-name=semiinf      # job name
#SBATCH --output="error.%x.%j"  # standard output and error are redirected to
				# <job name>_<job ID>.out
# for OpenMP jobs
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
# load modules
module purge
module load intel/2013_sp1 openmpi/1.8.8
module load python/3.6.3

# run parallel
srun python3 match_main.py > run.out

