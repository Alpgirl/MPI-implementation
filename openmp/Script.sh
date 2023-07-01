#!/bin/bash
#
#SBATCH --nodes=1		# Number of nodes
#SBATCH --ntasks-per-node=1	# Number of processes per node
#SBATCH --cpus-per-task=4	# Number of threads per task (process)
#SBATCH --job-name=main
#SBATCH --partition=compute
#SBATCH --comment="Pthreads_OpenMP"
#SBATCH --output=out.txt
#SBATCH --error=error.txt
#SBATCH --time=1:00:00

./integral
