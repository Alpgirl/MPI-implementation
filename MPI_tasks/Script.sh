#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --job-name=main
#SBATCH --output=out.txt
#SBATCH --error=error
#SBATCH --time=1:00:00

mpirun ./integral
