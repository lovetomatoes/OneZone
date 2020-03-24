#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH --partition=C032M0128G
#SBATCH --qos=high
#SBATCH -A hpc1806187118
#SBATCH -J PYevol_wli
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00

module load gsl/2.4 anaconda/3.7.1

python evol.py 
