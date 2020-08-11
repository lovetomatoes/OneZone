#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH --partition=C032M0128G
#SBATCH --qos=low
#SBATCH -A hpc1806187118
#SBATCH -J 10Jc
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks=100
###SBATCH --cpus-per-task=12
#SBATCH --time=10:56:00

module load gsl/2.4
module load gcc/7.2.0
module load mpich/3.2.1-gcc-4.8.5

make
mpirun -n 100 ./main

# mk && m cannot be executed, bash_profile alias only by login shells