#!/bin/bash
#SBATCH -n 20 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3200 # Runtime in minutes
#SBATCH -c 1 # cpus-per-task
#SBATCH -p bigmem # Partition to submit to
# SBATCH --gres=gpu:1 
# SBATCH -p gpu_requeue # Partition to submit to
#SBATCH --mem=200000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o hostname_%j.out # Standard out goes to this file
#SBATCH -e hostname_%j.err # Standard err goes to this filehostname

module load R/4.0.2-fasrc01
Rscript gwas_data_example1_slurm.R