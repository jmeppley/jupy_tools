#!/bin/bash
#SBATCH -N 1               # 1 node
#SBATCH -c 40
#SBATCH --export=PATH
#SBATCH -o slurm/snake.%A.out 
#SBATCH -e slurm/snake.%A.out 

J=$SLURM_CPUS_PER_TASK

echo -n "SNAKEMAKE VERSION: "
snakemake --version
echo " "
echo "snakemake -c $J $@"
echo " "

snakemake -c $J $@
