#!/bin/bash
#SBATCH -N 1               # 1 node
#SBATCH -c 1
#SBATCH -o slurm/snake.%A.out 
#SBATCH -e slurm/snake.%A.out 

# Use this script if your snakeake command has --cluster in it to leave -j unchanged
N=$SLURM_CPUS_PER_TASK

echo snakemake --local-cores $N "$(printf " %q" "$@")"
echo " "

snakemake --local-cores $N "$@"
