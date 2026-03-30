#!/bin/bash
#SBATCH -J brakersnake
#SBATCH -n 10
#SBATCH -t 72:00:00
#SBATCH -p "normal"
#SBATCH --mem=5G
#SBATCH -o logs/test.%A.out
#SBATCH -e logs/test.%A.err

snakemake --cores 10 -j 9 --snakefile workflow/Snakefile --workflow-profile ./profiles/slurm --rerun-incomplete



