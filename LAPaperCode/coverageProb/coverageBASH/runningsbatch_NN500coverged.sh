#!/bin/bash
#SBATCH --job-name=NN500covergedBASH
#SBATCH --mail-type=ALL
#SBATCH --mail-user=george.glen@ufl.edu
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=100gb
#SBATCH --time=20-00:00:00
#SBATCH --output=NN500covergedBASH_%j.log
#SBATCH --qos=epi
#SBATCH --account=epi

module load R/4.2
Rscript /home/george.glen/cmd/TMBPaper_MCcoverage/NN500covergedBASH.R
