#!/bin/sh


#SBATCH --job-name=SOWFA_s4
#SBATCH --partition=genoa
#SBATCH --time=00:30:00
#SBATCH --ntasks=48
#SBATCH --output JobStatistics_%j.log
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=
#SBATCH --output=logs4_%j.log


. ./prepTurb.bash