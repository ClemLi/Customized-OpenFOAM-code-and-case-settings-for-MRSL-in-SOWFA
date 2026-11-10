#!/bin/sh

#SBATCH --job-name=rec_s2
#SBATCH --partition=genoa
#SBATCH --time=12:00:00
#SBATCH --ntasks=192
#SBATCH --output JobStatistics_%j.log
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=
#SBATCH --output=logs2_%j.log
#SBATCH --nodes=1
#SBATCH --exclusive


source $HOME/.bash_profile
cores=$SLURM_NTASKS
source /home/ndangi/cobra_ENV2.bash
SOWFA6env

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp
export OMPI_MCA_btl_tcp_if_include=ib0
export OMPI_MCA_oob_tcp_if_include=ib0
export WS_HUB=6.5

echo "Working directory is" $SLURM_SUBMIT_DIR
echo "Job name is" $SLURM_JOB_NAME
echo "Submit time is" $(squeue -u $USER -o '%30j %20V' | grep -e $SLURM_JOB_NAME | awk '{print $2}')
echo "Starting OpenFOAM job at: " $(date)
echo "using" $cores "core(s)"
# ******************************************** USER INPUT ******************************************** #


reconTime=37800

echo "done basic stuff"

check_success() {
    if [ $? -ne 0 ]; then
        echo "Command failed. Exiting."
        exit 1
    fi
}


echo "going to recons"

reconstructPar -time $reconTime > log.reconstruct.$reconTime 2>&1

z0=$(foamDictionary -entry "z0" -value setUp)
sed -i 's/-nan/0/g' $reconTime/Rwall
foamDictionary -entry "boundaryField.lower.z0" -set $z0 $reconTime/Rwall


reconTime=39600

echo "going to recons"

reconstructPar -time $reconTime > log.reconstruct.$reconTime 2>&1

z0=$(foamDictionary -entry "z0" -value setUp)
sed -i 's/-nan/0/g' $reconTime/Rwall
foamDictionary -entry "boundaryField.lower.z0" -set $z0 $reconTime/Rwall



reconTime=41400

echo "going to recons"

reconstructPar -time $reconTime > log.reconstruct.$reconTime 2>&1

z0=$(foamDictionary -entry "z0" -value setUp)
sed -i 's/-nan/0/g' $reconTime/Rwall
foamDictionary -entry "boundaryField.lower.z0" -set $z0 $reconTime/Rwall


# **************************************************************************************************** #
