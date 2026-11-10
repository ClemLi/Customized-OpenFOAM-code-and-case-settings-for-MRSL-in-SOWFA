#!/bin/sh

#SBATCH --job-name=SOWFA_pre
#SBATCH --partition=genoa
#SBATCH --time=01:00:00
#SBATCH --ntasks=192
#SBATCH --output JobStatistics_%j.log
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=
#SBATCH --output=logs_pre_%j.log
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

#echo "Working directory is" $SLURM_SUBMIT_DIR
#echo "Job name is" $SLURM_JOB_NAME
#echo "Submit time is" $(squeue -u $USER -o '%30j %20V' | grep --e $SLURM_JOB_NAME | awk '{print $2}')
#echo "Starting OpenFOAM job at: " $(date)
#echo "using" $cores "core(s)"


cp -r ./../precursor_stage_2_1_Master/* ./.
cp -r ./../../precursor_stage_1/precursor_stage_1_NBL/21600.essentials ./21600.essentials

# ******************************************** USER INPUT ******************************************** #
startTime=21600                     # Start time
stabState=neutral               # neutral, stable, unstable, wrf
solver=superDeliciousVanilla
# **************************************************************************************************** #

# ************************************** COPY APPROPRIATE FILES ************************************** #
cp system/controlDict.$solver                          system/controlDict
cp setUp.$stabState                                    setUp
# **************************************************************************************************** #



# Function to check command success
check_success() {
    if [ $? -ne 0 ]; then
        echo "Command failed. Exiting."
        exit 1
    fi
}

# Build the mesh and check it
blockMesh -noFunctionObjects > log.1.blockMesh 2>&1
check_success

#checkMesh -noFunctionObjects > log.1.checkMesh.background 2>&1


echo "Ending OpenFOAM job of blockmesh at: " $(date)

# Check if job 1 completed successfully
if [ $? -eq 0 ]; then
 echo "going to job 2"
 sbatch precursorStage2Run.bash  # Submit job 2 upon successful completion of job 1
else
 echo "Job 1 encountered an error. Exiting without running job 2."
fi


# **************************************************************************************************** #
