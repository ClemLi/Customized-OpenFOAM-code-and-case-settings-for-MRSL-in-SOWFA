#!/bin/sh


#SBATCH --job-name=pre2_2NBL
#SBATCH --partition=genoa
#SBATCH --time=120:00:00
#SBATCH --ntasks=384
#SBATCH --output JobStatistics_%j.log
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=
#SBATCH --output=logs3_%j.log
#SBATCH --nodes=2
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
startTime=36000                # Start time (DO NOT change for restarted runs)
writeIntervalBeforeBdData=1800 # Averages computation starts with "BdData" times
endTimeBeforeBdData=$startTime 
startTimeBdData=$endTimeBeforeBdData
writeIntervalBdData=1800
endTimeBdData=41400
startAveraging=37800
timeStepSize=0.36
sliceDeltaT=0.36
solver=superDeliciousVanilla
initializer=setFieldsABL

stabState=neutral

cp -r ./../../precursor_stage_2_1/precursor_stage_2_1_NBL/$startTime.essentials ./.
cp -r ./../precursor_stage_2_2_Master/* ./.


rm -rf setUp
cp setUp.$stabState                                    setUp

cp -rf $startTime.essentials                           $startTime
# cp     $startTime/qwall.$stabState  $startTime/qwall
# rm     $startTime/qwall.*
# **************************************************************************************************** #

# ***************************** COPY APPROPRIATE FILES AND SET VARIABLES ***************************** #
cp -rf system/controlDict.$solver                          system/controlDict
# **************************************************************************************************** #

echo "done basic stuff"

check_success() {
    if [ $? -ne 0 ]; then
        echo "Command failed. Exiting."
        exit 1
    fi
}

blockMesh -noFunctionObjects > log.1.blockMesh 2>&1
check_success


# changeDictionary -dict system/changeDictionaryDict.updateBCs.qwall  -subDict dictionaryReplacement -time $startTime -enableFunctionEntries > log.changeDictionary.updateBCs.qwall.1 2>&1


# Run decomposePar if domain not yet decomposed
decomposePar -cellDist -force -time $startTime > log.2.decomposePar 2>&1
latestTime=$startTime


echo "now bddata"
# Split run to save boundaryData and averages
# continueTime=$(( $latestTime > $endTimeBeforeBdData ? $latestTime : $startTimeBdData ))
continueTime=$startTime


mpirun -np $cores $solver -parallel > log.3.$solver.saveBdData.startAt$continueTime 2>&1

check_success

echo "Ending OpenFOAM job OF ABL solver at: " $(date)

mv foam2run_${SLURM_JOBID}.log foam2run_startAt$continueTime


sbatch onlyRecon_s2.bash

sbatch submittingPrepTurb.bash



# **************************************************************************************************** #
