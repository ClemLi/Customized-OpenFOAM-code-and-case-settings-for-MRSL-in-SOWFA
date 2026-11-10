#!/bin/sh


#SBATCH --job-name=pre2_1NBL
#SBATCH --partition=genoa
#SBATCH --time=120:00:00
#SBATCH --ntasks=384
#SBATCH --output JobStatistics_%j.log
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=
#SBATCH --output=logs2_%j.log
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
startTime=21600                # Start time (DO NOT change for restarted runs)
writeIntervalBeforeBdData=1800 # Averages computation starts with "BdData" times
endTimeBeforeBdData=36000 
startTimeBdData=$endTimeBeforeBdData
writeIntervalBdData=400
endTimeBdData=36000
startAveraging=32400
timeStepSize=0.36
sliceDeltaT=72
solver=superDeliciousVanilla
initializer=setFieldsABL

stabState=neutral
rm setUp
cp setUp.$stabState                                    setUp
# rm -rf $startTime
cp -rf $startTime.essentials $startTime
# cp     $startTime/qwall.$stabState  $startTime/qwall
# rm     $startTime/qwall.*
# **************************************************************************************************** #

# ***************************** COPY APPROPRIATE FILES AND SET VARIABLES ***************************** #
cp -rf system/controlDict.$solver                          system/controlDict
# **************************************************************************************************** #


changeDictionary -dict system/changeDictionaryDict.updateBCs.qwall  -subDict dictionaryReplacement -time $startTime -enableFunctionEntries > log.changeDictionary.updateBCs.qwall.1 2>&1


echo "done basic stuff"

check_success() {
    if [ $? -ne 0 ]; then
        echo "Command failed. Exiting."
        exit 1
    fi
}

latestTime=$startTime

# echo "going todecomp"

# Run decomposePar if domain not yet decomposed
decomposePar -cellDist -force -time $startTime > log.2.decomposePar 2>&1
latestTime=$startTime


check_success


check_success
echo "going to super"
# Split run to get to developed-flow stage


mpirun -np $cores $solver -parallel > log.2.$solver.startAt$latestTime 2>&1


check_success
echo "going to recons"
# Make the precursor data ready for future mapping 

foamDictionary -entry "writeFormat" -set ascii -disableFunctionEntries system/controlDict
reconstructPar -time -latestTime -fields '(U T k p_rgh kappat nut qwall Rwall)'> log.3.reconstructPar 2>&1
# The z0 gets recomposed as 0, sometimes with '-nan` values. This is a known bug. Workaround below
z0=$(foamDictionary -entry "z0" -value setUp)
sed -i 's/-nan/0/g' $startTimeBdData/Rwall
foamDictionary -entry "boundaryField.lower.z0" -set $z0 $startTimeBdData/Rwall
foamDictionary -entry "writeFormat" -set ascii -disableFunctionEntries system/controlDict


mv $startTimeBdData ./$startTimeBdData.essentials

check_success
echo "Done recon"

echo "Ending OpenFOAM job OF ABL solver at: " $(date)

mv foam2run_${SLURM_JOBID}.log foam2run_startAt$continueTime

reconstructPar > log.2.reconstructParAllTime 2>&1

z0=$(foamDictionary -entry "z0" -value setUp)
sed -i 's/-nan/0/g' $startTimeBdData/Rwall
foamDictionary -entry "boundaryField.lower.z0" -set $z0 $startTimeBdData/Rwall
foamDictionary -entry "writeFormat" -set ascii -disableFunctionEntries system/controlDict


# **************************************************************************************************** #
