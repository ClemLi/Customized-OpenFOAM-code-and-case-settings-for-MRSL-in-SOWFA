#!/bin/sh


#SBATCH --job-name=pre_stage1
#SBATCH --partition=genoa
#SBATCH --time=120:00:00
#SBATCH --ntasks=384
#SBATCH --output JobStatistics_%j.log
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=Y.Li-18@tudelft.nl
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
startTime=0                     # Start time (DO NOT change for restarted runs)
writeIntervalBeforeBdData=1800 # Averages computation starts with "BdData" times
endTimeBeforeBdData=21600 
startTimeBdData=$endTimeBeforeBdData
writeIntervalBdData=400
endTimeBdData=21600
timeStepSize=0.36
solver=superDeliciousVanilla
initializer=setFieldsABL

stabState=neutral
cp system/controlDict.$solver                          system/controlDict
rm setUp
cp setUp.$stabState                                    setUp
rm -rf $startTime
cp -rf $startTime.original $startTime
cp -rf $startTime/qwall.$stabState  $startTime/qwall
rm     $startTime/qwall.*
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
echo "going todecomp"

# Run decomposePar if domain not yet decomposed
if [ ! -d "processor0" ]; then
    decomposePar -cellDist -force > log.2.decomposePar 2>&1
    latestTime=$startTime
else
    latestTime=$(foamListTimes -processor -latestTime -withZero -noFunctionObjects | tail -1)
fi

check_success

echo "time for setfieldsabl"
# Run the flow field initializer
if [ $initializer = setFieldsABL ] && [ ! -f log.2.$initializer ];  then
        # srun --mpi=pmix -n $cores --cpu_bind=cores $initializer -parallel
    # mpirun --oversubscribe -np $cores $initializer -parallel > log.2.$initializer 2>&1
        mpirun -np $cores $initializer -parallel > log.2.$initializer 2>&1
fi

check_success
echo "going to super"
# Split run to get to developed-flow stage
if [ $latestTime -lt $endTimeBeforeBdData ]; then
    mpirun -np $cores $solver -parallel > log.2.$solver.startAt$latestTime 2>&1
fi

check_success
echo "going to recons"
# Make the precursor data ready for future mapping 
if [ -d processor0/$startTimeBdData ] && [ ! -f log.3.reconstructPar ]; then
   foamDictionary -entry "writeFormat" -set ascii -disableFunctionEntries system/controlDict
    reconstructPar -time $startTimeBdData -fields '(U T k p_rgh kappat nut qwall Rwall)'> log.3.reconstructPar 2>&1
    # The z0 gets recomposed as 0, sometimes with '-nan` values. This is a known bug. Workaround below
    z0=$(foamDictionary -entry "z0" -value setUp)
    sed -i 's/-nan/0/g' $startTimeBdData/Rwall
    foamDictionary -entry "boundaryField.lower.z0" -set $z0 $startTimeBdData/Rwall
    foamDictionary -entry "writeFormat" -set ascii -disableFunctionEntries system/controlDict
else
    touch "WARNING_4"
    echo "Case not reconstructed and first time of boundaryData not fixed. Check if endTime dir exists."
fi

mv 21600 ./21600.essentials

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
