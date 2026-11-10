#!/bin/sh
#

#SBATCH --job-name=preMRSL
#SBATCH --partition=genoa
#SBATCH --time=24:00:00
#SBATCH --ntasks=192
#SBATCH --output JobStatistics_%j.log
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=
#SBATCH --output=logs2_%j.log
#SBATCH --nodes=1
#SBATCH --exclusive

#Change ntasks and mem for final job, this is for trail.

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
echo "Running on $SLURM_NNODES nodes."
# ******************************************** USER INPUT ******************************************** #
stabState=neutral

cp -r ./../suscessor_Master/* ./.
export precursorDir='./../../precursor_stage_2_2/precursor_stage_2_2_NBL'


rm setUp
cp setUp.$stabState                                    setUp



export startTime=36000                 # Start time
export endTime=43200                 # Start time
export updateBCType=1                  # Boolean for whether or not the boundary condition types will be updated over 
                                #    what is in the initial conditions files. Leave it 0 for precursors
export inflowDir='southWest'    # For inflow/outflow cases, specify the inflow direction.  Choices are 'west',
                                #   'east', 'south', 'west', 'southWest', 'northWest', 'southEast', and
                                #   'northEast'.  There is a 'cyclic' option too in case you need to change
                                #    back to cyclic or simply update the number of boundary face entries.

export parallel=0
export cores=1
export refinementLevels=0 # CHANGE THIS; for trial run its 0 to save time. typically on basis of epsilon>=2*grid size


# Function to check command success
check_success() {
    if [ $? -ne 0 ]; then
        echo "Command failed. Exiting."
        exit 1
    fi
}


refineMeshLocal()
{
   i=1
   while [ $i -le $1 ]
   do
      echo "   -Performing level $i local refinement with topoSet/refineHexMesh"
      echo "      *selecting cells to refine..."
      topoSet -dict system/topoSetDict.local.$i > log.topoSet.local.$i 2>&1

      echo "      *refining cells..."
      refineHexMesh local -overwrite > log.refineHexMesh.local.$i 2>&1

      let i=i+1
   done
}

# Global refinement performed in parallel.
refineMeshGlobal()
{
   i=1
   while [ $i -le $1 ]
   do
      echo "   -Performing level $i global refinement with refineMesh"
      echo "      *refining cells..."
      mpirun -np $cores refineMesh -all -parallel -overwrite > log.refineMesh.global$

      let i=i+1
   done
}



cp -rf system/controlDict.superDeliciousVanilla            system/controlDict




rm -rf constant/polyMesh
# Build the mesh and check it
blockMesh -noFunctionObjects > log.1.blockMesh 2>&1
check_success

# refineMeshLocal $refinementLevels
# check_success

mapFields $precursorDir -sourceTime 36000 -mapMethod cellPointInterpolate -consistent > log.mapFields 2>&1
check_success


cp -rf $precursorDir/postProcessing/givenSourceT constant/
cp -rf $precursorDir/postProcessing/givenSourceU constant/ 

topoSet    > log.1.topoSet 2>&1

cd constant
ln -s $precursorDir/drivingData/boundaryData ./boundaryData
ln -s $precursorDir/postProcessing/sources ./sources
cp $precursorDir/postProcessing/givenSource* .
cd ../

changeDictionary -dict system/changeDictionaryDict.updateBCs.$inflowDir  -subDict dictionaryReplacement -time $startTime -enableFunctionEntries > log.changeDictionary.updateBCs.$inflowDir.1 2>&1


renumberMesh -overwrite > log.renumberMesh 2>&1

decomposePar -cellDist -force > log.decomposePar 2>&1
check_success


sbatch runturb.bash

