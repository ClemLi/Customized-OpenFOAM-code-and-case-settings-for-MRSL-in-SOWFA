#!/bin/bash

module purge

# Prevent old GCCcore 10.3.0 libstdc++ from hijacking runtime
export LD_LIBRARY_PATH=$(echo "${LD_LIBRARY_PATH:-}" | tr ':' '\n' \
  | grep -v '/home/ndangi/.local/easybuild' \
  | paste -sd ':' -)

module load 2025
module load SciPy-bundle/2025.06-gfbf-2025a

# source $HOME/.bash_profile
# cores=$SLURM_NTASKS
# source /home/ndangi/cobra_ENV2.bash
# SOWFA6env

# export OMPI_MCA_pml=ob1
# export OMPI_MCA_btl=self,tcp
# export OMPI_MCA_btl_tcp_if_include=ib0
# export OMPI_MCA_oob_tcp_if_include=ib0
# export WS_HUB=6.5


module load 2025
module load SciPy-bundle/2025.06-gfbf-2025a


# Prepare for turbine simulations
mkdir drivingData
cp -r constant drivingData/.
mv postProcessing/boundaryData drivingData/

cd postProcessing

# cp /home/ndangi/sourceHistoryRead.py .
# cp /home/ndangi/sourceData.py .

cp ./../sourceHistoryRead.py .
cp ./../sourceData.py .

# unset PYTHONPATH
# export PYTHONNOUSERSITE=1

# Run the python script
python sourceHistoryRead.py

# Return to the parent directory
cd ../


