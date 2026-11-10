# Customized-OpenFOAM-code-and-case-settings-for-MRSL-in-SOWFA

Case settings and source code for "Effects of atmospheric thermal stratification on wake aerodynamics of a regenerative wind farm unit". 	
https://doi.org/10.5194/wes-2026-116

################################

Educational and research use of the material is permitted. Please cite the following reference when using or extending the code, figures, or data:

Y. Li, A. Sciacchitano, and W. Yu, 2026 (Accepted by Wind Energy Science),  Effects of atmospheric thermal stratification on wake aerodynamics of a regenerative wind farm unit.

This GitHub repository contains a subset of the "caseSettingsAndCodes" folder from the dataset archived on 4TU.ResearchData. It has been created to facilitate convenient access to and reuse of the code. See https://doi.org/10.4121/ef3ef4c3-1158-4e8b-a872-4432cf124087 for further documentation.

The customized library in "sourceCodeMRSL" can also serve as a general-purpose actuator disk model with force-smearing capabilities. It demonstrates excellent parallel scalability with negligible computational overhead, achieved by restricting the search domains used for velocity sampling and body-force application. Particular care has been taken to ensure that the body forces are correctly distributed, even when the smearing kernel spans multiple processor domains. However, actuator elements should not be placed exactly on computational cell faces or edges, as this may introduce ambiguity in determining which processor owns the actuator element, potentially leading to a segmentation fault during execution.

################################

This repository contains customized OpenFOAM source code libraries (see directory "sourceCodeMRSL"), case setup files (see directory "caseSettings"), and files related to SOWFA-6 (see directory "sowfa_6_related") and OpenFOAM-6 (see directory "openFOAM_6_related"). The files related to SOWFA-6 and OpenFOAM-6 are provided for reference and to document the exact code versions used in this work.

To reproduce the simulations, the customized OpenFOAM libraries must first be compiled by executing the files named "Allwmake" in the "src" directories under "sourceCode". Before compilation, users must specify the output paths for the generated ".so" files in "./src/Make/files".

After compiling the customized libraries and configuring SOWFA-6, the precursor and successor simulations can be executed using the case settings in "caseSettings". Each case directory contains HPC job scripts used in this study. Some modifications to these scripts may be required depending on the local SOWFA-6/OpenFOAM-6 installation and HPC environment. For the current work, the simulations were performed on Snellius (Snellius: de Nationale Supercomputer). For reference, the environment configuration file "cobra_ENV2.bash" is located in "sowfa_6_related".

To start the stage-1 precursor, navigate to "caseSettings/precursor_stage_1/precursor_stage_1_NBL" and submit "preProcess.bash". If successful, the script automatically submits "precursorStage1Run.bash". After completion, a directory named "21600.essentials" should be generated, containing the data required for the stage-2 precursor simulations.

For the stage-2_1 precursor, navigate to any directory under "caseSettings/precursor_stage_2_1" except "precursor_stage_2_1_Master", and submit "preProcess.bash". If successful, the script automatically submits "precursorStage2Run.bash". After completion, a directory named "36000.essentials" should be generated. The stage-2_1 precursor modifies the bottom boundary condition "qwall" to trigger the development of CBL, NBL, or SBL conditions depending on the selected directory. Users should ensure that the directory structure is preserved because files are copied using relative paths.

For the stage-2_2 precursor, navigate to any directory under "caseSettings/precursor_stage_2_2" except "precursor_stage_2_2_Master", and submit "runBDsampling.bash". If successful, the script automatically submits "onlyRecon_s2.bash" and "submittingPrepTurb.bash". This stage samples the boundary fields used as inflow conditions for the successor simulations. Each completed case requires approximately 500 GB of storage.

Finally, for the successor simulations, navigate to any directory under "caseSettings/suscessors" except "suscessor_Master", and submit "preRunturb.bash". If successful, the script automatically submits "runturb.bash" followed by "onlyRecon_s2.bash". After completion, a directory named "41400" should be generated. Each completed successor case requires approximately 1 TB of storage, however, one may reduce the storage requirement to around 250 GB by muting probe sampling in "controlDict". Note that the case directories are named consistently with the case labels provided in Table 1 in the article.

