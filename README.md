HIRSCHHEYDT_STOFER_KERY


[![DOI](https://zenodo.org/badge/493551256.svg)](https://zenodo.org/badge/latestdoi/493551256)


This repository contains all the code used for the data simulation & analysis, and the code necessary to reproduce the figures and numbers presented in the article von Hirschheydt, G., Stofer, S., Kéry, M. “Mixed” occupancy designs: when do additional single-visit data improve the inferences from standard multi-visit models?


#  REPRODUCIBILITY OF RESULTS  ##
#-------------------------------#

All our results should be reproducible with the provided code.
Below, we describe the structure and content of this repository.
The file "workflow.Rmd" guides you through the R-files in the correct order to:
1. simulate the data with the same seed we used
2. fit occupancy models to the simulated data
3. exclude non-valid estimates (& counting the proportion of valid estimates)
4. fit linear mixed models to the standard errors of the estimates
5. reproduce the figures we created for the manuscript
6. calculate the bias of estimates and reproduce figures we presented in the supplementary materials.
IMPORTANT: The simulation of the data and the fitting of the occupancy models takes a long time: 6 R-files which each take between 6-10 days to run! Fitting linear mixed models to the estimates takes another 2-3 hours.
If you do not want to take the time to reproduce the simulation and occupancy model results, but you would nevertheless like to explore the output of ours, please contact Gesa von Hirschheydt so the data package (ca. 10 GB) can be sent to you.


#  STRUCTURE OF THIS REPOSITORY  ##
#---------------------------------#

# FILE: Hirschheydt_Stofer_Kery.Rproj
This file stores the information about this repository which was set up as an R-project in RStudio. You can open the project by clicking on this file. Opening the project in this way will automatically define this repository as the working directory of the R-session.
If you are not using RStudio or you do not want to work with R-projects, you must make sure that you set this repository as the working directory before you run any of the R-code.

# FILE: workflow.Rmd
This file explains in which order the provided R-files must be executed in order to reproduce the results we presented in the manuscript.

# FOLDER: 0_Singularity_container
This folder contains the ... and the ... file which we used to run the R-files in the folder '1_Simulate_data_and_fit_occupancy_models' on a High Performance Cluster.

# FOLDER: 1_Simulate_data_and_fit_occupancy_models
This folder contains 6 R-files. Each file was used to simulate data and fit occupancy models to the simulated data under one of 6 scenarios.
Running the R-files within this folder will create one additional .RData file within the folder. The new file contains the estimates from the occupancy models fitted to the simulated data. The simulated data itself is not stored.
Note: Each of the R-files in this folder takes between 6 and 10 days to run!

# FOLDER: 2_Fit_linear_models_to_SEs
This folder contains 1 R-file which was used to remove all non-valid estimates and bring all valid estimates into a single .RData file, count the proportion of valid estimates, fit linear mixed regression models to the standard errors of the parameter estimates and store the slope of these regressions, and estimate bias of parameter estimates.
Running the R-file within this folder will create one very large .RData file within the folder. The new file contains all the cleaned estimates, proportions of valid estimates, and the calculated bias.
Note: The R-file in this folder takes approximately 2 hours to run.

# FOLDER: 3_Make_figures_and_extract_numbers
This folder contains 2 R-files which were used to create the figures published in the manuscript and the supplementary materials and to extract the numbers presented in the manuscript.
Running the R-files within this folder will create 9 .tif files within the folder.
Note: The R-file in this folder takes max. 10 minutes to run.

