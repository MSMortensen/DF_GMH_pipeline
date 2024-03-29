#!/bin/bash

##############################################################################
# This file contains the information and commands necessary to  analyse Ion-
# Torrent 16S rRNA gene sequencing data with standard settings for DTU FOOD Gut
# Microbio group using the DF-2022.1 conda environment. 
##############################################################################
##############################################################################
###                         HOW TO USE THIS SCRIPT                         ###
##############################################################################
#
# This file contains the commands necessary to run the current 16S rRNA gene 
# sequencing pipeline at Research Group for Gut, Microbes and Health (DTU FOOD).
#   NOTE: The file README.md contains installation instructions needed.
#
# All variables and settings for this pipeline is defined in "settings.sh"
#   NOTE: I recommend renaming the file for each analysis
#   REMEMBER: Update step 1 to point to the correct file
#
##############################################################################
###                                ALWAYS                                  ###
##############################################################################
#
#   - Create a <PROJECT_FOLDER> (in the pipeline folder)
#   - Run step 0
#       If resuming analyses remember to run step 0 before any script
#   - Run step 1-4
#   - Pipeline output is exported into an .RData file
#       The name of the file is written in the last line from step 4
#
##############################################################################
###                       MULTIPLE SEQUENCING RUNS                         ###
##############################################################################       
#
#   After finishing the prior analyses for all runs
#   - Run Step 5
#       This combines all .rds files found in $out_dir
#   - Pipeline output is exported into an .RData file
#       The name of the file is written in the last line from step 5
#
##############################################################################
###                               PIPELINE                                 ###
##############################################################################  

### STEP 0: PREPARATION
# Load environment
conda activate DF-2022.1.1

# Move to project folder
cd ~/Github/DF_GMH_pipeline/<PROJECT_FOLDER>

### STEP 1: SETTINGS AND FILES
# Copy files to <PROJECT_FOLDER>
cp ../analysis/* .

# Create a copy of the settings file for each run and make sure to change relevant variables

### STEP 2: VERIFICATION
# Check necessary folders exist and create if necessary. Also verifies that files are available. 
# > If multiple runs copy command and edit name of the settings file 
bash ../scripts/Check_if_ready.sh -s settings_A.sh
bash ../scripts/Check_if_ready.sh -s settings_B.sh

### STEP 3A: DEMULTIPLEX AND TRIM READS (incl QC)
# > If multiple runs copy command and edit name of the settings file
# > If multiple projects on one chip run only once and use one of the settings files
bash ../scripts/DemultiplexAndTrim.sh -s settings_A.sh

## STEP 3B: COPY FILES TO PROJECT FOLDER
# > If multiple projects on one chip run only once
bash ../scripts/RenameAndMove.sh -s settings_A.sh
bash ../scripts/RenameAndMove.sh -s settings_B.sh

### STEP 4: RUN DADA2
# With standard settings
# > If multiple runs copy command and edit name of the settings file 
Rscript --vanilla ../scripts/RunDADA2.R -s settings_A.sh
Rscript --vanilla ../scripts/RunDADA2.R -s settings_B.sh

### STEP 5: MERGE RUNS
# Can also just be run with the settings from one of the runs being merged ("Rscript --vanilla Merged_Analysis.R -s settings_A.sh")
Rscript --vanilla ../scripts/Merged_Analysis.R -o output -p <PROJECT_NAME>  -r ../DB -n 0