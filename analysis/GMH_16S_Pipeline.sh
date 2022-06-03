# Created:  2022-05-20 (masmo)
# Updated:   
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
#   NOTE: The file READ_ME.txt contains installation instructions needed.
#
# All variables and settings for this pipeline is defined in "Analysis_settings.sh"
#   NOTE: I recommend renaming the file for each analysis
#   REMEMBER: Update step 1 to point to the correct file
#
# # To run the analysis always
#   - Run step 0
#       If resuming analyses remember to run step 0 and step 1 before any script
# 
##############################################################################
###                         SINGLE SEQUENCING RUN                          ###
##############################################################################
# 
#   - Run step 1-4
#   - Import the produced .RData file into R and continue the analysis
#       The name of the file is written in the last line from step 4
#
##############################################################################
###                       MULTIPLE SEQUENCING RUNS                         ###
##############################################################################       
#
#   - Create a copy of Analysis_settings for each run
#   - Ensure that the following IF the content and name of the settings file.
#   - Run step 0
#       If resuming analyses remember to run step 0 and step 1 before any script
#   - Run step 1-4 for each run
#       REMEMBER TO UPDATE THE NAME OF THE FILE IN STEP 1
#   - Run Step 5
#       This combines all .rds files found in $out_dir
#   - Import the produced .RData file into R and continue the analysis
#       The name of the file is written in the last line from step 5

## As this is a combination of multiple analyses on the sample data I will only run STEP 3 once and then just copy the files
### STEP 0: PREPARATION
# Load environment
conda activate DF-2022.1

### STEP 1: SETTINGS
# Create a copy of the settings file and make sure to change relevant variables

### STEP 2: VERIFICATION
# Check necessary folders and create if necessary. Also verifies that files are available.
bash ../scripts/Check_if_ready.sh -s settings.sh

### STEP 3: DEMULTIPLEX AND TRIM READS (incl QC)
bash ../scripts/DemultiplexAndTrim.sh -s settings.sh

### STEP 4: RUN DADA2
# With standard settings
Rscript --vanilla ../scripts/RunDADA2.R -s settings.sh

### STEP 5: MERGE RUNS
# Can also just be run with the settings from one of the runs being merged ("Rscript --vanilla Merged_Analysis.R -s settings.sh")
Rscript --vanilla ../scripts/Merged_Analysis.R -o output -r ../DB -n 0
