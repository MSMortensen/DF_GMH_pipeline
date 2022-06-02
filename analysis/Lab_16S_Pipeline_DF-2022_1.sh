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
#       The name of the file is written in the last line from step 4


### STEP 0: PREPARATION
# Load environment
conda activate DF-2022.1

##### RUN WITH MOCK
### STEP 1: SETTINGS
export ANALYSIS_FILE=Analysis_settings.sh 

### STEP 2: VERIFICATION
# Check necessary folders and create if necessary. Also verifies that files are available.
bash ./Check_if_ready.sh

### STEP 3: DEMULTIPLEX AND TRIM READS (incl QC)
bash ./DemultiplexAndTrim.sh

### STEP 4: RUN DADA2
Rscript --vanilla RunDADA2.R

##### RUN WITHOUT MOCK
### STEP 1: SETTINGS
export ANALYSIS_FILE=Analysis_settings_nm.sh 

### STEP 2: VERIFICATION
# Check necessary folders and create if necessary. Also verifies that files are available.
bash ./Check_if_ready.sh

### STEP 3: DEMULTIPLEX AND TRIM READS (incl QC)
bash ./DemultiplexAndTrim.sh
rm trimmed/NoMock.unknown.fastq

### STEP 4: RUN DADA2
Rscript --vanilla RunDADA2.R

##### RUN WITHOUT pooling
### STEP 1: SETTINGS
export ANALYSIS_FILE=Analysis_settings_nopool.sh 

### STEP 2: VERIFICATION
# Check necessary folders and create if necessary. Also verifies that files are available.
bash ./Check_if_ready.sh

### STEP 3: DEMULTIPLEX AND TRIM READS (incl QC)
bash ./DemultiplexAndTrim.sh
rm trimmed/NoPool.unknown.fastq

### STEP 4: RUN DADA2
Rscript --vanilla RunDADA2.R

##### RUN WITHOUT MOCK and without pooling
### STEP 1: SETTINGS
export ANALYSIS_FILE=Analysis_settings_nopoolnomock.sh 

### STEP 2: VERIFICATION
# Check necessary folders and create if necessary. Also verifies that files are available.
bash ./Check_if_ready.sh

### STEP 3: DEMULTIPLEX AND TRIM READS (incl QC)
bash ./DemultiplexAndTrim.sh
rm trimmed/NopoolNomock.unknown.fastq

### STEP 4: RUN DADA2
Rscript --vanilla RunDADA2.R

### STEP 5: MERGE RUNS
Rscript --vanilla Merged_Analysis.R
