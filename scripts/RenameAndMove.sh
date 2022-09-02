#!/bin/bash

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "This script copies the demultiplexed files to a project"
   echo "folder and renames them based on the current project"
   echo
   echo "Syntax: RenameAndMove.sh -s "SETTINGS_FILE" [options] [-h]"
   echo "options:"
   echo "  -s     Name of settings file to use"
   echo "  -h     Print this Help."
   echo
}

############################################################
############################################################
# Main program                                             #
############################################################
############################################################
############################################################
# Process the input options. Add options as needed.        #
############################################################
# Get the options
while getopts ":hs:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      s) # Enter a name
         ANALYSIS_FILE=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

if [ -z "$ANALYSIS_FILE" ]; then
    echo "A settings file must be provided"
    exit 2;
fi
##############################################################################
###                         LOAD SETTINGS FROM FILE                        ###
##############################################################################
### load settings
source $ANALYSIS_FILE

##############################################################################
###                           DEMULTIPLEX SAMPLES                          ###
##############################################################################
echo "1) Copy files"
# Demultiplex using cutadapt

while read TAG NAME; do 
    cp $TRIM/$SEQ_RUN.$TAG.fastq $in_dir/$PROJECT_NAME.$SEQ_RUN.$NAME.fastq
done < $SAMPLE_FILE

##############################################################################
###                             SAVE SETTINGS                              ###
##############################################################################
echo "2) Save settings"
# Create settings file
cp $TRIM/$SEQ_RUN.settings $in_dir/$PROJECT_NAME.$SEQ_RUN.settings

printf "##############################################################################\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
printf "###                           RENAMING SETTINGS                            ###\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
printf "##############################################################################\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
printf "\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
printf "\tINPUT:\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
printf "\t\tRun name: $SEQ_RUN\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
printf "\t\tProject name: $PROJECT_NAME\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
printf "\t\tSample file: $SAMPLE_FILE\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
printf "\tOUTPUT:\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
printf "\t\tFastq files: $in_dir\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
printf "\n" >> $in_dir/$PROJECT_NAME.$SEQ_RUN.settings
