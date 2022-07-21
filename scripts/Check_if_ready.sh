#!/bin/bash

# Created:  2022-05-20  masmo
# Updated:  2022-05-31  masmo   Added option for using sample names
#           2022-06-02  masmo   Corrected printed name of final phyloseq file
#
#
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "This script verifies that the necessary files and folders are available for the pipeline,"
   echo "creates folders and files if needed, and prints the relevant folders and locations"
   echo
   echo "Syntax: Check_if_ready.sh -s "SETTINGS_FILE" [options] [-h]"
   echo "options:"
   echo "  -s     Name of settings file to check"
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
###                           CHECK SETTINGS FILE                          ###
##############################################################################
if [ ! -f $ANALYSIS_FILE ]
    then 
        echo "WARNING: The file $ANALYSIS_FILE is missing, without this file the pipeline cannot run"
        exit 0
    else 
        source $ANALYSIS_FILE
        echo "The analyses settings are read from $ANALYSIS_FILE and stored in $out_dir/${RUN_NAME}_settings.txt"
fi

##############################################################################
###                          ANALYSIS INFORMATION                          ###
##############################################################################
# General Info
echo "This analysis if for run $SEQ_NAME as part of the project $PROJECT_NAME. "
if [ $ANALYSIS == "PARTIAL" ]
    then
        echo "A partial analysis will be run. The script 'RunDADA2.R' will stop after asv calling and $out_dir/$PROJECT_NAME.$SEQ_NAME.seqtab_out.rds must be included as input to 'Merge_Runs.R'"
    else 
        echo "A full analysis will be run. The script 'RunDADA2.R' will produce a phyloseq object named 'phy' which can be loaded into R by loading the .RData file ($out_dir/$PROJECT_NAME.$SEQ_NAME.phyloseq_object.RData)"
fi

##############################################################################
###                        CONFIRM FOLDERS AND FILES                       ###
##############################################################################
# Create folders and check files are available
echo "WARNINGS:"
if [[ ! -e $TRIM ]]; then mkdir $TRIM; elif [[ ! -d $TRIM ]]; then echo "$TRIM already exists but is not a directory"; else echo "Any files in the folder $TRIM might be overwritten"; fi
if [[ ! -e $DEMUX ]]; then mkdir $DEMUX; elif [[ ! -d $DEMUX ]]; then echo "$DEMUX already exists but is not a directory"; else echo "Any files in the folder $DEMUX might be overwritten"; fi
if [[ ! -e $OUT ]]; then mkdir $OUT; elif [[ ! -d $OUT ]]; then echo "$OUT already exists but is not a directory"; else echo "Any files in the folder $OUT might be overwritten"; fi
if [ ! -f "$INPUT" ]; then echo "The input file ($INPUT) is missing."; fi
if [ "$TRIM" != "$in_dir" ]; then if [[ ! -e $in_dir ]]; then mkdir $in_dir; elif [[ ! -d $in_dir ]]; then echo "$in_dir already exists but is not a directory"; else echo "Any files in the folder $in_dir might be overwritten"; fi; fi

PLOC=$(echo $PWD | sed 's![^/]*$!!')
if [ ! -f "$PLOC/scripts/RunDADA2.R" ]; then echo "The Rscript (RunDADA2.R) is missing. please copy to $PLOC/scripts/\n"; fi
if [ $reference_dir == "../DB" ]; then export reference_dir=${PLOC}DB; fi
if [ ! -f "$reference_dir/rdp_train_set_18.fa.gz" ]; then echo "The database to assign taxonomy is missing. please copy to $reference_dir"; fi
if [ ! -f "$reference_dir/rdp_species_assignment_18.fa.gz" ]; then echo "The database to assign species is missing. please copy to $reference_dir"; fi
if [ ! -f "$reference_dir/$INDECES" ]; then echo "The index file to demultiplex the fastq file is missing. please copy to $reference_dir"; fi

if [ $SAMPLE_FILE == "NULL" ]; then echo "No sample file has been provided to rename demultiplexed files"; elif [ ! -f $SAMPLE_FILE ]; then echo "$SAMPLE_FILE is missing. This will cause renaming to fail and use the RunDADA2 script as if the \"SAMPLE FILE\" variable was not set"; fi