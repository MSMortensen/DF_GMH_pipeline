#!/bin/bash

# Created:  2022-05-20 (masmo)
# Updated:  2022-06-02  masmo   Updated cutadapt settings to find mock
#                               Script now excludes reads not matching indeces
##############################################################################
# 
# 
# 
##############################################################################
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "This script demultiplexes and trims IonTorrent 16S rRNA gene sequencing data"
   echo "with standard settings for DTU FOOD Gut Microbio group using the DF-2022.1 "
   echo "conda environment. Standard use is to run with settings from file."
   echo
   echo "Syntax: DemultiplexAndTrim.sh -s "SETTINGS_FILE" [options] [-h]"
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

if [ $nthreads==0 ]; then CORES=$(($(nproc)-1)); else CORES=$nthreads; fi

##############################################################################
###                           DEMULTIPLEX SAMPLES                          ###
##############################################################################
echo "1) Demultiplex samples"
# Demultiplex using cutadapt
cutadapt -e 0.15 --no-indels -g file:$INDECES --cores=$CORES -o "$DEMUX/$RUN_NAME.{name}.fastq" $INPUT

echo "Demultiplexed files are stored in: $DEMUX"
# Initiate trimming report
echo "sample status  in_reads        in_bp   too_short       too_long        too_many_n      out_reads       w/adapters      qualtrim_bp     out_bp" > $OUT/$RUN_NAME.cutadapt_report.txt

##############################################################################
###                               TRIM READS                               ###
##############################################################################
echo "2) Trim reads"
# Trim adapters from demultiplexed reads
# At the same time we trim the length of the reads to be between 110 and 180 bp 
# Primers: "5' PBU_Fw CCTACGGGAGGCAGCAG" and "3' PBR_Rev CCAGCAGCCGCGGTAAT"

for FILE in $(ls $DEMUX/$RUN_NAME.*.fastq); do 
    NAME=$(echo $FILE | sed "s/$DEMUX/$TRIM/g")
    cutadapt -e 1 -g $FP -a $RP -n2 -m $MINLEN -M $MAXLEN --cores=4 --discard-untrimmed --report=minimal -o $NAME $FILE | tail -n1 | sed -e "s|^|$NAME\t|g" >> $OUT/$RUN_NAME.cutadapt_report.txt
done

echo "Report for trimming is stored in: $OUT/$RUN_NAME.cutadapt_report.txt"
rm $TRIM/$RUN_NAME.unknown.fastq

##############################################################################
###                               RUN FASTQC                               ###
##############################################################################
echo "3) Run FastQC"
# Run fastqc to analyse sequence quality
cat $TRIM/$RUN_NAME.*.fastq > $RUN_NAME.fastq
fastqc -o $OUT -t 4 $RUN_NAME.fastq
rm $RUN_NAME.fastq

echo "Fastqc report post trimming is stored in: $OUT/${RUN_NAME}_fastqc.html"
