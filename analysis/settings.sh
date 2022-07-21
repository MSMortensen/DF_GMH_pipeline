#!/bin/bash

##############################################################################
###                         HOW TO USE THIS SCRIPT                         ###
##############################################################################
#
# This script defines the settings used to run the current 16S rRNA gene GMH  
# Ion Torrent 16S rRNA gene pipeline. 
#   General settings should be updated for each run
#   Advanced settings should only be updated after careful consideration.
# 
#	REMEMBER:	Fastq files from multiple runs can be demultiplexed in the
#               same folder, IF the sequencing run names differ
#
##############################################################################
###                          DEMULTIPLEX SETTINGS                          ###
##############################################################################
#
### SEQUENCING RUN NAME. This name will be used for all output files.
#	USE ONLY ALPHANUMERIC CHARACTERS
export SEQ_RUN="Run220707" 
#                               
### INPUT FILE
# Name of the input fastq file
#	NOTE:   Must be updated and it is recommend to use the full path
export INPUT="/mnt/c/KaliFiles/Julius/R_2022_07_07_13_02_07_user_GSS5-0533-96-Bodil_Madsen_Julius.fastq"
#
##############################################################################
###                            NAMING SETTINGS                             ###
##############################################################################
#
#                               
### PROJECT NAME
# Define the run name. This name will be used for all output files.
#	USE ONLY ALPHANUMERIC CHARACTERS
export PROJECT_NAME="Julius" 
#
### SAMPLE NAMES
# Name of file with sample names and tags
#   NOTE:   This file does not need to contain all samples, only the ones that
#           matches samples to include. b
#   NOTE:   If not used set to NULL
export SAMPLE_FILE="Julius_sample_tags.tsv"
#
##############################################################################
###                             DADA2 SETTINGS                             ###
##############################################################################
#
### PIPELINE
# Choose if the full DADA2 script should be run or only the ASV calling.
#	NOTE:   Default is to run the full analysis. 
#	        If set to "PARTIAL", the analysis will finish after ASV calling.
export ANALYSIS="FULL"
#
##############################################################################
###                     ADVANCED DEMULTIPLEX SETTINGS                      ###
##############################################################################
#
### FOLDERS
# Folder where demultiplexed files will be stored
#	NOTE:   These files can be removed after analyses
export DEMUX="demux"
#
# Folder where trimmed files will be stored
#	NOTE:   These files can be removed after analyses
#           IF updated remember to update input folder for DADA2
export TRIM="${SEQ_RUN}_fastq"
#
# Folder where CutAdapt output files will be stored
export OUT="output"
#
### INDEX FILE
# The setting for the index file will automatically be added here by STEP 2
export INDECES="../DB/Indeces.fasta"
#
### PCR PRIMERS
# These should not be changed if sequencing 16S on IonTorrent
#   Forward primer: "5' PBU_Fw CCTACGGGAGGCAGCAG"
export FP="CCTACGGGAGGCAGCAG"
#   Reverse primer: "3' PBR_Rev CCAGCAGCCGCGGTAAT"
export RP="CCAGCAGCCGCGGTAAT"
#
### CUTADAPT TRIMMING
# These values are linked to the primerset used
# Minimum read length after removal of primers
export MINLEN=110
# Maximin read length after removal of primers
export MAXLEN=180
#
# The number of processor threads to use.
#	NOTE:	Special values: 0 - detect available cores and use all except one
export nthreads=0 
#
##############################################################################
###                        ADVANCED DADA2 SETTINGS                         ###
##############################################################################
#
### FOLDERS
# Input path for DADA2
#	NOTE:   If a sample file is provided input fastq files will be moved
if [ -z $SAMPLE_FILE ]; then export in_dir=$TRIM; elif [ ! -f $SAMPLE_FILE ]; then export in_dir=$TRIM; else export in_dir="${PROJECT_NAME}_fastq"; fi
#
# Path where DADA2 output files will be stored
#	NOTE: 	This should normally match the output folder for CutAdapt
#			Any runs that has to be merged must have the same output folder
export out_dir=$OUT
#
# Path where DADA2 stored filtered fastq files
#	NOTE: 	These files can be removed after analyses
export filter_dir="./filtered"
#
# Path to the reference database.
#	NOTE:	The reference database should not change between runs
#			Therefore, this should always be the same for you
export reference_dir="../DB" 
#
# Position at which to truncate reads.
#       NOTE:   Reads shorter than truncLen will be discarded.
#                       Special values: 0 - no truncation or length filtering
export truncLen=0
#
# Number of nucleotides to remove from the start of each read
#       NOTE:   Should be less than truncLen for obvious reasons
export trimLeft=0
#
# Reads with expected errors higher than maxEE are discarded
export maxEE=2
#
# Minimim quality score for trimming
#       NOTE:   Reads are truncated at the first instance of quality score truncQ
#                       If the read is then shorter than truncLen, it is discarded
export truncQ=2
#
# Remove reads with length greater than maxLen. maxLen is enforced on the raw reads.
#       NOTE:   Default Inf - no maximum
export maxLen=Inf
#
# The method used to pool (or not) samples during denoising.
#       Valid options:
#               independent: No pooling, samples are denoised indpendently.
#               pseudo: Samples are 'pseudo-pooled' for denoising.
#               pool: (Default) Samples are pooled for denoising.
export poolMethod="pool"
#
# The method used to remove chimeras.
#       Valid options are:
#               none: No chimera removal is performed.
#               pooled: (Default) All reads are pooled prior to chimera detection.
#               consensus: Chimeras are detect in samples individually, and a consensus decision is made for each sequence variant
export chimeraMethod="pooled"
#
# The minimum abundance of potential 'parents' of a sequence being tested as chimeric,
# expressed as a fold-change versus the abundance of the sequence being tested.
#       NOTE:   Values should be greater than or equal to 1
#                       (i.e. parents should be more abundant than the sequence being tested)
export minParentFold=1.0
#
# The minimum number of reads to learn the error model from.
#       NOTE:   Special values: 0 - Use all input reads
export nreads_learn=1000000
#
# The cost of gaps in homopolymer regions (>=3 repeated bases).
#       NOTE:   Default is -1
export HOMOPOLYMER_GAP_PENALTY=-1
#
# BAND_SIZE for Needleman-Wunsch alignment
#       NOTE:   Default value is 32
#       NOTE:   If set to a netagive number banding is turned off (runs full Needleman-Wunsch)
#       NOTE:   Banded Needleman-Wunsch alignments can be diabled by setting this variable to: NULL
export BAND_SIZE=32
