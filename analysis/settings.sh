#!/bin/bash

# Created:  2022-05-20  masmo
# Updated:  2022-05-31  masmo   Added option for using sample names
#
##############################################################################
###                         HOW TO USE THIS SCRIPT                         ###
##############################################################################
#
# This file defines the settings used to run the current 16S rRNA gene.
# To use alternative settings create a copy of this file instead of changing here
#	REMEMBER:	Step 1 of the pipeline must be updated accordingly.
#				If the RUN_NAME is changed the updated analysis can
#				be run in the same folder
#
# This script has two types of settings:
# 	1) General settings - 	Will often have to be changed between analyses
#	2) Advanced settings - 	Are optimised for our setup and should not be
#							modified without carefull consideration
#
##############################################################################
###                            GENERAL SETTINGS                            ###
##############################################################################
#
### PIPELINE
# Choose if the full DADA2 script should be run or only the ASV calling.
#	Default is to run the full analysis. 
#	If set to "PARTIAL", the analysis will finish after ASV calling.
export ANALYSIS="PARTIAL"
#
### RUN NAME
# Define the run name. This name will be used for all output files.
#	USE ONLY ALPHANUMERIC CHARACTERS
export RUN_NAME="standard" 
#
### USE SAMPLE NAMES
# Define if generic sample names are used or run specific names. If not using
# generid names, then a tab-separated file with a column of the used tags and 
# sample names must be given
#   VALUES: false   (default - )
#           true    (Next variable MUST be set)
export USE_SAMPLE_NAMES=true
#
# Name of file with sample names and tags
#   NOTE:   This file does not need to contain all samples, only the ones that
#           matches samples to include.
#   NOTE:   This variable is only used if USE_SAMPLE_NAMES is true
export SAMPLE_FILE="sample_tags.tsv"
#
### FILES AND FOLDERS
# Name of the input fastq file
#	NOTE: Must be updated
export INPUT="R_2021_10_05_12_05_48_user_GSS5-0533-66-Martin_Frederik_Laursen_run_10.fastq"
#
# Path where demultiplexed files will be stored
#	NOTE: These files can be removed after analyses
export DEMUX="demux"
#
# Path where trimmed files will be stored
#	NOTE: These files can be removed after analyses
export TRIM="trimmed"
#
# Input path for DADA2
#	NOTE: This should normally be the folder with the trimmed files
export in_dir="trimmed"
#
# Path where CutAdapt output files will be stored
export OUT="output"
#
# Path where DADA2 output files will be stored
#	NOTE: 	This should normally match the output folder for CutAdapt
#			Any runs that has to be merged must have the same output folder
export out_dir="output"
#
# Path where DADA2 stored filtered fastq files
#	NOTE: 	These files can be removed after analyses
export filter_dir="./filtered"
#
# Path to the reference database.
#	NOTE:	The reference database should not change between runs
#			Therefore, this should always be the same for you
export reference_dir="~/DB" 
#
##############################################################################
###                           ADVANCED SETTINGS                            ###
##############################################################################
#
### PCR PRIMERS
# These should not be changed if sequencing 16S on IonTorrent
# Forward primer: "5' PBU_Fw CCTACGGGAGGCAGCAG"
export FP="CCTACGGGAGGCAGCAG"
# Reverse primer: "3' PBR_Rev CCAGCAGCCGCGGTAAT"
export RP="CCAGCAGCCGCGGTAAT"
#
### CUTADAPT TRIMMING
# These values are linked to the primerset used
# Minimum read length after removal of primers
export MINLEN=110
# Maximin read length after removal of primers
export MAXLEN=180
#
### DADA2 SETTINGS
# These settings are adapted to running DADA2 for IonTorrent sequencing
#
# Position at which to truncate reads.
#	NOTE:	Reads shorter than truncLen will be discarded. 
#			Special values: 0 - no truncation or length filtering
export truncLen=0
#
# Number of nucleotides to remove from the start of each read
# 	NOTE:	Should be less than truncLen for obvious reasons
export trimLeft=0 
#
# Reads with expected errors higher than maxEE are discarded
export maxEE=2
#
# Minimim quality score for trimming
#	NOTE:	Reads are truncated at the first instance of quality score truncQ
#			If the read is then shorter than truncLen, it is discarded
export truncQ=2
#
# Remove reads with length greater than maxLen. maxLen is enforced on the raw reads.
#	NOTE:	Default Inf - no maximum
export maxLen=Inf 
# 
# The method used to pool (or not) samples during denoising. 
#	Valid options: 
#		independent: No pooling, samples are denoised indpendently. 
#		pseudo: Samples are 'pseudo-pooled' for denoising. 
#		pool: (Default) Samples are pooled for denoising.
export poolMethod="pool"
#
# The method used to remove chimeras. 
#	Valid options are: 
#		none: No chimera removal is performed. 
#		pooled: (Default) All reads are pooled prior to chimera detection. 
#		consensus: Chimeras are detect in samples individually, and a consensus decision is made for each sequence variant
export chimeraMethod="pooled" 
#
# The minimum abundance of potential 'parents' of a sequence being tested as chimeric, 
# expressed as a fold-change versus the abundance of the sequence being tested. 
#	NOTE:	Values should be greater than or equal to 1 
#			(i.e. parents should be more abundant than the sequence being tested)
export minParentFold=1.0
#
# The number of processor threads to use.
#	NOTE:	Special values: 0 - detect available cores and use all except one
#			This value will also be used for CutAdapt
export nthreads=0 
#
# The minimum number of reads to learn the error model from. 
#	NOTE:	Special values: 0 - Use all input reads
export nreads_learn=1000000
#
# The cost of gaps in homopolymer regions (>=3 repeated bases).
#	NOTE:	Default is -1
export HOMOPOLYMER_GAP_PENALTY=-1 
#
# BAND_SIZE for Needleman-Wunsch alignment
#	NOTE:	Default value is 32
#	NOTE:	If set to a netagive number banding is turned off (runs full Needleman-Wunsch)
#	NOTE:	Banded Needleman-Wunsch alignments can be diabled by setting this variable to: NULL
export BAND_SIZE=32 
#
##############################################################################
###                          AUTOMATICALLY ADDED                           ###
##############################################################################
#
### INDEX FILE
# The setting for the index file will automatically be added here by STEP 2
