#! ~/miniconda3/envs/DF-2022.1/bin/Rscript

# Created:  2022-05-20  masmo
# Updated:  2022-05-31  masmo Fixed problems loading settings
##############################################################################
# This R script runs with standard settings for DTU FOOD GutMicrobio group 
# using the DF-2022.1 conda environment. Standard use is to run with settings 
# from file, but the script can be run with settings from commandline
##############################################################################
##############################################################################
###                        DESCRIPTION OF ARGUMENTS                        ###
##############################################################################
# NOTE: All numeric arguments should be zero or positive.
# NOTE: All numeric arguments save maxEE are expected to be integers.
# NOTE: ALL ARGUMENTS DEVIATING FROM STANDARD MUST BE NAMED!
# NOTE: IF $ANALYSIS_FILE is found in the working directory, 
#		all other arguments will be ignored

#Prep
cat(R.version$version.string, "\n")
errQuit <- function(mesg, status=1) { message("Error: ", mesg); q(status=status) }
args <- commandArgs(TRUE)

##############################################################################
###                     	   LOAD LIBRARIES                              ###
###############################################################################
cat("0) Load libraries and variables\n")
suppressMessages(suppressWarnings(library(methods)))
suppressMessages(suppressWarnings(library(dada2)))
suppressMessages(suppressWarnings(library(RcppParallel)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(phyloseq)))
suppressMessages(suppressWarnings(library(DECIPHER)))
suppressMessages(suppressWarnings(library(phangorn)))

cat("PACKAGES USED TO CALL ASVS:\nDADA2:", as.character(packageVersion("dada2")), "|",
    "Rcpp:", as.character(packageVersion("Rcpp")), "|",
    "RcppParallel:", as.character(packageVersion("RcppParallel")), "|",
    "parallel:", as.character(packageVersion("parallel")), "|",
    "ggplot2:", as.character(packageVersion("ggplot2")), "|",
    "optparse:", as.character(packageVersion("optparse")), "\n")
cat("ADDITIONAL PACKAGES USED IF FULL PIPELINE IS RUN:\nphyloseq", as.character(packageVersion("phyloseq")), "|",
    "DECIPHER:", as.character(packageVersion("DECIPHER")), "|",
    "phangorn:", as.character(packageVersion("phangorn")), "\n")

##############################################################################
###                     CREATE DEFAULT ARGUMENT LIST                       ###
##############################################################################
# import and format arguments
option_list = list(
  make_option(c("-s", "--settings_file"), type="character", default=NULL, help='If a settings file is provided here the system variable (ANALYSIS_FILE) will be ignored', metavar="number"),
  make_option(c("-A", "--ANALYSIS"), type="character", default="FULL", help='When set to "PARTIAL", the analysis will be stopped after ASV calling as all following analyses are performed by the script Merge_Runs.R', metavar="number"),
  make_option(c("-R", "--RUN_NAME"), type="character", default="./trimmed", help="File path to directory with the .fastq files to be processed", metavar="character"),
  make_option(c("-i", "--in_dir"), type="character", default="./trimmed", help="File path to directory with the .fastq files to be processed", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default="./out", help="File path for all output", metavar="character"),
  make_option(c("-f", "--filter_dir"), type="character", default="./filtered", help="File path to directory in which to write the filtered .fastq.gz files. These files are intermediate for the full workflow. Currently they remain after the script finishes.", metavar="character"),
  make_option(c("-r", "--reference_dir"), type="character", default="../DB", help="This must point to the folder with the reference database",metavar="character"),
  make_option(c("-t", "--truncLen"), type="numeric", default=0, help="The position at which to truncate reads. Reads shorter than truncLen will be discarded. Special values: 0 - no truncation or length filtering.", metavar="number"),
  make_option(c("-T", "--trimLeft"), type="numeric", default=0, help="The number of nucleotides to remove from the start of each read. Should be less than truncLen for obvious reasons.", metavar="number"),
  make_option(c("-e", "--maxEE"), type="numeric", default=2, help="Reads with expected errors higher than maxEE are discarded.", metavar="number"),
  make_option(c("-Q", "--truncQ"), type="numeric", default=2, help="Reads are truncated at the first instance of quality score truncQ If the read is then shorter than truncLen, it is discarded.", metavar="number"),
  make_option(c("-L", "--maxLen"), type="numeric", default=Inf, help="Remove reads with length greater than maxLen. maxLen is enforced on the raw reads. Default Inf - no maximum.", metavar="number"),
  make_option(c("-p", "--poolMethod"), type="character", default="pool", help="The method used to pool (or not) samples during denoising. Valid options are: independent: No pooling, samples are denoised indpendently. pseudo: Samples are 'pseudo-pooled' for denoising. pool: (Default) Samples are pooled for denoising", metavar="character"),
  make_option(c("-c", "--chimeraMethod"), type="character", default="pooled", help="The method used to remove chimeras. Valid options are: none: No chimera removal is performed. pooled: All reads are pooled prior to chimera detection. consensus: Chimeras are detect in samples individually, and a consensus decision is made for each sequence variant", metavar="character"),
  make_option(c("-F", "--minParentFold"), type="numeric", default=1.0, help="The minimum abundance of potential 'parents' of a sequence being tested as chimeric, expressed as a fold-change versus the abundance of the sequence being tested. Values should be greater than or equal to 1 (i.e. parents should be more abundant than the sequence being tested).", metavar="number"),
  make_option(c("-n", "--nthreads"), type="numeric", default=0, help="The number of threads to use. Special values: 0 - detect available cores and use all except one.", metavar="number"),
  make_option(c("-N", "--nreads_learn"), type="numeric", default=1000000, help="The minimum number of reads to learn the error model from. Special values: 0 - Use all input reads.", metavar="number"),
  make_option(c("-H", "--HOMOPOLYMER_GAP_PENALTY"), type="numeric", default=-1, help="The cost of gaps in homopolymer regions (>=3 repeated bases). Default is -1", metavar="number"),
  make_option(c("-B", "--BAND_SIZE"), type="numeric", default=32, help="When set, banded Needleman-Wunsch alignments are performed. The default value of BAND_SIZE is 32. Setting BAND_SIZE to a negative number turns off banding (i.e. full Needleman-Wunsch).", metavar="number")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

##############################################################################
###                        IMPORT ARGUMENTS FROM FILE                      ###
##############################################################################
# Import settings from settings file if provided
if (!is.null(opt$settings_file) | file.exists(Sys.getenv("ANALYSIS_FILE"))){
  if (is.null(opt$settings_file)){
    tmp <- read.table(Sys.getenv("ANALYSIS_FILE"), header = F, comment.char = "#")
  } else tmp <- read.table(opt$settings_file, header = F, comment.char = "#")
  vars <- as.data.frame( t( stringr::str_split(tmp[,2], "=", simplify = T ) ) )
  colnames(vars) <- vars[1,]
  vars <- vars[2,]

  suppressWarnings(for (i in 1:ncol(vars)) {
    if (!is.na(as.numeric(vars[,i]))) vars[,i] <-as.numeric(vars[,i]) 
  })

  # replace settings in opt
  for (i in 1:(length(opt)-1)){
    if (names(opt)[i] != "settings_file") opt[i][[1]] <- vars[names(opt)[i]][[1]]
  }

  rm(tmp, vars)
}
for (i in 1:(length(opt)-1)){
  if (is.character(opt[i][[1]])) opt[i][[1]] <- gsub("\\\"","",opt[i][[1]])
}

##############################################################################
###                            VALIDATE ARGUMENTS                          ###
##############################################################################
# Convert nthreads to the logical/numeric expected by dada2
if(opt$nthreads < 0) {
  errQuit("nthreads must be non-negative.")
} else if(opt$nthreads == 0) {
  multithread <- detectCores()-1 # detect and use all
} else if(opt$nthreads == 1) {
  multithread <- FALSE
} else {
  multithread <- nthreads
}

# Convert poolMethod
if(opt$poolMethod == "pool") opt$poolMethod <- TRUE else if (opt$poolMethod == "pseudo") opt$poolMethod <- "pseudo" else opt$poolMethod <- FALSE

# Output files are to be filenames (not directories) and are to be
# removed and replaced if already present.
if(!dir.exists(opt$out_dir)) {
      dir.create(opt$out_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}
cat("Output folder:", opt$out_dir, "\n")

tmp <- list.files(opt$out_dir, pattern=paste0(opt$RUN_NAME), full.names=FALSE)
tmp <- tmp[!grepl("html|zip|txt",tmp)]
if(length(tmp) > 0) errQuit("The planned output files already exists")

##############################################################################
###                           TRIM AND FILTER                              ###
##############################################################################
# Input directory is expected to contain .fastq file(s) that have not yet been 
# filtered and globally trimmed to the same length.
if(!dir.exists(opt$in_dir)) {
  errQuit("Input directory does not exist.")
} else {
  unfilts <- list.files(opt$in_dir, pattern=".fastq$", full.names=TRUE) # Finds fastq files
  unfilts <- unfilts[grepl(opt$RUN_NAME, unfilts)] # keeps only the ones named with the defined RUN_NAME
  if(length(unfilts) == 0) {
    errQuit("No input files with the expected filename format found.")
  }
}

# Trim and filter
cat("1) Filtering\n")
filts <- file.path(opt$filter_dir, basename(unfilts))
out <- suppressWarnings(filterAndTrim(unfilts, filts, truncLen=opt$truncLen, trimLeft=opt$trimLeft,
                                      maxEE=opt$maxEE, truncQ=opt$truncQ, rm.phix=TRUE,
                                      multithread=multithread, maxLen=opt$maxLen))
cat(ifelse(file.exists(filts), ".", "x"), sep="")
cat("\n")

if(sum(file.exists(filts)) == 0) { # All reads were filtered out
  errQuit("No reads passed the filter (was truncLen longer than the read length?)", status=2)
}

##############################################################################
###                          LEARN ERROR RATES                             ###
##############################################################################
# Dereplicate enough samples to get nreads.learn total reads
cat("2) Learning Error Rates\n")
err <- suppressWarnings(learnErrors(filts, nreads=opt$nreads_learn, multithread=multithread,
                   HOMOPOLYMER_GAP_PENALTY=opt$HOMOPOLYMER_GAP_PENALTY, BAND_SIZE=opt$BAND_SIZE))

# Just as a sanity check
ggsave(filename = paste(opt$RUN_NAME,"dada2_ErrorModel.pdf", sep = "."),
  plot = plotErrors(err, nominalQ=TRUE),
  device = "pdf",
  path = opt$out_dir,
  width = 18,
  height = 15,
  units = "cm",
  dpi = 300)

cat("Dada2 error model:", file.path(opt$out_dir, paste(opt$RUN_NAME,"dada2_ErrorModel.pdf", sep = ".")),"\n")

##############################################################################
###                         PERFORM ASV CALLING                            ###
##############################################################################
# Run dada2 with learned error rates
cat("3) Denoise samples \n")
drp <- vector("list", length(filts))
for(j in seq(length(filts))) {
  drp[[j]] <- derepFastq(filts[[j]])
  if (j < length(filts)) cat(".") else cat(".\n")
}

dadaFs <- dada(drp, err=err, multithread=multithread, pool=opt$poolMethod, HOMOPOLYMER_GAP_PENALTY=opt$HOMOPOLYMER_GAP_PENALTY, BAND_SIZE=opt$BAND_SIZE)

### Make sequence table
seqtab <- makeSequenceTable(dadaFs)
cat("Dataset contains ", ncol(seqtab), " ASVs from ",nrow(seqtab), " samples\n", sep = "")

# Inspect distribution of sequence lengths
dat <- data.frame(number = 1:ncol(seqtab), seq_length = nchar(getSequences(seqtab)))
ggsave(filename = paste(opt$RUN_NAME,"ASV_sequence_lengths.pdf",sep = "."),
  plot = ggplot(dat, aes(x = seq_length)) + geom_histogram(binwidth = 1),
  device = "pdf",
  path = opt$out_dir,
  width = 18,
  height = 15,
  units = "cm",
  dpi = 300)

cat("Histogram of ASV lengths:", file.path(opt$out_dir, paste(opt$RUN_NAME,"ASV_sequence_lengths.pdf", sep = ".")),"\n")

# Remove chimeras
cat("4) Remove chimeras (method = ", opt$chimeraMethod, ")\n", sep="")
if(opt$chimeraMethod %in% c("pooled", "consensus")) {
  seqtab.nochim <- removeBimeraDenovo(seqtab, method=opt$chimeraMethod, minFoldParentOverAbundance=opt$minParentFold, multithread=multithread)
} else { # No chimera removal, copy seqtab to seqtab.nochim
  seqtab.nochim <- seqtab
}

### Remove ASVs with sequences shorter than 75 bp
keep <- nchar(colnames(seqtab.nochim)) >= 75
cat(sum(!keep), " ASVs were removed because they were shorter than 75 bp.\n", sum(seqtab.nochim[,!keep]), " reads (",sum(seqtab.nochim[,!keep])/sum(seqtab.nochim),"%) were removed\n")
seqtab.nochim <- seqtab.nochim[,keep]

##############################################################################
###               REPORT READ FRACTIONS THROUGH PIPELINE                   ###
##############################################################################
# Create report
cat("5) Report read numbers through the pipeline\n")
# Handle edge cases: Samples lost in filtering; One sample
track <- cbind(out, matrix(0, nrow=nrow(out), ncol=2))
colnames(track) <- c("input", "filtered", "denoised", "non-chimeric")
passed.filtering <- track[,"filtered"] > 0
track[passed.filtering,"denoised"] <- rowSums(seqtab)
track[passed.filtering,"non-chimeric"] <- rowSums(seqtab.nochim)
write.table(track, file.path(opt$out_dir,paste(opt$RUN_NAME,"sample_reads.tsv", sep = ".")), sep="\t", row.names=TRUE, col.names=NA,  quote=FALSE)
cat("Report stored in: ",file.path(opt$out_dir,paste(opt$RUN_NAME,"sample_reads.tsv", sep = ".")),"\n",sep = "")

##############################################################################
###                 EXPORT AND STOP IF PARTIAL ANALYSIS                    ###
##############################################################################
# Save seqtab
sample.names <- stringr::str_remove(basename(filts),"\\.fastq")
seq.out <- seqtab.nochim
row.names(seq.out) <- paste(opt$RUN_NAME, sample.names, sep = "_")
write.table(seq.out,file.path(opt$out_dir,paste(opt$RUN_NAME,"seqtab_out.tsv",sep = ".")), sep = "\t",quote = FALSE)
saveRDS(seq.out, file.path(opt$out_dir,paste(opt$RUN_NAME,"seqtab_out.rds",sep = ".")))
cat("seqtab file (needed to merge runs) stored in: ",file.path(opt$out_dir,paste(opt$RUN_NAME,"seqtab_out.rds", sep = ".")),"\n",sep = "")

# Stop the script if running partial analysis
if (opt$ANALYSIS == "PARTIAL") message("ASV calling complete, remaining analyses will be performed by the script Merged_Analysis.R"); q(status=0) 

##############################################################################
###                     ASSIGN TAXONOMY AND PHYLOGENY                      ###
##############################################################################
# Assign taxonomy
cat("6) Assign taxonomy\n")
taxa <- assignTaxonomy(seqtab.nochim, file.path(opt$reference_dir, "rdp_train_set_18.fa.gz"), multithread=TRUE, verbose=TRUE)
taxa.plus <- addSpecies(taxa, file.path(opt$reference_dir, "rdp_species_assignment_18.fa.gz"), verbose=TRUE,allowMultiple = T)
asv.names <- paste("ASV", stringr::str_pad(1:ncol(seqtab.nochim),4, pad = "0"), sep = "_")
row.names(taxa.plus) <- asv.names
write.table(taxa.plus, file = file.path(opt$out_dir,paste(opt$RUN_NAME,"taxonomy.txt", sep = ".")), row.names = asv.names)
cat("Taxonomy stored in: ",file.path(opt$out_dir,paste(opt$RUN_NAME,"taxonomy.txt", sep = ".")),"\n",sep = "")

### Construct Phylogenetic Tree
cat("7) Create phylogenetic tree\n")
# Extract sequences from DADA2 output
sequences<-getSequences(seqtab.nochim)
names(sequences)<-paste("ASV", stringr::str_pad(1:ncol(seqtab.nochim),4, pad = "0"), sep = "_")

# Run Sequence Alignment (MSA) using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA, verbose = FALSE)

# Change sequence alignment output into a phyDat structure
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

# Create distance matrix
dm <- dist.ml(phang.align)

# Perform Neighbor joining
treeNJ <- NJ(dm) # Note, tip order != sequence order

# Internal maximum likelihood
fit = pml(treeNJ, data=phang.align)

# negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

##############################################################################
###                         CREATE PHYLOSEQ OBJECT                         ###
##############################################################################
# Create phyloseq object
cat("8) Create phyloseq object\n")
# Rename ASVs and samples
asv_table <- seqtab.nochim
colnames(asv_table) <- paste("ASV", stringr::str_pad(1:ncol(seqtab.nochim),4, pad = "0"), sep = "_")
row.names(asv_table) <- paste(opt$RUN_NAME, sample.names, sep = "_")

# Create metadata file
map <- data.frame(
  row.names = row.names(asv_table),
  Tag = stringr::str_split(sample.names,"\\.",simplify=TRUE)[,2],
  Run = opt$RUN_NAME,
  reads = rowSums(asv_table),
  ASVs = rowSums(asv_table > 0)
  )

phy <- phyloseq(otu_table(asv_table, taxa_are_rows=FALSE), tax_table(taxa.plus),phy_tree(fitGTR$tree), refseq(DNAStringSet(sequences)), sample_data(map))

# Currently, the phylogenetic tree is not rooted Though it is not necessary here, you will need to root the tree if you want to calculate any phylogeny based diversity metrics (like Unifrac)
set.seed(711) # As there is some randomness in this process, setting the seed ensures reproducibility
phy_tree(phy) <- root(phy_tree(phy), sample(taxa_names(phy), 1), resolve.root = TRUE)
is.rooted(phy_tree(phy))

##############################################################################
###                              CREATE PLOTS                              ###
##############################################################################
# Create a plot of richness compared to read_count
ggsave(filename = paste(opt$RUN_NAME,"ASVs_vs_reads.pdf", sep = "."),
  plot = ggplot(map, aes(x = reads, y = ASVs)) + geom_point() + geom_smooth(method=lm, se=FALSE, col='red', size=2),
  device = "pdf",
  path = opt$out_dir,
  width = 18,
  height = 15,
  units = "cm",
  dpi = 300)

cat("Plot comparing read depth and number of ASVs stored in: ",file.path(opt$out_dir,paste(opt$RUN_NAME,"ASVs_vs_reads.pdf", sep = ".")),"\n",sep = "")

##############################################################################
###                        EXPORT PHYLOSEQ OBJECT                          ###
##############################################################################
# Export data
save(phy, file = file.path(opt$out_dir,paste(opt$RUN_NAME,"phyloseq_object.RData", sep = ".")))
cat("Final phyloseq object (phy) stored in: ",file.path(opt$out_dir,paste(opt$RUN_NAME,"phyloseq_object.RData", sep = ".")),"\n",sep = "")

q(save="no",status=0)
