#! ~/miniconda3/envs/DF-2022.1/bin/Rscript

# Created:  2022-05-20 (masmo)
# Updated: 
##############################################################################
# This R script merges multiple runs with standard settings for DTU FOOD Gut
# Microbio group using the DF-2022.1 conda environment. Standard use is to run 
# with settings from file, but the script can be run with settings from commandline
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
###                     CREATE DEFAULT ARGUMENT LIST                       ###
##############################################################################
#load libraries
library("optparse",quietly = TRUE)

# import and format arguments
option_list = list(
  make_option(c("-o", "--out_dir"), type="character", default="./out", help="File path for all output", metavar="character"),
  make_option(c("-r", "--reference_dir"), type="character", default="~/DB", help="This file path has data from the prior script and this together",metavar="character"),
  make_option(c("-n", "--nthreads"), type="numeric", default=0, help="The number of threads to use. Special values: 0 - detect available cores and use all except one.", metavar="number")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

##############################################################################
###                     	   LOAD LIBRARIES                              ###
##############################################################################
suppressMessages(suppressWarnings(library(methods)))
suppressMessages(suppressWarnings(library(dada2)))
suppressMessages(suppressWarnings(library(RcppParallel)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(ggplot2)))

cat("PACKAGES USED TO CALL ASVS:\nDADA2:", as.character(packageVersion("dada2")), "|",
    "Rcpp:", as.character(packageVersion("Rcpp")), "|",
    "RcppParallel:", as.character(packageVersion("RcppParallel")), "|",
    "parallel:", as.character(packageVersion("parallel")), "|",
    "ggplot2:", as.character(packageVersion("ggplot2")), "\n")

suppressMessages(suppressWarnings(library(DECIPHER)))
suppressMessages(suppressWarnings(library(phangorn)))
suppressMessages(suppressWarnings(library(phyloseq)))
cat("ADDITIONAL PACKAGES USED TO CREATE PHYLOSEQ OBJECT:\nphyloseq", as.character(packageVersion("phyloseq")), "|",
  "DECIPHER:", as.character(packageVersion("DECIPHER")), "|",
  "phangorn:", as.character(packageVersion("phangorn")), "\n")
  
##############################################################################
###                        IMPORT ARGUMENTS FROM FILE                      ###
##############################################################################
# Import settings from "Analysis_settings.sh" if it is present
if (file.exists(Sys.getenv("ANALYSIS_FILE"))){
  
  tmp <- read.table(Sys.getenv("ANALYSIS_FILE"), header = F, comment.char = "#")
  vars <- as.data.frame( t( stringr::str_split(tmp[,2], "=", simplify = T ) ) )
  colnames(vars) <- vars[1,]
  vars <- vars[2,]

  suppressWarnings(for (i in 1:ncol(vars)) {
    if (!is.na(as.numeric(vars[,i]))) vars[,i] <-as.numeric(vars[,i]) 
  })

  vars.list <- as.list(vars)

  # replace settings in opt
  for (i in 1:(length(opt)-1)){
    opt[i][[1]] <- vars.list[names(opt)[i]][[1]]
  }

  rm(tmp, vars,vars.list)
}

for (i in 1:(length(opt)-1)){
  if (is.character(opt[i][[1]])) opt[i][[1]] <- gsub("\\\"","",opt[i][[1]])
}

 # replace reference dir setting
if (opt$reference_dir == "~/DB") {
    opt$reference_dir <- gsub("~",Sys.getenv("HOME"),opt$reference_dir)
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

##############################################################################
###                           IMPORT AND MERGE                             ###
##############################################################################
# All .rds files in the out_dir will be used
cat("1) merging runs\n")
incl.runs <- list.files(opt$out_dir, pattern="seqtab_out.rds$", full.names=TRUE)

cat("The following files are being imported: ",incl.runs,sep = "\n")

full.seqtab <- mergeSequenceTables(tables = incl.runs)

##############################################################################
###                     ASSIGN TAXONOMY AND PHYLOGENY                      ###
##############################################################################
# Assign taxonomy
cat("2) Assign taxonomy\n")
taxa <- assignTaxonomy(full.seqtab, file.path(opt$reference_dir, "rdp_train_set_18.fa.gz"), multithread=TRUE, verbose=TRUE)
taxa.plus <- addSpecies(taxa, file.path(opt$reference_dir, "rdp_species_assignment_18.fa.gz"), verbose=TRUE,allowMultiple = T)
asv.names <- paste("ASV", stringr::str_pad(1:ncol(full.seqtab),4, pad = "0"), sep = "_")
row.names(taxa.plus) <- asv.names
write.table(taxa.plus, file = file.path(opt$out_dir,"full_taxonomy.txt"), row.names = asv.names)
cat("Taxonomy stored in: ",file.path(opt$out_dir,"full_taxonomy.txt"),"\n",sep = "")

### Construct Phylogenetic Tree
cat("3) Create phylogenetic tree\n")
# Extract sequences from DADA2 output
sequences<-getSequences(full.seqtab)
names(sequences)<-paste("ASV", stringr::str_pad(1:ncol(full.seqtab),4, pad = "0"), sep = "_")

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
cat("4) Create phyloseq object\n")
# Rename ASVs and samples
asv_table <- full.seqtab
colnames(asv_table) <- paste("ASV", stringr::str_pad(1:ncol(full.seqtab),4, pad = "0"), sep = "_")

#row.names(asv_table) <- paste(opt$RUN_NAME, sample.names, sep = "_")

# Create metadata file
map <- data.frame(
  row.names = row.names(asv_table),
  Tag = stringr::str_split(row.names(asv_table),"\\.",simplify=TRUE)[,2],
  Run = stringr::str_split(row.names(asv_table),"_",simplify=T)[,1],
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
ggsave(filename = "full.ASVs_vs_reads.pdf",
  plot = ggplot(map, aes(x = reads, y = ASVs, color = Run)) + geom_point() + geom_smooth(method=lm, se=FALSE, size=2),
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
save(phy, file = file.path(opt$out_dir,"full.phyloseq_object.RData"))
cat("Final phyloseq object (phy) stored in: ",file.path(opt$out_dir,"full.phyloseq_object.RData"),"\n",sep = "")

q(save="no",status=0)