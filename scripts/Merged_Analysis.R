#! ~/miniconda3/envs/DF-2022.1/bin/Rscript

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
###                     	   LOAD LIBRARIES                              ###
##############################################################################
cat("0) Load libraries and variables\n")
suppressMessages(suppressWarnings(library(methods)))
suppressMessages(suppressWarnings(library(dada2)))
suppressMessages(suppressWarnings(library(RcppParallel)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(DECIPHER)))
suppressMessages(suppressWarnings(library(phangorn)))
suppressMessages(suppressWarnings(library(phyloseq)))
suppressMessages(suppressWarnings(library(optparse)))

cat("PACKAGES USED:\n",
    "DADA2:", as.character(packageVersion("dada2")), "|",
    "ggplot2:", as.character(packageVersion("ggplot2")), "|", 
    "phyloseq", as.character(packageVersion("phyloseq")), "|",
    "parallel:", as.character(packageVersion("parallel")), "|",
    "Rcpp:", as.character(packageVersion("Rcpp")), "\n",
    "RcppParallel:", as.character(packageVersion("RcppParallel")), "|",
    "DECIPHER:", as.character(packageVersion("DECIPHER")), "|",
    "phangorn:", as.character(packageVersion("phangorn")), "|",
    "optparse:", as.character(packageVersion("optparse")), "\n")
##############################################################################
###                     CREATE DEFAULT ARGUMENT LIST                       ###
##############################################################################
# import and format arguments
option_list = list(
  make_option(c("-s", "--settings_file"), type="character", default=NULL, help='If a settings file is provided here the system variable (ANALYSIS_FILE) will be ignored', metavar="number"),
  make_option(c("-o", "--out_dir"), type="character", default="output", help="File path for all output", metavar="character"),
  make_option(c("-p", "--PROJECT_NAME"), type="character", default=NULL, help="Name of project to merge, files without this variable in name is not included", metavar="character"),
  make_option(c("-r", "--reference_dir"), type="character", default="../DB", help="This file path has data from the prior script and this together",metavar="character"),
  make_option(c("-n", "--nthreads"), type="numeric", default=0, help="The number of threads to use. Special values: 0 - detect available cores and use all except one.", metavar="number")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

  
##############################################################################
###                        IMPORT ARGUMENTS FROM FILE                      ###
##############################################################################
# Import settings from settings file if provided
if (!is.null(opt$settings_file)) {
  tmp <- readLines(opt$settings_file)
  tmp <- tmp[grepl("export",tmp)]
  var.tmp <- stringr::str_split(tmp[!grepl("in_dir",tmp)], " ", simplify = T )[,2]
  vars <- as.data.frame( stringr::str_split(var.tmp, "=", simplify = T ) )
  row.names(vars) <- vars[,1]
  vars$V1 <- NULL
  
  # Fix the out_dir variable
  if (vars["out_dir",] == "$OUT") vars["out_dir",] <- vars["OUT",] 

  # replace settings in opt
  for (i in 1:(length(opt))){
    if (names(opt)[i] != "settings_file") opt[i][[1]] <- ifelse(is.numeric(opt[i][[1]]), as.numeric(vars[names(opt)[i],]),vars[names(opt)[i],])
  }

  rm(tmp, vars)
}

# Remove unwanted characters
for (i in 1:(length(opt))){
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

##############################################################################
###                           IMPORT AND MERGE                             ###
##############################################################################
# Indentify all .rds files in the out_dir
cat("1) merging runs\n")
incl.runs <- list.files(opt$out_dir, pattern="seqtab_out.rds$", full.names=TRUE)

# Subset list to files from project if setting is provided
if (!is.null(opt$PROJECT_NAME)){
  incl.runs <- incl.runs[grepl(opt$PROJECT_NAME, incl.runs)]
}
cat("The following files are being imported: ",incl.runs,sep = "\n")

# Merge runs
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

# Export taxonomy
if (!is.null(opt$PROJECT_NAME)) {
  out_name <- paste(opt$PROJECT_NAME,"merged.taxonomy.txt",sep = ".")
} else out_name <- "merged.taxonomy.txt"

write.table(taxa.plus, file = file.path(opt$out_dir,out_name), row.names = asv.names)
cat("Taxonomy stored in: ",file.path(opt$out_dir,out_name),"\n",sep = "")
} else {
  write.table(taxa.plus, file = file.path(opt$out_dir,"merged.taxonomy.txt"), row.names = asv.names)
  cat("Taxonomy stored in: ",file.path(opt$out_dir,"merged.taxonomy.txt"),"\n",sep = "")
}

### Construct Phylogenetic Tree
cat("3) Create phylogenetic tree\n")
# Extract sequences from DADA2 output
sequences <- getSequences(full.seqtab)
names(sequences) <- asv.names

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
# Rename ASVs
asv_table <- full.seqtab
colnames(asv_table) <- asv.names

# Store samples names
sample.names <- row.names(asv_table)

# Create metadata file
if (ncol(stringr::str_split(sample.names,"\\.",simplify=TRUE)) == 3){
  map <- data.frame(
    row.names = sample.names,
    Sample = stringr::str_split(sample.names,"\\.",simplify=TRUE)[,3],
    Project = stringr::str_split(sample.names,"\\.",simplify=TRUE)[,1],
    Seq_run = stringr::str_split(sample.names,"\\.",simplify=TRUE)[,2],
    reads = rowSums(asv_table),
    ASVs = rowSums(asv_table > 0)
  )
} else {
  map <- data.frame(
    row.names = sample.names,
    reads = rowSums(asv_table),
    ASVs = rowSums(asv_table > 0)
  )
}

phy <- phyloseq(otu_table(asv_table, taxa_are_rows=FALSE), tax_table(taxa.plus),phy_tree(fitGTR$tree), refseq(DNAStringSet(sequences)), sample_data(map))

# Currently, the phylogenetic tree is not rooted Though it is not necessary here, you will need to root the tree if you want to calculate any phylogeny based diversity metrics (like Unifrac)
set.seed(711) # As there is some randomness in this process, setting the seed ensures reproducibility
phy_tree(phy) <- root(phy_tree(phy), sample(taxa_names(phy), 1), resolve.root = TRUE)

##############################################################################
###                              CREATE PLOTS                              ###
##############################################################################
# Create a plot of richness compared to read_count
if (!is.null(opt$PROJECT_NAME)) {
  out_name <- paste(opt$PROJECT_NAME,"merged.ASVs_vs_reads.pdf",sep = ".")
} else out_name <- "merged.ASVs_vs_reads.pdf"

if (sum(c("Seq_run","Project") %in% sample_variables(phy)) == 2) {
  ggsave(filename = out_name,
    plot = ggplot(sample_data(phy), aes(x = reads, y = ASVs, color = Seq_run)) + 
      geom_point() + 
      geom_smooth(method=lm, se=FALSE, size=2) + 
      facet_wrap("Project"),
    device = "pdf",
    path = opt$out_dir,
    width = 18,
    height = 15,
    units = "cm",
    dpi = 300)
  cat("Plot comparing read depth and number of ASVs stored in: ",file.path(opt$out_dir,out_name),"\n",sep = "")
} 

##############################################################################
###                        EXPORT PHYLOSEQ OBJECT                          ###
##############################################################################
if (!is.null(opt$PROJECT_NAME)) {
  out_name <- paste(opt$PROJECT_NAME,"merged.phyloseq_object.RData",sep = ".")
} else out_name <- "merged.phyloseq_object.RData"

save(phy, file = file.path(opt$out_dir,out_name))
cat("Final phyloseq object (phy) stored in: ",file.path(opt$out_dir,out_name),"\n",sep = "")

q(save="no",status=0)
