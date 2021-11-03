## Julien Wollbrett Jun 21, 2019
## This script loads expression at gene level for all libraries of one species, and sums the signal. 

## Usage:
## R CMD BATCH --no-save --no-restore '--args rna_seq_sample_excluded="rna_seq_sample_excluded.txt" all_libraries_dir="path/to/all_libraries" sum_by_species_dir="path/to/output/dir" species_id=9606' sum_by_species.R sum_by_species.Rout
## rna_seq_sample_excluded - one column text file without header containing the name of all libraries not to take into account to calculate the sum by species
## all_libraries_dir       - path to dir containing one dir for each library that has to be summed. Each library dir contains the abundance.tsv file created by kallisto.
## tx2gene_file            - path to the file containing mapping transcript to gene
## gene2biotype_file       - path to the file containing mapping gene to gene biotype
## sum_by_species_dir      - dir where summed abundance file is created
## species_id              - NCBI species id

## Session info
print(sessionInfo())

## reading in arguments provided in command line
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed in command line
## checking if all necessarily arguments were provided properly
if (!exists("all_libraries_dir")){ stop("all_libraries_dir not defined") }
if (!exists("sum_by_species_dir")){ stop("sum_by_species_dir not defined") }
if (!exists("tx2gene_file")){ stop("tx2gene mapping file not defined") }
if (!exists("gene2biotype_file")){ stop("gene2biotype mapping file not defined") }
if (!exists("species_id")){ stop("species_id not defined") }

################################## FUNCTIONS #######################################

## calculate effect_length for each transcript by using weighted.mean...
calculate_effect_length <- function(counts, effec_length){
  myWeightedMean <- c()
  ## transpose the matrix (this means each column is a transcript ID and each row is a library)
  rawcounts <- t(counts)
  raw_effeclength <- t(effec_length)
  
  for (i in 1:ncol(raw_effeclength)) {
    if (sum(rawcounts[,i]) == 0 ){
      ## provide the same weight for a transcript in case the est_count is always zero
      rawcounts[,i] <- ifelse(rawcounts[,i] == 0, 1)
      weightedMean <-  weighted.mean(x=raw_effeclength[,i], w=rawcounts[,i])
      myWeightedMean[i] <- weightedMean
    } else {
      weightedMean <-  weighted.mean(x=raw_effeclength[,i], w=rawcounts[,i])
      myWeightedMean[i] <- weightedMean
    }
  }
  return(myWeightedMean)
}

estCount_to_tpm <- function(est_count, effec_length){
  rate <- log(est_count) - log(effec_length)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

estCount_to_fpkm <- function(est_count, effec_length){
  N <- sum(est_count)
  exp( log(est_count) + log(1e9) - log(effec_length) - log(N) )
}

#####################################################################################

## Create list of libraries for which abundance has to be summed 
sample_excluded <- NULL;
if (exists("rna_seq_sample_excluded")){ 
  if( file.exists(rna_seq_sample_excluded) ){
  	sample_excluded <- read.table(rna_seq_sample_excluded, h=FALSE, sep="\t", comment.char="")
  } else {
    stop("rna_seq_sample_excluded file ", rna_seq_sample_excluded, " does not exist.")
  }
}

if (file.exists(all_libraries_dir)){
  library_dir_name <- as.data.frame(list.dirs(path = all_libraries_dir,
    full.names = FALSE, recursive = FALSE))
  names(library_dir_name) <- c("libraryId")
  if(!is.null(sample_excluded)) {
    library_dir_name <- subset(library_dir_name, subset = !(libraryId %in% sample_excluded[,1]))
  }
} else {
  stop("directory ", all_libraries_dir, " does not exist.")
}

message("Summing data for species ", speciesid)


tx2gene <- read.table(file = file.path(tx2gene_file), header = TRUE, sep = "\t")
colnames(tx2gene) <- c("target_id", "gene_id")
gene2biotype <- read.table(file = file.path(gene2biotype_file), header = TRUE, sep = "\t")
#combine transcript_id, gene_id and gene biotype in same dataframe
tx2gene <- merge(tx2gene, gene2biotype, by.x = "gene_id", by.y = "id")
# create template of data.frame that will contain summed counts and ponderated mean of eff_length
summed <- read.table(file = file.path(kallisto_count_folder, library_dir_name$libraryId[1], "/abundance.tsv"), header = TRUE, sep = "\t")
summed <- merge(x = summed , y = tx2gene, by = "target_id")
#sanity check that transcript are always in the same order in abundance files
summed <- summed[order(summed$target_id),]
## detect number of transcripts and number of libraries to create matrix with proper size
number_transcripts <- nrow(summed)
numberLib <- nrow(library_dir_name)
# create data frame at proper size in order not to add new columns for each library
effec_length_info <- as.data.frame(matrix(nrow=number_transcripts, ncol=numberLib))
raw_counts_info <- as.data.frame(matrix(nrow=number_transcripts, ncol=numberLib))

numLibs = 0
  # parse abundance files generated by kallisto for each library and generate data 
  # frames containing effective length and raw counts for each 
for(library_id in library_dir_name$libraryId){


  # For each library use the file at transcriptID level
  file <- file.path(kallisto_count_folder, libraryId, "abundance.tsv")
  if (file.exists(file)){
    cat("  Reading the ", libraryId, "\n")
    numLibs = numLibs + 1
    ## read transcript level data for each library
    kallisto_transcript_current <- read.table(file, h=T, sep="\t")
    #sanity check that transcript are always in the same order in abundance files
    kallisto_transcript_current <- kallisto_transcript_current[order(kallisto_transcript_current$target_id),]
    raw_counts_info[,numLibs] <- kallisto_transcript_current$est_counts
    effec_length_info[,numLibs] <- kallisto_transcript_current$eff_length
  } else {
    stop("  No data found from", libraryId, "\n")
  }
}

if (numLibs == 1){
  ## if only one library provide output directly from gene level
  file2 <- file.path(kallisto_count_folder, sampleInfo$libraryId[sampleInfo$speciesId == species][1],
    "gene_level_abundance.tsv")
  summed <- read.table(file2, h=T, sep="\t")
} else {
  ## sum the raw counts for all libraries
  summed$est_counts <- apply(raw_counts_info, 1, sum)

  ## each row is the weighted.mean of eff_length for each transcriptID
  summed$eff_length <- as.numeric(as.character(calculate_effect_length(raw_counts_info, effec_length_info)))

  # calculate effect_length for each transcript by using weighted.mean...
  ## re-calculate TPM and FPKM after collect the weighted.mean of eff_length for each transcriptID and after sum the est_count for each transcriptID
  summed$tpm <- as.numeric(as.character(estCount_to_tpm(summed$est_counts, summed$eff_length)))
  summed$fpkm <- as.numeric(as.character(estCount_to_fpkm(summed$est_counts, summed$eff_length)))
   
  ## select intergenic regions
  intergenic_regions <- summed[summed$type == "intergenic", c("target_id", "est_counts", "tpm", "fpkm", "type", "biotype")]
  names(intergenic_regions)[1] <- "gene_id"

  ## Sum est_counts, TPM and FPKM for genic regions
  sumGenic <- summed[summed$type == "genic", c("gene_id", "est_counts", "tpm", "fpkm", "type", "biotype")]
  sumGenic <- aggregate(sumGenic[,2:4], list(sumGenic$gene_id), sum)
  names(sumGenic)[1] <- "gene_id"
  sumGenic$type <-  rep("genic", times=length(sumGenic$gene_id))
  select_genic <- dplyr::filter(summed, type == "genic")
  select_biotype_genic = distinct(select_genic, gene_id, .keep_all= TRUE)
  sumGenic$biotype <- select_biotype_genic$biotype
  ## Final Table with genic and intergenic information
  summed <- rbind(sumGenic, intergenic_regions)
}

## Export file with summed gene abundance data
cat("  Exporting file with summed gene abundance data for species", species_id, "\n")
write.table(summed, file = file.path(sum_by_species_dir, 
    paste0("summed_abundance_", species_id,".tsv")), 
  quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


  