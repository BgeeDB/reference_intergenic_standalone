## Julien Wollbrett Jun 21, 2019
## This script loads expression at gene level for all libraries of one species, and sums the signal. 

## Usage:
## R CMD BATCH --no-save --no-restore '--args rna_seq_sample_excluded="rna_seq_sample_excluded.txt" all_libraries_dir="path/to/all_libraries" sum_by_species_dir="path/to/output/dir" species_id=9606' sum_by_species.R sum_by_species.Rout
## rna_seq_sample_excluded - one column text file with the name of all libraries not to take into account to calculate the sum by species
## all_libraries_dir    - path to dir containing one dir for each library that has to be summed. Each library dir contains the abundance.tsv file created by kallisto.
## sum_by_species_dir   - dir where summed abundance file is created
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
if (!exists("species_id")){ stop("species_id not defined") }

#####################################################################################

## Create list of libraries for which abundance has to be summed 
sample_excluded <- NULL;
if (exists("rna_seq_sample_excluded")){ 
  if( file.exists(rna_seq_sample_excluded) ){
  	sample_excluded <- read.table(rna_seq_sample_excluded, h=FALSE, sep="\t", comment.char="")
  }
}

if (file.exists(all_libraries_dir)){
  library_dir_name <- as.data.frame(list.dirs(path = all_libraries_dir,
    full.names = FALSE, recursive = FALSE))
  names(library_dir_name) <- c("libraryId")
  if(!is.null(sample_excluded)) {
    library_dir_name <- subset(library_dir_name, subset = !(libraryId %in% sample_excluded$V1))
  }
}

numLibs = 0
for(library_id in library_dir_name$libraryId){
  gene_abundance_file <- file.path(all_libraries_dir, library_id, "gene_level_abundance.tsv")
  if (file.exists(gene_abundance_file)){
    cat("  Summing data from", library_id, "\n")
    numLibs = numLibs + 1
    ## read gene level data
    kallisto_gene_counts <- read.table(gene_abundance_file, h=T, sep="\t")
    if (numLibs == 1){
      summed <- kallisto_gene_counts
    } else {
      summed$counts <- summed$counts + kallisto_gene_counts$counts
      summed$abundance <- summed$abundance + kallisto_gene_counts$abundance
    }
  } else {
    warning("  No data found from", library_id, "\n")
  }
}

## Export file with summed gene abundance data
cat("  Exporting file with summed gene abundance data for species", species_id, "\n")
write.table(summed, file = file.path(sum_by_species_dir, 
    paste0("summed_abundance_", species_id,".tsv")), 
  quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


  