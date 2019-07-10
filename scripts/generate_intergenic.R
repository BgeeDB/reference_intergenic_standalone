## Julien Wollbrett Jun 21, 2019
## This script summarize abundance of one library at gene level

## Usage:
## R CMD BATCH --no-save --no-restore '--args species_id=species_id sum_and_classification_file_path="path/to/sum/classification/file.tsv" gaussian_choice_file_path="path/gaussian/choice.tsv" transcriptome_file_path="/path/to/transcriptome/file.fa" output_dir_path="/output/dir/"' generate_intergenic.R generated_intergenic.Rout
## species_id                       - NCBI species ID
## sum_and_classification_file_path - path to the file summing abundance for all libraries and containing deconvoluted gaussian number
## gaussian_choice_file_path        - path to the file gaussian_choice.tsv
## transcriptome_file_path          - path to transcriptome file
## output_dir_path                  - path to the folder where result are written
## generate_non_ref_intergenic      - (optional) boolean allowing to also generate non reference intergenic sequences (FALSE by default)


## Session info
print(sessionInfo())

library(Biostrings)

## reading in arguments provided in command line
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
    for( i in 1:length(cmd_args) ){
        eval(parse(text=cmd_args[i]))
    }
}

## checking if all necessary arguments were passed in command line
if (!exists("species_id")){ stop("species_id not defined") }
if (!exists("sum_and_classification_file_path")){ stop("sum_and_classification_file_path not defined") }
if (!exists("gaussian_choice_file_path")){ stop("gaussian_choice_file_path not defined") }
if (!exists("transcriptome_file_path")){ stop("transcriptome_file_path not defined") }
if (!exists("output_dir_path")){ stop("output_dir_path not defined") }

## checking optional arguments
non_ref_intergenic <- FALSE
if (exists("generate_non_ref_intergenic")) {non_ref_intergenic <- generate_non_ref_intergenic}


#######################################################################

# load transcriptome and gtf
transcriptome <- readDNAStringSet(transcriptome_file_path)

# load gaussian choice information
gaussian_choice <- read.table(file = gaussian_choice_file_path, header = TRUE, sep = "\t")

# load summed abundance and corresponding gaussian classification
sum_and_classification <- read.table(file = sum_and_classification_file_path, header = TRUE, sep = "\t")

if( gaussian_choice$numberGaussiansIntergenic < gaussian_choice$selectedGaussianIntergenic ) {
    stop("The number of selected intergenic gaussians must be lower or equal to the total number of intergenic gaussians")
}
if(nrow(gaussian_choice) != 1) {
    stop("The gaussian_choice.tsv file should only contain 2 rows. one corresponding to header and one corresponding to guassian choice data specific to your species.")
}
# keep all intergenic regions
all_intergenic_regions <- sum_and_classification[sum_and_classification$type == "intergenic",]

# define max abundance for reference intergenic regions
max_intergenic <- max(all_intergenic_regions$abundance[
    all_intergenic_regions$classification %in% paste0("intergenic_", 
                                                      gaussian_choice$selectedGaussianIntergenic)])

# keep only intergenic transcriptome information
names(transcriptome) <- unlist(lapply(strsplit(as.character(names(transcriptome)), " "), '[[', 1))

# generate ref intergenic sequences
ref_intergenic_regions <- all_intergenic_regions[all_intergenic_regions$abundance <= max_intergenic,]
ref_intergenic_sequences <- transcriptome[ref_intergenic_regions$id]
writeXStringSet(x = ref_intergenic_sequences, 
                filepath = file.path(output_dir_path, paste0("ref_intergenic_",species_id,".fa.gz")), 
                compress = TRUE)

# generate non reference intergenic sequences
if (non_ref_intergenic) {
    non_ref_intergenic_regions <- all_intergenic_regions[all_intergenic_regions$abundance > max_intergenic,]
    non_ref_intergenic_sequences <- transcriptome[non_ref_intergenic_regions$id]
    writeXStringSet(x = non_ref_intergenic_sequences, 
                    filepath = file.path(output_dir_path, paste0("non_ref_intergenic_",species_id,".fa.gz")), 
                    compress = TRUE)
}


