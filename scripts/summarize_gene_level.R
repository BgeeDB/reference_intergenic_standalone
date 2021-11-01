## Julien Wollbrett Jun 21, 2019
## This script summarize abundance of one library at gene level

## Usage:
## R CMD BATCH --no-save --no-restore '--args annotation_file_path="path/to/annotation.gtf" abundance_kallisto_file_path="path/to/abundance.tsv" abundance_gene_level_dir_path="path/to/output/dir"' summarize_gene_level.R summarize_gene_level.Rout
## annotation_file_path           - path to the annotation file
## abundance_kallisto_file_path   - path to the abundance.tsv file created by kallisto
## tx2gene_file_path              - path to the tx2gene mapping file. This file was created during the prepare_GTF step and map a transcript_id to its gene_id
## gene2biotype_file_path         - path to the gene2biotype mapping file. This file was created during the prepare_GTF step and map a gene_id to its biotype
## abundance_gene_level_dir_path  - path to the folder where the file containing the abundance summarized at gene level will be saved 
## ignore_tx_version              - (optional) Argument allowing to remove transcript version (e.g FBtr0306541.3 will become FBtr0306541)

## Session info
print(sessionInfo())

library(tximport)
library(rtracklayer)
library(GenomicFeatures)

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
if (!exists("annotation_file_path")){ stop("annotation_file_path not defined") }
if (!exists("abundance_kallisto_file_path")){ stop("abundance_kallisto_file_path not defined") }
if (!exists("tx2gene_file_path")){ stop("tx2gene_file_path not defined. This file has been generated in the prepare_GTF step") }
if (!exists("gene2biotype_file_path")){ stop("gene2biotype_file_path not defined. This file has been generated in the prepare_GTF step") }
if (!exists("abundance_gene_level_dir_path")){ stop("abundance_gene_level_dir_path not defined") }


# detect if transcript version has to be removed
ignoreTxVersion <- FALSE
if (exists("ignore_tx_version")) {
  ignoreTxVersion <- TRUE
  message("transcript version will be ignored.")
}

# generate mapping between transcript_id and gene_id
tx2gene <- read.table(tx2gene_file_path, sep="\t", header = TRUE)

# generate the mapping between gene_id and biotype (also contains the type column)
gene2biotype <- read.table(gene2biotype_file_path, sep="\t", header = TRUE)


# use tximport to generate the transcript to gene mapping
tximportObject <- tximport(abundance_kallisto_file_path, type = "kallisto", tx2gene = tx2gene, 
                txOut = FALSE, ignoreTxVersion = ignoreTxVersion)
tx_df <- as.data.frame(tximportObject)
tx_df$id <- rownames(tx_df)
tx_df$countsFromAbundance <- NULL
abundance <- merge(tx_df, gene2biotype, by = "id", all = FALSE)
write.table(x = abundance, file = file.path(abundance_gene_level_dir_path, "gene_level_abundance.tsv") , sep = "\t",
          col.names  = TRUE, row.names = FALSE)
