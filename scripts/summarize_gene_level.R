## Julien Wollbrett Jun 21, 2019
## This script summarize abundance of one library at gene level

## Usage:
## R CMD BATCH --no-save --no-restore '--args annotation_file_path="path/to/annotation.gtf" abundance_kallisto_file_path="path/to/abundance.tsv" abundance_gene_level_dir_path="path/to/output/dir"' summarize_gene_level.R summarize_gene_level.Rout
## annotation_file_path           - path to the annotation file
## abundance_kallisto_file_path   - path to the abundance.tsv file created by kallisto
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
if (!exists("abundance_gene_level_dir_path")){ stop("abundance_gene_level_dir_path not defined") }

################################# Functions ###############################

create_tx2gene <- function(annotation_file_path) {
  txdb <- makeTxDbFromGFF(file = annotation_file_path)
  k <- biomaRt::keys(txdb, keytype = "TXNAME")
  tx2gene <- as.data.frame(biomaRt::select(txdb, k, "GENEID", "TXNAME"))
  return(tx2gene)
}

create_gene2biotype <- function(annotation_file_path) {
  column_names <- c("id", "biotype", "type")
  annotation_object <- rtracklayer::import(annotation_file_path)
  annotation_df = as.data.frame(annotation_object)
  annotation_gene <- annotation_df[annotation_df$source != "intergenic",]
  annotation_gene <- as.data.frame(unique(cbind(annotation_gene$gene_id, 
                                              annotation_gene$gene_biotype)))
  annotation_gene[, 3] <- "genic"
  names(annotation_gene) <- column_names

  annotation_intergenic <- annotation_df[annotation_df$source == "intergenic",]
  annotation_intergenic <- as.data.frame(unique(cbind(annotation_intergenic$gene_id, 
                                                    annotation_intergenic$gene_biotype)))
  annotation_intergenic[, 3] <- "intergenic"
  names(annotation_intergenic) <- column_names
  gene_to_biotype <- rbind(annotation_gene, annotation_intergenic)
  return(gene_to_biotype)
}

#######################################################################

# detect if transcript version has to be removed
ignoreTxVersion <- FALSE
if (exists("ignore_tx_version")) {
  ignoreTxVersion <- TRUE
}

# generate mapping between transcript_id and gene_id
tx2gene <- create_tx2gene(annotation_file_path)

# generate the mapping between gene_id and biotype (also contains the type column)
gene2biotype <- create_gene2biotype(annotation_file_path)


# use tximport to generate the transcript to gene mapping
tximportObject <- tximport(abundance_kallisto_file_path, type = "kallisto", tx2gene = tx2gene, 
                txOut = FALSE, ignoreTxVersion = ignoreTxVersion)
tx_df <- as.data.frame(tximportObject)
tx_df$id <- rownames(tx_df)
tx_df$countsFromAbundance <- NULL
abundance <- merge(tx_df, gene2biotype, by = "id", all = FALSE)
write.table(x = abundance, file = file.path(abundance_gene_level_dir_path, "gene_level_abundance.tsv") , sep = "\t",
          col.names  = TRUE, row.names = FALSE)
