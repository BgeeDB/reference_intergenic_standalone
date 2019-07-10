## Julien Wollbrett Jun 24, 2019
## This script is based on the deconvolution script used by the Bgee pipeline and written by Julen Roux.
## It deconvolutes intergenic gaussians and generate :
##    - one file similar to the sum_by_species.tsv file with one more column describing which deconvoluted gaussian each intergenic region corresponds to,
##    - one tsv file called gaussian_choice.tsv where empty fields have to be manually fill in to generate reference intergenique sequences,
##    - one density plot showing density of TPM abundance for genic, protein coding, intergenic and all deconvoluted intergenic gaussians. 


## Usage:
## R CMD BATCH --no-save --no-restore '--args sum_by_species_dir="path/to/sum_by_species.tsv" species_id=species_id output_dir_path="output/dir/path/"' rna_seq_sum_by_species.R rna_seq_sum_by_species.Rout
## sum_by_species_dir  		  - directory where summed abundance information is stored
## species_id                 - NCBI species ID
## output_dir_path            - path to the output dir

## Session info
print(sessionInfo())

## Load library used to deconvolute intergenic gaussians
library(mclust)

## reading in arguments provided in command line
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
    for( i in 1:length(cmd_args) ){
        eval(parse(text=cmd_args[i]))
    }
}

## checking if all necessary arguments were passed in command line
if (!exists("sum_by_species_dir")){ stop("sum_by_species_dir not defined") }
if (!exists("output_dir_path")){ stop("output_dir_path not defined") }
if (!exists("species_id")){ stop("species_id not defined") }

#####################################################################

## write header of file gaussian_choice_by_species.txt
cat("numberGaussiansIntergenic\tselectedGaussianIntergenic\tselectionSideIntergenic\tcomment\n", file = file.path(output_dir_path, "gaussian_choice.txt"), sep = "\t")

## Deconvolute TPM intergenic and genic distributions
## As in Hebenstreit 2011 Mol Syst Biol: use clustering approach
## Mclust: Normal Mixture Modelling for Model-Based Clustering, Classification, and Density Estimation
## We do not chose the number of gaussians, and let mclust choose
cat("  Deconvoluting sub-distributions of intergenic regions\n")

## load summed abundance per species
summed <- read.csv(file.path(sum_by_species_dir,paste0("summed_abundance_", species_id,".tsv")), header = T, sep = "\t")

## Focus on regions with enough signal (remove TPM = 0 or very small)
summed_filtered <- summed[summed$abundance > 10^-6, ]


######################### denstiy plot without deconvolution ###########################

cat("  Plotting density of aggregated data\n")
## Density plot of summed data
pdf(file = paste0(output_dir_path, "/distribution_TPM_sum_", species_id, ".pdf"), width = 6, height = 5)
## par(mar=c(5,6,1,1)) ## bottom, left, top and right margins
## density of log2(TPM) of summed data
dens <- density(log2(na.omit(summed$abundance) + 10^-6))
## Subgroups densities. Visualization trick: we add an invisible set of points at x=-30, to make densities comparable
## genic regions
dens_genic <- density(c(rep(-30, times=sum(summed$type != "genic")), log2(summed$abundance[summed$type == "genic"] + 10^-6)))
## protein-coding genes only (had to take care of NAs strange behavior)
dens_coding <- density(c(rep(-30, times=sum(!summed$biotype %in% "protein_coding")), log2(summed$abundance[summed$biotype %in% "protein_coding"] + 10^-6)))
## intergenic
dens_intergenic <- density(c(rep(-30, times=sum(summed$type != "intergenic")), log2(summed$abundance[summed$type == "intergenic"] + 10^-6)))
## Plot whole distribution
plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-23, 21), lwd=2, main=paste0(as.character(species_id), " (", "TODO", " libraries)"), bty="n", axes=T, xlab="log2(TPM + 10^-6)")
## Add subgroups distributions (genic, intergenic, etc):
## genic
lines(dens_genic, col="firebrick3", lwd=2)
## protein-coding genes
lines(dens_coding, col="firebrick3", lwd=2, lty=2)
## intergenic
lines(dens_intergenic, col="dodgerblue3", lwd=2)
## legend
legend("topleft", c(paste0("all (", length(summed[,1]),")"), paste0("genic (", sum(summed$type == "genic"), ")"), paste0("coding (", sum(summed$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
dev.off()

############################################################################################

## open PDF device
pdf(file = file.path(output_dir_path, paste0("distribution_TPM_sum_deconvolution_", species_id, ".pdf")), width = 6, height = 5)

## Coding regions
mod1 = densityMclust(log2(summed_filtered$abundance[summed_filtered$biotype %in% "protein_coding"]))
plot(mod1, what = "BIC")
cat("    Protein-coding genes:\n")
print(summary(mod1, parameters = TRUE))
plot(mod1, what = "density", data = log2(summed_filtered$abundance[summed_filtered$biotype %in% "protein_coding"]), breaks = 100, xlab="log2(TPM) - protein-coding genes")

## Intergenic regions
mod2 = densityMclust(log2(summed_filtered$abundance[summed_filtered$type == "intergenic"]))
plot(mod2, what = "BIC")
cat("    Intergenic regions:\n")
print(summary(mod2, parameters = TRUE))
plot(mod2, what = "density", data = log2(summed_filtered$abundance[summed_filtered$type == "intergenic"]), breaks = 100, xlab="log2(TPM) - intergenic")

## Plot the density of the original data, and the density of regions classified to different gaussians
cat("  Plotting density of deconvoluted genic and intergenic regions\n")
dens <- density(log2(summed_filtered$abundance))
## Plot whole distribution
plot(dens, ylim=c(0, max(dens$y)*1.1), xlim=c(-7, 20), lwd=2, main=paste0(as.character(species_id), " (", "TODO", " libraries)"), bty="n", axes=T, xlab="log2(TPM)")

## protein-coding genes only (had to take care of NAs strange behavior)
dens_coding <- density(log2(summed_filtered$abundance[summed_filtered$biotype %in% "protein_coding"]))
## Normalize density for number of observations
dens_coding$y <- dens_coding$y * sum(summed_filtered$biotype %in% "protein_coding") / length(summed_filtered$abundance)
lines(dens_coding, col="firebrick3", lwd=2, lty=2)

## intergenic
dens_intergenic <- density(log2(summed_filtered$abundance[summed_filtered$type == "intergenic"]))
dens_intergenic$y <- dens_intergenic$y * sum(summed_filtered$type == "intergenic") / length(summed_filtered$abundance)
lines(dens_intergenic, col="dodgerblue3", lwd=2)
for (i in 1:mod2$G){
    ## if any point classified
    if (sum(mod2$classification == i) >= 2){
        dens_intergenic_sub <- density(log2(summed_filtered$abundance[summed_filtered$type == "intergenic"][mod2$classification == i]))
        ## y-axis scaling
        dens_intergenic_sub$y <- dens_intergenic_sub$y * length(summed_filtered$abundance[summed_filtered$type == "intergenic"][mod2$classification == i]) / length(summed_filtered$abundance)
        lines(dens_intergenic_sub, col=paste0("grey", trunc(100/(mod2$G+1))*i), lwd=2)
        ## Print gaussian number on plot: at location of max value of gaussian
        text(dens_intergenic_sub$x[dens_intergenic_sub$y == max(dens_intergenic_sub$y)], 0.005, labels = i, col=paste0("grey", trunc(100/(mod2$G+1))*i))
    }
}

## legend
legend("topleft", c(paste0("all (", length(summed_filtered[,1]),")"), paste0("coding (", sum(summed_filtered$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed_filtered$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "dodgerblue3"), lty=c(1, 2, 1), bty="n")
dev.off()

## Export file with summed data and classification of intergenic regions
cat("  Exporting aggregated data and classification of intergenic regions\n")
## Add new column to summed object
summed$classification <- NA
summed$classification[summed$abundance > 10^-6 & summed$biotype %in% "protein_coding"] <- paste("coding_", mod1$classification, sep="")
summed$classification[summed$abundance > 10^-6 & summed$type == "intergenic"] <- paste("intergenic_", mod2$classification, sep="")
write.table(summed, file = file.path(output_dir_path, paste0("sum_abundance_gene_classification_", species_id, ".tsv")), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

