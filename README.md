# Protocol to generate reference intergenic sequences:

Generation of reference intergenic sequences is a preliminary step for the generation of [present/absent gene expression calls](https://github.com/BgeeDB/bgee_pipeline/tree/master/pipeline/RNA_Seq#bgee-rna-seq-analysis-pipeline).  
For species present in Bgee these sequences have already been generated and it is possible to reuse them to generate present/absent calls for your own RNA-Seq libraries using the [BgeeCall R package](https://bioconductor.org/packages/release/workflows/html/BgeeCall.html).  
For species absent from Bgee, anyone can generate their own reference intergenic sequences, and then use them to generate present/absent calls. If there is interest, we will establish a platform to share these reference intergenic sequences with the community, allowing other BgeeCall users to reuse reference intergenic sequences. 

Before starting the generation of reference intergenic regions it is important to understand the following:

* One step of this approach sums gene expression level for different libraries. This summed expression is used to detect a subset of reference intergenic region from all intergenic regions automatically generated using both the annotations and the genome. The more libraries you have, the more precise your reference intergenic regions will be. Do not hesitate to use publicly available data in addition to your own, if possible.
* Bgee only integrates healthy and wild type transcriptomic data. The generation of reference intergenic regions has only been tested in this context. We cannot assure you that this approach will work if your libraries come from heterogeneous cancer or knock-out samples.
* The annotation file has to come from Ensembl or should be formatted like Ensembl GTF files (“exon” and “gene” features for column 3, “gene_id”, “transcript_id”, and “gene_biotype” attributes for column 9)
* This approach was not created for single cell RNA-Seq libraries; testing is ongoing.

# 1. Prerequisite
* [kallisto](https://pachterlab.github.io/kallisto/about) (4.4.0 or later)
* [topHat](https://ccb.jhu.edu/software/tophat/index.shtml) (2.1.1 or later): install bowtie2 and tophat; add tophat and gtf_to_fasta to your path
* [R](https://www.r-project.org/) (3.5.0 or later) 
* R packages [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html), [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html), [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html), [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), and [mclust](https://cran.r-project.org/web/packages/mclust/)
* Perl
* R scripts available [here](https://github.com/BgeeDB/reference_intergenic_standalone/tree/master/scripts)

Of course you will also need :
* the genome annotation (GTF file from Ensembl or similar)
* the genome sequence (FASTA file)
* your reads (FASTQ files)

# 2. Generate custom GTF file

Detects a set of regions in the genome without annotations in the GTF file. These regions, called “intergenic regions” will be added to the original gene annotations in a custom GTF file and are defined following these rules :
* at least 500 bp from any gene annotation;
* at least 1000 bp long;
* if longer than 20 kb only keep 10 kb from each side of the middle of the region;
* no long blocks of Ns (we recommend a limit of 100 bp of consecutive Ns), nor sequences containing more than a certain proportion of Ns (we recommend a limit of 5% of Ns).

<b>inputs : </b> annotation, genome, N blocks size, max proportion of N, output dir  
<b>output : </b> custom GTF with intergenic regions  
<b>command line :</b>  
```
R CMD BATCH --vanilla '--args gene_gtf_path="path/to/gtf_file.gtf" genome_fasta_path="path/to/genome_file.fa" N_block_size=100 N_proportion=0.05 output_gtf_path="path/to/output/custom/gtf_file.gtf"' prepare_GTF.R prepare_GTF.Rout
```

# 3. Generate custom transcriptome

Now that the custom GTF file has been created it is possible to use it to generate a custom transcriptome containing both genes and intergenic sequences.  
For this step Bgee uses TopHat gtf_to_fasta but it is also possible to use BedTools.

<b>input : </b> custom GTF, genome  
<b>output : </b> custom transcriptome  
<b>command lines :</b>  
```
gtf_to_fasta  path/to/genome_file.fa path/to/custom/gtf_file.gtf path/to/output/custom/transcriptome.fa  
perl -i -pe 's/^>\d+ +/>/' path/to/output/custom/transcriptome.fa
```

# 4. Create kallisto index

Generates kallisto indexes using the custom transcriptome previously generated.  
Depending on the read length of the RNA-Seq libraries 2 different indexes with 2 different kmer sizes should be created. If the read length of a RNA-Seq library is small it is better to quantify abundance using a transcriptome index with a smaller kmer size.  
Thus if you have only libraries with read lengths shorter than 50 bp, you need to create a transcriptome index with a kmer size of 15 bp; if you have only libraries with read lengths longer than 50 bp, you need to create a transcriptome index with a kmer size of 31 bp; and if you have libraries with both categories of read lengths, you need to create two indexes, with kmer sizes 15 and 31 bp.  

<b>input : </b> custom transcriptome  
<b>output : </b> kallisto transcriptome indexes  
<b>command lines : </b>  
* read size < 50 bp  
```
kallisto index -k 15 -i output/index/transcriptome_k15.idx path/to/custom/transcriptome.fa  
```
* read size >= 50 bp  
```
kallisto index -k 31 -i output/index/transcriptome.idx path/to/custom/transcriptome.fa  
```

# 5. Quantify abundances of transcripts and intergenic regions

Run kallisto quantification on each library with the previously generated index. Depending on the read size of the library the index generated using a kmer size of 15 bp or 31 bp should be used. Furthermore, kallisto should not be run with the same options for single-end and paired-end libraries.  

<b>input : </b> kallisto transcriptome index, RNA-Seq fastq files  
<b>output : </b> transcript and intergenic abundances  
<b>command lines (more information [here](https://pachterlab.github.io/kallisto/manual)) :  </b>  
* single-end read length < 50 bp  
```
kallisto quant -i index/transcriptome_k15.idx -o output/folder --single -l 180 -s 30 --bias path/to/run.fastq.gz  
```
* paired-end read length < 50 bp  
```
kallisto quant -i index/transcriptome_k15.idx -o output/folder  --bias path/to/run1_1.fastq.gz path/to/run1_2.fastq.gz  
```
* single-end read length >= 50 bp  
```
kallisto quant -i index/transcriptome.idx -o output/folder --single -l 180 -s 30 --bias path/to/run.fastq.gz  
```
* paired-end read length >= 50 bp  
```
kallisto quant -i index/transcriptome.idx -o output/folder  --bias path/to/run1_1.fastq.gz path/to/run1_2.fastq.gz  
```

# 6. Summarize abundance at gene level

kallisto generates abundance at the transcript level but in this approach we are interested to abundance at the gene level.
This can easily be achieved for each library using the tximport R package. This step has to be run for each library.
An optional argument called <i>ignore_tx_version</i> allows to ignore the version of the transcripts (part of the id after a dot e.g NM_000014.5)  

<b>input : </b> annotation, gene2biotype file, tx2gene file, abundances files generated by kallisto  
<b>output : </b> abundance at gene level  
<b>command line : </b>  
```
R CMD BATCH --no-save --no-restore '--args annotation_file_path="path/to/custom_annotation.gtf" abundance_kallisto_file_path="path/to/abundance.tsv" abundance_gene_level_dir_path="path/to/output/dir" tx2gene_file_path="path/to/file.tx2gene" gene2biotype_file_path="path/to/file.gene2biotype"' summarize_gene_level.R summarize_gene_level.Rout
```

# 7. Sum by species

Sums for each gene and intergenic region the abundance of reads from all libraries of the species.  
The script provided below assumes that all your libraries are in the same folder. It also assumes that one abundance file generated by kallisto is present in each library folder. It is possible to define a subset of RNA-Seq by creating a one column text file where each row correspond to one library directory that should not be taken into account to generate the sum.  
If your arborescence of folders is different you should adapt the script accordingly.  

<b>input: </b> abundance file of all libraries  
<b>output : </b> one file with summed gene abundances  
<b>command line : </b>  
```
R CMD BATCH --no-save --no-restore '--args rna_seq_sample_excluded="path/to/rna_seq_sample_excluded.txt" all_libraries_dir="path/to/all_libraries/dir/" sum_by_species_dir="path/to/output/dir/" species_id=9606' sum_by_species.R sum_by_species.Rout
```

# 8. Gaussian deconvolution

Deconvolutes the gaussian representing the summed abundance of intergenic sequences, and generate a density plot representing the distribution of their expression abundance.  
This density plot also represents the summed abundance of protein coding genes. This deconvolution is done using the Mclust R package.

<b>input : </b> summed abundance, NCBI species ID  
<b>output : </b> density plot, summed abundance with deconvoluted gaussian number, template of gaussian_choice.tsv file  
<b>command line : </b>  
```
R CMD BATCH --no-save --no-restore '--args sum_by_species_file_path="path/to/file.tsv" output_dir_path="output/dir" species_id=species_id' intergenic_deconvolution.R intergenic_deconvolution.Rout
```

# 9. MANUALLY choose gaussians representing reference intergenic

The density plot (last plot of the file distribution_TPM_sum_deconvolution_<i>speciesId</i>.pdf) created in the previous step is used to manually define which deconvoluted gaussians will be considered as reference intergenic regions.  
The global idea is to keep intergenic regions present in gaussians with low intersection to the gaussian representing protein coding genes.  
Practically this step allows to define a deconvoluted gaussian number, and the side of this gaussian (left or right) corresponding to the TPM threshold under which all intergenic regions will be considered as reference intergenic. There are examples explaining how to select the gaussian number in the [github of Bgee](https://github.com/BgeeDB/bgee_pipeline/tree/master/pipeline/RNA_Seq#presenceabsence-calls).  

<b>input : </b> density plot  
<b>output : </b> gaussian number and side of the gaussian used as threshold  
<b>example : </b>  

This density plot was generated using 14 RNA-Seq libraries of D. melanogaster. Mclust deconvoluted 3 intergenic gaussians. We now need to select which of these gaussians will be used to generate reference intergenic sequences. We decided to <b>keep as reference all intergenic sequences part of gaussian 1 and 2</b>. In order to do that we have to fill the gaussian_choice.tsv file created in the previous step. For the column <b>numberGaussiansIntergenic</b> we write <b>3</b>. For the column <b>selectedGaussianIntergenic</b> we write <b>2</b> and for the column <b>selectionSideIntergenic</b> we write <b>Left</b> meaning that the second deconvoluted gaussian will be at the left of our cutoff (see vertical red line). The last column, called comment should be filled in with the explanation of your gaussian choice (e.g Gaussian 3 intergenic very small and on the right: removed)

# 10. Generate reference intergenic sequences

Generates a FASTA file containing reference intergenic sequences (i.e. with TPM below the TPM threshold defined in the previous step on the sum of libraries).  
In this step it is also possible to generate non reference intergenic sequences. These sequences have a TPM comparable to protein coding sequences but are not present in the annotation files. They might be of interest as potential non annotated genes. To generate these sequences please add the argument <i>generate_non_ref_intergenic=TRUE</i> to the command line.

<b>input : </b> gaussian_choice.tsv file, summed abundance with deconvoluted intergenic gaussian number, custom transcriptome, NCBI species ID  
<b>output : </b> reference intergenic sequences, non reference intergenic sequences (optional)  
<b>command line : </b>  
```
R CMD BATCH --no-save --no-restore '--args species_id=species_id sum_and_classification_file_path="path/to/sum/classification/file.tsv" gaussian_choice_file_path="path/to/gaussian_choice.tsv" transcriptome_file_path="/path/to/custom_transcriptome.fa" output_dir_path="/output/dir/"' generate_intergenic.R generated_intergenic.Rout
```

# 11. Share your reference intergenic sequences (Work in progress)

Now that the reference intergenic sequences file was created it is possible to use it to generate present/absent gene expression calls with BgeeCall.  
As a member of the Bgee community, it is gratifying to let other members be aware of your work. It can even be better to allow the community to reuse your intergenic sequences. The Bgee team is currently testing different public dataset repositories (Zenodo, figshare, ...) in order to easily publish and retrieve the sequences. 


