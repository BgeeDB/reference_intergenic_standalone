# code used for generation of reference intergenic sequence of D. melanogaster (7227)

## 1. install all prerequisite

## 2. Prepare custom GTF
R CMD BATCH --vanilla '--args gene_gtf_path="../input/7227_ensembl_BDGP6_84.gtf" genome_fasta_path="../input/7227_BDGP6_dna_toplevel.fa" N_block_size=100 N_proportion=0.05 output_gtf_file_path="../output/custom_annotation.gtf"' prepare_GTF.R prepare_GTF.Rout 

## 3. Prepare custom transcriptome
# both tophat and gtf_to_fasta have to be in your path
tophat gtf_to_fasta  ../input/7227_BDGP6_dna_toplevel.fa ../output/custom_annotation.gtf ../output/custom_transcriptome.fa
perl -i -pe 's/^>\d+ +/>/' ../output/custom_transcriptome.fa

## 4. create kallisto indexes
kallisto index -k 15 -i ../output/transcriptome_k15.idx ../output/custom_transcriptome.fa
kallisto index -k 31 -i ../output/transcriptome.idx ../output/custom_transcriptome.fa

## 5. process gene expression abundance for each library
kallisto quant -i ../output/transcriptome.idx -o ../input/SRX018864/ --single -l 180 -s 30 --bias ../input/SRX018864/SRR039433.fastq.gz ../input/SRX018864/SRR039434.fastq.gz ../input/SRX018864/SRR039435.fastq.gz 
kallisto quant -i ../output/transcriptome.idx -o ../input/SRX018865/ --single -l 180 -s 30 --bias ../input/SRX018865/SRR039436.fastq.gz ../input/SRX018865/SRR039437.fastq.gz ../input/SRX018865/SRR039438.fastq.gz  
kallisto quant -i ../output/transcriptome_k15.idx -o ../input/SRX018866/ --single -l 180 -s 30 --bias ../input/SRX018866/SRR039439.fastq.gz ../input/SRX018866/SRR039440.fastq.gz ../input/SRX018866/SRR039441.fastq.gz ../input/SRX018866/SRR039442.fastq.gz ../input/SRX018866/SRR039443.fastq.gz ../input/SRX018866/SRR039444.fastq.gz ../input/SRX018866/SRR039445.fastq.gz
kallisto quant -i ../output/transcriptome_k15.idx -o ../input/SRX018867/ --single -l 180 -s 30 --bias ../input/SRX018867/SRR039446.fastq.gz ../input/SRX018867/SRR039447.fastq.gz ../input/SRX018867/SRR039448.fastq.gz ../input/SRX018867/SRR039449.fastq.gz ../input/SRX018867/SRR039450.fastq.gz ../input/SRX018867/SRR039451.fastq.gz ../input/SRX018867/SRR039452.fastq.gz
kallisto quant -i ../output/transcriptome.idx -o ../input/SRX054459/  --bias ../input/SRX054459/SRR166807_1.fastq.gz ../input/SRX054459/SRR166807_2.fastq.gz
kallisto quant -i ../output/transcriptome.idx -o ../input/SRX054460/  --bias ../input/SRX054460/SRR166808_1.fastq.gz ../input/SRX054460/SRR166808_2.fastq.gz
kallisto quant -i ../output/transcriptome.idx -o ../input/SRX054461/  --bias ../input/SRX054461/SRR166809_1.fastq.gz ../input/SRX054461/SRR166809_2.fastq.gz
kallisto quant -i ../output/transcriptome.idx -o ../input/SRX054462/  --bias ../input/SRX054462/SRR166810_1.fastq.gz ../input/SRX054462/SRR166810_2.fastq.gz
kallisto quant -i ../output/transcriptome_k15.idx -o ../input/SRX157603/ --single -l 180 -s 30 --bias ../input/SRX157603/SRR518259.fastq.gz
kallisto quant -i ../output/transcriptome_k15.idx -o ../input/SRX157604/ --single -l 180 -s 30 --bias ../input/SRX157604/SRR518260.fastq.gz
kallisto quant -i ../output/transcriptome_k15.idx -o ../input/SRX323030/ --single -l 180 -s 30 --bias ../input/SRX323030/SRR935153.fastq.gz
kallisto quant -i ../output/transcriptome_k15.idx -o ../input/SRX323031/ --single -l 180 -s 30 --bias ../input/SRX323031/SRR935154.fastq.gz
kallisto quant -i ../output/transcriptome_k15.idx -o ../input/SRX331966/  --bias ../input/SRX331966/SRR948304_1.fastq.gz ../input/SRX331966/SRR948304_2.fastq.gz
kallisto quant -i ../output/transcriptome_k15.idx -o ../input/SRX331967/  --bias ../input/SRX331967/SRR948305_1.fastq.gz ../input/SRX331967/SRR948305_2.fastq.gz

## 6. Summarize abundance at gene level for each library
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX018864/abundance.tsv" abundance_gene_level_dir_path="../input/SRX018864/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX018865/abundance.tsv" abundance_gene_level_dir_path="../input/SRX018865/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX018866/abundance.tsv" abundance_gene_level_dir_path="../input/SRX018866/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX018867/abundance.tsv" abundance_gene_level_dir_path="../input/SRX018867/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX054459/abundance.tsv" abundance_gene_level_dir_path="../input/SRX054459/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX054460/abundance.tsv" abundance_gene_level_dir_path="../input/SRX054460/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX054461/abundance.tsv" abundance_gene_level_dir_path="../input/SRX054461/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX054462/abundance.tsv" abundance_gene_level_dir_path="../input/SRX054462/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX157603/abundance.tsv" abundance_gene_level_dir_path="../input/SRX157603/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX157604/abundance.tsv" abundance_gene_level_dir_path="../input/SRX157604/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX323030/abundance.tsv" abundance_gene_level_dir_path="../input/SRX323030/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX323031/abundance.tsv" abundance_gene_level_dir_path="../input/SRX323031/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX331966/abundance.tsv" abundance_gene_level_dir_path="../input/SRX331966/"' summarize_gene_level.R summarize_gene_level.Rout
R CMD BATCH --no-save --no-restore '--args annotation_file_path="../output/custom_annotation.gtf" abundance_kallisto_file_path="../input/SRX331967/abundance.tsv" abundance_gene_level_dir_path="../input/SRX331967/"' summarize_gene_level.R summarize_gene_level.Rout

## 7. Sum by species
# run in the script folder
R CMD BATCH --no-save --no-restore '--args rna_seq_sample_excluded="../input/not_wanted_lib.tsv" all_libraries_folder="../input" sum_by_species_folder="../output" species_id=7227' sum_by_species.R sum_by_species.Rout

## 8. Deconvolution of intergenic gaussians
# 
R CMD BATCH --no-save --no-restore '--args sum_by_species_file_path="../output/summed_abundance_7227.tsv" output_dir_path="../output" species_id=7227' intergenic_deconvolution.R intergenic_deconvolution.Rout

## 9. manually modify ../output/gaussian_choice.tsv file

## 10. Generate reference and/or non-reference intergenic sequences 
# sun in the script folder
R CMD BATCH --no-save --no-restore '--args species_id=7227 sum_and_classification_file_path="../output/summed_abundance_7227.tsv" gaussian_choice_file_path="../output" transcriptome_file_path="" output_dir_path="../output/"' generate_intergenic.R generated_intergenic.Rout
