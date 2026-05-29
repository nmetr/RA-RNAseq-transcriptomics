setwd('~/Transcriptomics/')

install.packages('BiocManager')
library(BiocManager)
BiocManager::install('Rsubread')
BiocManager::install('Rsamtools')
BiocManager::install('DESeq2')
library(Rsubread)
library(Rsamtools)
library(DESeq2)

# mapping ----
buildindex(
  basename = 'reference_index/ref_human',
  reference = 'reference/GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 5000,
  indexSplit = TRUE)

align.ctrl1 <- align(index = 'ref_human', readfile1 = 'raw_data/SRR4785819_1_subset40k.fastq', readfile2 = 'raw_data/SRR4785819_2_subset40k.fastq', output_file = 'alignments/ctrl1.BAM')
align.ctrl2 <- align(index = 'ref_human', readfile1 = 'raw_data/SRR4785820_1_subset40k.fastq', readfile2 = 'raw_data/SRR4785820_2_subset40k.fastq', output_file = 'alignments/ctrl2.BAM')
align.ctrl3 <- align(index = 'ref_human', readfile1 = 'raw_data/SRR4785828_1_subset40k.fastq', readfile2 = 'raw_data/SRR4785828_2_subset40k.fastq', output_file = 'alignments/ctrl3.BAM')
align.ctrl4 <- align(index = 'ref_human', readfile1 = 'raw_data/SRR4785831_1_subset40k.fastq', readfile2 = 'raw_data/SRR4785831_2_subset40k.fastq', output_file = 'alignments/ctrl4.BAM')
align.arth1 <- align(index = 'ref_human', readfile1 = 'raw_data/SRR4785979_1_subset40k.fastq', readfile2 = 'raw_data/SRR4785979_2_subset40k.fastq', output_file = 'alignments/arth1.BAM')
align.arth2 <- align(index = 'ref_human', readfile1 = 'raw_data/SRR4785980_1_subset40k.fastq', readfile2 = 'raw_data/SRR4785980_2_subset40k.fastq', output_file = 'alignments/arth2.BAM')
align.arth3 <- align(index = 'ref_human', readfile1 = 'raw_data/SRR4785986_1_subset40k.fastq', readfile2 = 'raw_data/SRR4785986_2_subset40k.fastq', output_file = 'alignments/arth3.BAM')
align.arth4 <- align(index = 'ref_human', readfile1 = 'raw_data/SRR4785988_1_subset40k.fastq', readfile2 = 'raw_data/SRR4785988_2_subset40k.fastq', output_file = 'alignments/arth4.BAM')

sample_table <- read.csv('sample_table.csv', row.names = 1)
sample_names <- c(row.names(sample_table))

lapply(sample_names, function(s) {sortBam(file = paste0('alignments/', s, '.BAM'), destination = paste0('alignments/', s, '.sorted'))})
lapply(sample_names, function(s) {indexBam(file = paste0('alignments/', s, '.sorted.bam'))})

# reads per gene ----
sample_files <- paste('alignments/', sample_names, '.BAM', sep = '')
sample_files

count_matrix <- featureCounts(
  files = sample_files,
  annot.ext = 'reference/GCF_000001405.40_GRCh38.p14_annotation.gtf',
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = 'gene', 
  GTF.attrType = 'gene_id',
  useMetaFeatures = TRUE
)

counts <- count_matrix$counts
colnames(counts) <- sample_names
head(counts)
str(counts)

write.csv(counts, file = 'results/count_matrix_RA.txt', row.names = TRUE)