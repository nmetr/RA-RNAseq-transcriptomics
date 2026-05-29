# setup ----
setwd('~/Transcriptomics/')

install.packages('BiocManager')
library(BiocManager)
BiocManager::install('DESeq2')
BiocManager::install('KEGGREST')
BiocManager::install('EnhancedVolcano')
BiocManager::install('pathview')
library(DESeq2)
library(KEGGREST)
library(EnhancedVolcano)
library(pathview)

# data import ----
sample_table <- read.csv('sample_table.csv', row.names = 1)
sample_names <- c(row.names(sample_table))

full_counts <- read.table('count_matrix_RA.txt')
colnames(full_counts) <- sample_names

# deseq analysis ----
deseq_dataset <- DESeqDataSetFromMatrix(countData = full_counts,
                                        colData = sample_table,
                                        design = ~ RA)

ds <- DESeq(deseq_dataset)
resultsNames(ds) # checking control group
ds_results <- results(ds)
ds_results

write.table(ds_results, file = 'DESeq_results.csv', row.names = TRUE, col.names = TRUE)

sum(ds_results$padj < 0.05, na.rm = TRUE) # checking [p < 0.05] vs. [p < 0.05 & |fc| > 1]
sum(ds_results$padj < 0.05 & abs(ds_results$log2FoldChange) > 1, na.rm = TRUE)

ds_results_sig <- ds_results[!is.na(ds_results$padj) & ds_results$padj < 0.05, ]

highest_fold_change <- ds_results_sig[order(ds_results_sig$log2FoldChange, decreasing = TRUE), ]
lowest_fold_change <- ds_results_sig[order(ds_results_sig$log2FoldChange, decreasing = FALSE), ]
head(highest_fold_change, n = 11) # determining xlim for plotting
head(lowest_fold_change, n = 11)

# volcano plot ----
EnhancedVolcano(ds_results,
                lab = rownames(ds_results),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05) +
  scale_x_continuous(
    limits = c(-13, 13),
    breaks = seq(-14, 14, by = 2))

dev.copy(png, 'figures/Volcano_plot.png', 
         width = 12,
         height = 10,
         units = 'in',
         res = 300)
dev.off()

# kegg pathway analysis ----
foldchange_df <- data.frame(ds_results$log2FoldChange, row.names = rownames(ds_results))
colnames(foldchange_df) <- 'log2FoldChange'

fc_vector <- ds_results$log2FoldChange
names(fc_vector) <- rownames(ds_results)

pathview(
  gene.data = fc_vector,
  pathway.id = 'hsa05323',
  species = 'hsa',
  gene.idtype = 'SYMBOL',     
  limit = list(gene = 5))

# gene ontology ----
BiocManager::install('goseq')
library(goseq)

sig_data <- ifelse(abs(fc_vector) > 1, 1, 0)

head(fc_vector)
head(sig_data)

BiocManager::install("EnsDb.Hsapiens.v79")  # genes: symbol -> ensembl
library(EnsDb.Hsapiens.v79)
geneSymbols <- names(sig_data)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys = geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
head(geneIDs)

sig_data_ensembl <- sig_data
names(sig_data_ensembl) <- geneIDs$GENEID
sig_data_ensembl <- sig_data_ensembl[!is.na(sig_data_ensembl)]
head(sig_data_ensembl,n=32)

pwf <- nullp(sig_data_ensembl, 'hg19', 'ensGene') # weigh by length

GO <- goseq(pwf, 'hg19', 'ensGene')

head(GO)
nrow(GO) # how many enriched GO terms? - 21429

GO_sig <- GO$category[GO$over_represented_pvalue<0.01] # how many of these are p < 0.01 ? - 83
head(GO_sig)
length(GO_sig)

library(GO.db) # output GO to txt file
capture.output(for(go in GO_sig[1:83]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="gene_ontology_result.txt")

GO_p0.01 <- GO[GO$category %in% GO_sig, ]

# visualization
GO_p0.01_BP <- GO_p0.01[GO_p0.01$ontology == 'BP', ]
GO_p0.01_BP$term <- reorder(GO_p0.01_BP$term, GO_p0.01_BP$numDEInCat / GO_p0.01_BP$numInCat)

ggplot(GO_p0.01_BP, aes(x = numDEInCat / numInCat, y = term)) +
  geom_point(aes(size = -log10(over_represented_pvalue))) +
  labs(
    x = "Gene ratio",
    y = "Biological process",
    size = 'Over-represented\np-value (-log10)') +
  theme_minimal()

dev.copy(png, 'figures/GO_plot.png', 
         width = 12,
         height = 8,
         units = 'in',
         res = 300)
dev.off()
