## Pasilla RNA-seq Differential Expression Analysis
## DESeq2 workflow: count matrix -> differential expression -> visualization

library(DESeq2)
library(ggplot2)
library(pheatmap)

## ---- Load count matrix from featureCounts output ----
counts <- read.table("results/pasilla_counts.txt",
			header = TRUE, row.names = 1, skip = 1)

## first 5 columns are featureCounts metadata (Chr, Start, End, Strand, Length)
counts <- counts[, 6:9]
colnames(counts) <- c("untreated_1", "untreated_2", "treated_1", "treated_2")

## ---- Sample metadata ----
sample_info <- data.frame(
	condition = c("untreated", "untreated", "treated", "treated"),
	row.names = colnames(counts)
)

## ---- Build DESeq2 dataset ----
dds <- DESeqDataSetFromMatrix(countData = counts,
				colData = sample_info,
				design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "untreated")

## ---- Run differential expression analysis ----
dds <- DESeq(dds)

## ---- Filter low-expression genes ----
dds_filtered <- dds[rowSums(counts(dds)) > 1, ]
dds_filtered <- DESeq(dds_filtered)

## ---- Threshold-based hypothesis testing ----
## Instead of post-hoc filtering (testing against log2FC = 0, then manually
## filtering |log2FC| > 1), the null hypothesis itself is shifted to test
## directly against the fold-change threshold. This is the statistically
## proper way to test "is this gene's fold change significantly greater
## than 1 (in either direction)", rather than just "is it different from
## zero and happens to be bigger than 1 after the fact".
res <- results(dds_filtered,
		lfcThreshold = 1,
		altHypothesis = "greaterAbs",
		alpha = 0.05)

## ---- Save full results ----
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
write.csv(res_df, "results/DESeq2_results.csv", row.names = FALSE)

## ---- Identify significant genes ----
## padj already reflects the lfcThreshold-based test above, so only the
## significance cutoff needs to be applied here
res_df$significant <- ifelse(!is.na(res_df$padj) &
				res_df$padj < 0.05,
				"Significant", "Not Significant")

sig_genes <- res_df[res_df$significant == "Significant", ]
write.csv(sig_genes, "results/Significant_Genes.csv", row.names = FALSE)

## ---- Top 20 DEGs by adjusted p-value ----
top20 <- sig_genes[order(sig_genes$padj), ][1:20, ]
write.csv(top20, "results/Top20_DEGs.csv", row.names = FALSE)

## ===================================================
## Volcano Plot
## ===================================================
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
	geom_point(size = 1.5, alpha = 0.7) +
	geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
	geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
	scale_color_manual(values = c("Significant" = "#D7263D", "Not Significant" = "gray70")) +
	labs(title = "Pasilla RNA-seq Volcano Plot",
		x = "Log2 Fold Change",
		y = "-Log10 Adjusted P-value") +
	theme_minimal()

ggsave("results/volcano_plot.png", width = 8, height = 6, dpi = 150)

## ===================================================
## PCA Plot
## ===================================================
vsd <- vst(dds_filtered, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
	geom_point(size = 4) +
	labs(x = paste0("PC1: ", percent_var[1], "% variance"),
		y = paste0("PC2: ", percent_var[2], "% variance"),
		title = "PCA Plot - Pasilla RNA-seq") +
	theme_minimal()

ggsave("results/PCA_plot.png", width = 7, height = 6, dpi = 150)

## ===================================================
## Heatmap of top differentially expressed genes
## ===================================================
top_genes <- rownames(top20)
heatmap_matrix <- assay(vsd)[top_genes, ]
heatmap_matrix <- heatmap_matrix - rowMeans(heatmap_matrix)

annotation_col <- as.data.frame(colData(dds_filtered)["condition"])

pheatmap(heatmap_matrix,
	annotation_col = annotation_col,
	show_rownames = TRUE,
	filename = "results/Heatmap.png",
	width = 7, height = 8)

## ===================================================
## MA Plot
## ===================================================
png("results/MA_plot.png", width = 800, height = 600, res = 150)
plotMA(res, ylim = c(-5, 5), main = "MA Plot - Pasilla RNA-seq")
dev.off()

## ---- Session summary ----
cat("Total genes analyzed:", nrow(res_df), "\n")
cat("Total significant genes:", nrow(sig_genes), "\n")
cat("Upregulated:", sum(sig_genes$log2FoldChange > 0), "\n")
cat("Downregulated:", sum(sig_genes$log2FoldChange < 0), "\n")
