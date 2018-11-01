## Preamble: workspace
setwd("~/Desktop/PhD/Experimental data/Transcriptome differential expression analysis/Workspaces_STAR/")
library(limma)
library(edgeR)
library(calibrate)
library(dplyr)
library(ggplot2)
library(gplots)
library(readr)
library(tibble)
library(magrittr)
library(ggbio)
library(biomaRt)
library(stringr)
library(reshape2)
library(ggrepel)
library(Glimma)
library(pheatmap)
library(GenomicRanges)
library(GenomeInfoDb)
library(biovizBase)
library(Rsubread)

targets_all <- readTargets(file = "Targets.txt")
celltype_all <- factor(targets_all$CellType)
design_all <- model.matrix(~celltype_all)

## Perform read counts on all samples
fc <- featureCounts(files = targets_all$OutputFile, annot.inbuilt = "hg19", isPairedEnd = TRUE, nthreads = 4)
dge_all_unfiltered <- DGEList(counts = fc$counts, genes = fc$annotation[,c("GeneID", "Length")], group = celltype_all)

## Filter lowly expressed genes (at least 10cpm, 2 libraries)
isexpr <- rowSums(cpm(dge_all_unfiltered) > 10) >= 2
dge_all <- dge_all_unfiltered[isexpr,]
dim(dge_all_unfiltered)-dim(dge_all)

## Calculate normalisation factors
dge_all <- calcNormFactors(dge_all)

## Hierarchical clustering dendrogram for pooled samples
dgecpm_all <- cpm(dge_all, log = TRUE)
colnames(dgecpm_all) <- targets_all[,1]
dist_all <- dist(t(dgecpm_all), method = "euclidean")
hclust_all <- as.dendrogram(hclust(dist_all, method = "ward.D2"))
tiff(filename = "hclust_all.tiff", width = 1920, height = 1080, pointsize = 28)
plot(hclust_all, ylab = "Distance measure", main="Unsupervised clustering: all samples", col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", col.axis = "#F38630", type = "triangle")
axis(side = 2, at = seq(0, 150, 10), col = "#F38630", labels = FALSE, lwd = 2)
dev.off()

## MDS plot for pooled samples
tiff(filename = "MDS_all.tiff", width = 1920, height = 1080, pointsize = 28)
plotMDS(dge_all, labels = celltype_all, main="MDS plot: all samples")
dev.off()

## Library size bar plot
tiff(filename = "libsize_all.tiff", width = 960, height = 540, pointsize = 28)
ggplot(dge_all$samples, aes(x = dge_all$samples$group, y = lib.size/1e06, fill = group)) +
  geom_bar(stat = "identity") +
  labs(y="Library Size (millions)", x="Cell line sample", fill="Group") + theme_bw() +
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.position = 'none')
dev.off()

## Load ncbi symbol names
ncbi.L1 <- readLines("/Users/benjamin.leow/Desktop/PhD/Experimental data/Transcriptome differential expression analysis/Workspaces_STAR/Homo_sapiens.gene_info", n=1)
ncbi.colname <- unlist(strsplit(substring(ncbi.L1, 10, 234), ' '))
ncbi <- read.delim("/Users/benjamin.leow/Desktop/PhD/Experimental data/Transcriptome differential expression analysis/Workspaces_STAR/Homo_sapiens.gene_info", skip=1, header=FALSE, stringsAsFactors=FALSE)
colnames(ncbi) <- ncbi.colname

## Annotate gene names
match <- match(dge_all$genes$GeneID, ncbi$GeneID)
dge_all$genes$Chr <- ncbi$chromosome[match]
dge_all$genes$Symbol <- ncbi$Symbol[match]
dge_all$genes$Strand <- NULL




## Preamble: 10/15/25nM Das vs naive/DMSO Ctrls (suffix: _10_15_25, targets: c(1,2,3,6,7), celltype: c(2,2,2,1,1))
targets_10_15_25 <- targets_all[c(1,2,3,6,7),]
targets_10_15_25$CellType <- c(2,2,2,1,1)
celltype_10_15_25 <- factor(targets_10_15_25$CellType)
design_10_15_25 <- model.matrix(~celltype_10_15_25)
dge_10_15_25 <- dge_all[,c(1,2,3,6,7)]
dge_10_15_25$samples$group <- targets_10_15_25$CellType

## MDS plot: 10/15/25nM Das vs naive/DMSO Ctrls
tiff(filename = "MDS_10_15_25.tiff", width = 1920, height = 1080, pointsize = 28)
plotMDS(dge_10_15_25, labels = targets_all$CellType[c(1,2,3,6,7)], main = "MDS plot: 10/15/25nM Das vs naive/DMSO Ctrls")
dev.off()

## Hierarchical clustering dendrogram: 10/15/25nM Das vs naive/DMSO Ctrls
dgecpm_10_15_25 <- cpm(dge_10_15_25, log = TRUE)
colnames(dgecpm_10_15_25) <- targets_all[c(1,2,3,6,7),1]
dist_10_15_25 <- dist(t(dgecpm_10_15_25), method = "euclidean")
hclust_10_15_25 <- as.dendrogram(hclust(dist_10_15_25, method = "ward.D2"))
tiff(filename = "hclust_10_15_25.tiff", width = 1920, height = 1080, pointsize = 28)
plot(hclust_10_15_25, ylab = "Distance measure", main="Unsupervised clustering: 10/15/25nM Das vs naive/DMSO Ctrls", col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", col.axis = "#F38630", type = "triangle")
axis(side = 2, at = seq(0, 150, 10), col = "#F38630", labels = FALSE, lwd = 2)
dev.off()

## Differentially expressed genes: 10/15/25nM Das vs naive/DMSO Ctrls
dge_10_15_25 <- estimateDisp(dge_10_15_25, design = design_10_15_25)
fit_10_15_25 <- glmQLFit(dge_10_15_25, design_10_15_25)
qlf.2vs1_10_15_25 <- glmQLFTest(fit_10_15_25, coef = 2)
test_10_15_25 <- exactTest(dge_10_15_25)
tags_10_15_25 <- topTags(test_10_15_25, n=Inf)
write.table(tags_10_15_25, file = "DE_10_15_25.txt", sep = "\t")

GSEA_10_15_25 <- tags_10_15_25$table
GSEA_10_15_25$stat <- ifelse(GSEA_10_15_25$logFC<0, GSEA_10_15_25$PValue*-1, GSEA_10_15_25$PValue)
GSEA_10_15_25 <- GSEA_10_15_25[,c(4,9)]
GSEA_10_15_25$stat <- 1/GSEA_10_15_25$stat
write.table(GSEA_10_15_25, file = "GSEA_10_15_25.rnk", sep = "\t", row.names = FALSE, quote=FALSE)

logFC_10_15_25 <- tags_10_15_25$table
logFC_10_15_25 <- logFC_10_15_25[,c(4,5)]
write.table(logFC_10_15_25, file = "logFC_10_15_25.rnk", sep = "\t", row.names = FALSE, quote=FALSE)

## Volcano plot: 10/15/25nM Das vs naive/DMSO Ctrls
tiff(filename = "volcano_10_15_25.tiff", width = 1920, height = 1080, pointsize = 28)
with(tags_10_15_25$table, plot(logFC, -log10(PValue), pch=20, cex=0.2, main="Volcano plot: 10/15/25nM Das vs naive/DMSO Ctrls"))
with(subset(tags_10_15_25$table, FDR<0.005), points(logFC, -log10(PValue), pch=20, cex=0.2, col="red"))
with(subset(tags_10_15_25$table, abs(logFC)>2), points(logFC, -log10(PValue), pch=20,cex=0.2, col="orange"))
with(subset(tags_10_15_25$table, FDR<0.005 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20,cex=0.5, col="green"))
with(subset(tags_10_15_25$table, FDR<0.005 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=Symbol, cex=0.5))
dev.off()
write.table(subset(tags_10_15_25$table, FDR<0.005 & abs(logFC)>2), file = "volcano_10_15_25.txt", sep = "\t")

## Heatmap (top diff): 10/15/25nM Das vs naive/DMSO Ctrls
counts_10_15_25 <- as.data.frame(dge_10_15_25$counts)
GeneID <- row.names(counts_10_15_25)
counts_10_15_25 <- cbind(counts_10_15_25, GeneID)
counts_10_15_25 <- merge(tags_10_15_25$table, counts_10_15_25, by="GeneID")
counts_10_15_25 <- counts_10_15_25[order(counts_10_15_25$PValue),]
Symbol <- counts_10_15_25$Symbol
counts_10_15_25$GeneID <- NULL
counts_10_15_25$Length <- NULL
counts_10_15_25$Chr <- NULL
counts_10_15_25$Symbol <- NULL
counts_10_15_25$logFC <- NULL
counts_10_15_25$logCPM <- NULL
counts_10_15_25$PValue <- NULL
counts_10_15_25$FDR <- NULL
counts_10_15_25 <- t(counts_10_15_25)
counts_10_15_25 <- scale(counts_10_15_25)
counts_10_15_25 <- t(counts_10_15_25)
row.names(counts_10_15_25) <- Symbol
tiff(filename = "heatmap_10_15_25.tiff", width = 1920, height = 2160, pointsize = 5)
heatmap.2(counts_10_15_25, col=greenred(75), labCol = dge_all$samples$group[c(1,2,3,6,7)], trace="none", density="none", srtCol = 20, cexCol = 1.5, cexRow = .7, revC = TRUE, main = "Differential expression: 10/15/25nM Das vs naive/DMSO Ctrls")
dev.off()
counts_10_15_25 <- head(counts_10_15_25, n=100L)
tiff(filename = "heatmap100_10_15_25.tiff", width = 1920, height = 2160, pointsize = 30)
heatmap.2(counts_10_15_25, col=greenred(75), labCol = dge_all$samples$group[c(1,2,3,6,7)], trace="none", density="none", srtCol = 20, cexCol = 1.5, cexRow = .7, revC = TRUE, main = "Differential expression: 10/15/25nM Das vs naive/DMSO Ctrls")
dev.off()




## Preamble: 50/200nM Das vs naive/DMSO Ctrls (suffix: _50_200, targets: c(4,5,6,7), celltype: c(2,2,1,1))
targets_50_200 <- targets_all[c(4,5,6,7),]
targets_50_200$CellType <- c(2,2,1,1)
celltype_50_200 <- factor(targets_50_200$CellType)
design_50_200 <- model.matrix(~celltype_50_200)
dge_50_200 <- dge_all[,c(4,5,6,7)]
dge_50_200$samples$group <- targets_50_200$CellType

## MDS plot: 50/200nM Das vs naive/DMSO Ctrls
tiff(filename = "MDS_50_200.tiff", width = 1920, height = 1080, pointsize = 28)
plotMDS(dge_50_200, labels = targets_all$CellType[c(4,5,6,7)], main = "MDS plot: 50/200nM Das vs naive/DMSO Ctrls")
dev.off()

## Hierarchical clustering dendrogram: 50/200nM Das vs naive/DMSO Ctrls
dgecpm_50_200 <- cpm(dge_50_200, log = TRUE)
colnames(dgecpm_50_200) <- targets_all[c(4,5,6,7),1]
dist_50_200 <- dist(t(dgecpm_50_200), method = "euclidean")
hclust_50_200 <- as.dendrogram(hclust(dist_50_200, method = "ward.D2"))
tiff(filename = "hclust_50_200.tiff", width = 1920, height = 1080, pointsize = 28)
plot(hclust_50_200, ylab = "Distance measure", main="Unsupervised clustering: 50/200nM Das vs naive/DMSO Ctrls", col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", col.axis = "#F38630", type = "triangle")
axis(side = 2, at = seq(0, 150, 10), col = "#F38630", labels = FALSE, lwd = 2)
dev.off()

## Differentially expressed genes: 50/200nM Das vs naive/DMSO Ctrls
dge_50_200 <- estimateDisp(dge_50_200, design = design_50_200)
fit_50_200 <- glmQLFit(dge_50_200, design_50_200)
qlf.2vs1_50_200 <- glmQLFTest(fit_50_200, coef = 2)
test_50_200 <- exactTest(dge_50_200)
tags_50_200 <- topTags(test_50_200, n=Inf)
write.table(tags_50_200, file = "DE_50_200.txt", sep = "\t")

GSEA_50_200 <- tags_50_200$table
GSEA_50_200$stat <- ifelse(GSEA_50_200$logFC<0, GSEA_50_200$PValue*-1, GSEA_50_200$PValue)
GSEA_50_200 <- GSEA_50_200[,c(4,9)]
GSEA_50_200$stat <- 1/GSEA_50_200$stat
write.table(GSEA_50_200, file = "GSEA_50_200.rnk", sep = "\t", row.names = FALSE, quote=FALSE)

logFC_50_200 <- tags_50_200$table
logFC_50_200 <- logFC_50_200[,c(4,5)]
write.table(logFC_50_200, file = "logFC_50_200.rnk", sep = "\t", row.names = FALSE, quote=FALSE)

## Volcano plot: 50/200nM Das vs naive/DMSO Ctrls
tiff(filename = "volcano_50_200.tiff", width = 1920, height = 1080, pointsize = 28)
with(tags_50_200$table, plot(logFC, -log10(PValue), pch=20, cex=0.2, main="Volcano plot: 50/200nM Das vs naive/DMSO Ctrls"))
with(subset(tags_50_200$table, FDR<0.005), points(logFC, -log10(PValue), pch=20, cex=0.2, col="red"))
with(subset(tags_50_200$table, abs(logFC)>2), points(logFC, -log10(PValue), pch=20,cex=0.2, col="orange"))
with(subset(tags_50_200$table, FDR<0.005 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20,cex=0.5, col="green"))
with(subset(tags_50_200$table, FDR<0.005 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=Symbol, cex=0.5))
dev.off()
write.table(subset(tags_50_200$table, FDR<0.005 & abs(logFC)>2), file = "volcano_50_200.txt", sep = "\t")

## Heatmap (top diff): 50/200nM Das vs naive/DMSO Ctrls
counts_50_200 <- as.data.frame(dge_50_200$counts)
GeneID <- row.names(counts_50_200)
counts_50_200 <- cbind(counts_50_200, GeneID)
counts_50_200 <- merge(tags_50_200$table, counts_50_200, by="GeneID")
counts_50_200 <- counts_50_200[order(counts_50_200$PValue),]
Symbol <- counts_50_200$Symbol
counts_50_200$GeneID <- NULL
counts_50_200$Length <- NULL
counts_50_200$Chr <- NULL
counts_50_200$Symbol <- NULL
counts_50_200$logFC <- NULL
counts_50_200$logCPM <- NULL
counts_50_200$PValue <- NULL
counts_50_200$FDR <- NULL
counts_50_200 <- t(counts_50_200)
counts_50_200 <- scale(counts_50_200)
counts_50_200 <- t(counts_50_200)
row.names(counts_50_200) <- Symbol
tiff(filename = "heatmap_50_200.tiff", width = 1920, height = 2160, pointsize = 5)
heatmap.2(counts_50_200, col=greenred(75), labCol = dge_all$samples$group[c(4,5,6,7)], trace="none", density="none", srtCol = 20, cexCol = 1.5, cexRow = .7, revC = TRUE, main = "Differential expression: 50/200nM Das vs naive/DMSO Ctrls")
dev.off()
counts_50_200 <- head(counts_50_200, n=100L)
tiff(filename = "heatmap100_50_200.tiff", width = 1920, height = 2160, pointsize = 30)
heatmap.2(counts_50_200, col=greenred(75), labCol = dge_all$samples$group[c(4,5,6,7)], trace="none", density="none", srtCol = 20, cexCol = 1.5, cexRow = .7, revC = TRUE, main = "Differential expression: 50/200nM Das vs naive/DMSO Ctrls")
dev.off()




## Preamble: Clone A8/C10/G5 vs naive/DMSO Ctrls (suffix: _A8_C10_G5, targets: c(6,7,8,10,12), celltype: c(1,1,2,2,2))
targets_A8_C10_G5 <- targets_all[c(6,7,8,10,12),]
targets_A8_C10_G5$CellType <- c(1,1,2,2,2)
celltype_A8_C10_G5 <- factor(targets_A8_C10_G5$CellType)
design_A8_C10_G5 <- model.matrix(~celltype_A8_C10_G5)
dge_A8_C10_G5 <- dge_all[,c(6,7,8,10,12)]
dge_A8_C10_G5$samples$group <- targets_A8_C10_G5$CellType

## MDS plot: Clone A8/C10/G5 vs naive/DMSO Ctrls
tiff(filename = "MDS_A8_C10_G5.tiff", width = 1920, height = 1080, pointsize = 28)
plotMDS(dge_A8_C10_G5, labels = targets_all$CellType[c(6,7,8,10,12)], main = "MDS plot: Clone A8/C10/G5 vs naive/DMSO Ctrls")
dev.off()

## Hierarchical clustering dendrogram: Clone A8/C10/G5 vs naive/DMSO Ctrls
dgecpm_A8_C10_G5 <- cpm(dge_A8_C10_G5, log = TRUE)
colnames(dgecpm_A8_C10_G5) <- targets_all[c(6,7,8,10,12),1]
dist_A8_C10_G5 <- dist(t(dgecpm_A8_C10_G5), method = "euclidean")
hclust_A8_C10_G5 <- as.dendrogram(hclust(dist_A8_C10_G5, method = "ward.D2"))
tiff(filename = "hclust_A8_C10_G5.tiff", width = 1920, height = 1080, pointsize = 28)
plot(hclust_A8_C10_G5, ylab = "Distance measure", main="Unsupervised clustering: Clone A8/C10/G5 vs naive/DMSO Ctrls", col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", col.axis = "#F38630", type = "triangle")
axis(side = 2, at = seq(0, 150, 10), col = "#F38630", labels = FALSE, lwd = 2)
dev.off()

## Differentially expressed genes: Clone A8/C10/G5 vs naive/DMSO Ctrls
dge_A8_C10_G5 <- estimateDisp(dge_A8_C10_G5, design = design_A8_C10_G5)
fit_A8_C10_G5 <- glmQLFit(dge_A8_C10_G5, design_A8_C10_G5)
qlf.2vs1_A8_C10_G5 <- glmQLFTest(fit_A8_C10_G5, coef = 2)
test_A8_C10_G5 <- exactTest(dge_A8_C10_G5)
tags_A8_C10_G5 <- topTags(test_A8_C10_G5, n=Inf)
write.table(tags_A8_C10_G5, file = "DE_A8_C10_G5.txt", sep = "\t")

GSEA_A8_C10_G5 <- tags_A8_C10_G5$table
GSEA_A8_C10_G5$stat <- ifelse(GSEA_A8_C10_G5$logFC<0, GSEA_A8_C10_G5$PValue*-1, GSEA_A8_C10_G5$PValue)
GSEA_A8_C10_G5 <- GSEA_A8_C10_G5[,c(4,9)]
GSEA_A8_C10_G5$stat <- 1/GSEA_A8_C10_G5$stat
write.table(GSEA_A8_C10_G5, file = "GSEA_A8_C10_G5.rnk", sep = "\t", row.names = FALSE, quote=FALSE)

logFC_A8_C10_G5 <- tags_A8_C10_G5$table
logFC_A8_C10_G5 <- logFC_A8_C10_G5[,c(4,5)]
write.table(logFC_A8_C10_G5, file = "logFC_A8_C10_G5.rnk", sep = "\t", row.names = FALSE, quote=FALSE)

## Volcano plot: Clone A8/C10/G5 vs naive/DMSO Ctrls
tiff(filename = "volcano_A8_C10_G5.tiff", width = 1920, height = 1080, pointsize = 28)
with(tags_A8_C10_G5$table, plot(logFC, -log10(PValue), pch=20, cex=0.2, main="Volcano plot: Clone A8/C10/G5 vs naive/DMSO Ctrls"))
with(subset(tags_A8_C10_G5$table, FDR<0.005), points(logFC, -log10(PValue), pch=20, cex=0.2, col="red"))
with(subset(tags_A8_C10_G5$table, abs(logFC)>2), points(logFC, -log10(PValue), pch=20,cex=0.2, col="orange"))
with(subset(tags_A8_C10_G5$table, FDR<0.005 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20,cex=0.5, col="green"))
with(subset(tags_A8_C10_G5$table, FDR<0.005 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=Symbol, cex=0.5))
dev.off()
write.table(subset(tags_A8_C10_G5$table, FDR<0.005 & abs(logFC)>2), file = "volcano_A8_C10_G5.txt", sep = "\t")

## Heatmap (top diff): Clone A8/C10/G5 vs naive/DMSO Ctrls
counts_A8_C10_G5 <- as.data.frame(dge_A8_C10_G5$counts)
GeneID <- row.names(counts_A8_C10_G5)
counts_A8_C10_G5 <- cbind(counts_A8_C10_G5, GeneID)
counts_A8_C10_G5 <- merge(tags_A8_C10_G5$table, counts_A8_C10_G5, by="GeneID")
counts_A8_C10_G5 <- counts_A8_C10_G5[order(counts_A8_C10_G5$PValue),]
Symbol <- counts_A8_C10_G5$Symbol
counts_A8_C10_G5$GeneID <- NULL
counts_A8_C10_G5$Length <- NULL
counts_A8_C10_G5$Chr <- NULL
counts_A8_C10_G5$Symbol <- NULL
counts_A8_C10_G5$logFC <- NULL
counts_A8_C10_G5$logCPM <- NULL
counts_A8_C10_G5$PValue <- NULL
counts_A8_C10_G5$FDR <- NULL
counts_A8_C10_G5 <- t(counts_A8_C10_G5)
counts_A8_C10_G5 <- scale(counts_A8_C10_G5)
counts_A8_C10_G5 <- t(counts_A8_C10_G5)
row.names(counts_A8_C10_G5) <- Symbol
tiff(filename = "heatmap_A8_C10_G5.tiff", width = 1920, height = 2160, pointsize = 5)
heatmap.2(counts_A8_C10_G5, col=greenred(75), labCol = dge_all$samples$group[c(6,7,8,10,12)], trace="none", density="none", srtCol = 20, cexCol = 1.5, cexRow = .7, revC = TRUE, main = "Differential expression: Clone A8/C10/G5 vs naive/DMSO Ctrls")
dev.off()
counts_A8_C10_G5 <- head(counts_A8_C10_G5, n=100L)
tiff(filename = "heatmap100_A8_C10_G5.tiff", width = 1920, height = 2160, pointsize = 30)
heatmap.2(counts_A8_C10_G5, col=greenred(75), labCol = dge_all$samples$group[c(6,7,8,10,12)], trace="none", density="none", srtCol = 20, cexCol = 1.5, cexRow = .7, revC = TRUE, main = "Differential expression: Clone A8/C10/G5 vs naive/DMSO Ctrls")
dev.off()




## Preamble: Clone B9/E2 vs naive/DMSO Ctrls (suffix: _B9_E2, targets: c(6,7,9,11), celltype: c(1,1,2,2))
targets_B9_E2 <- targets_all[c(6,7,9,11),]
targets_B9_E2$CellType <- c(1,1,2,2)
celltype_B9_E2 <- factor(targets_B9_E2$CellType)
design_B9_E2 <- model.matrix(~celltype_B9_E2)
dge_B9_E2 <- dge_all[,c(6,7,9,11)]
dge_B9_E2$samples$group <- targets_B9_E2$CellType

## MDS plot: Clone B9/E2 vs naive/DMSO Ctrls
tiff(filename = "MDS_B9_E2.tiff", width = 1920, height = 1080, pointsize = 28)
plotMDS(dge_B9_E2, labels = targets_all$CellType[c(6,7,9,11)], main = "MDS plot: Clone B9/E2 vs naive/DMSO Ctrls")
dev.off()

## Hierarchical clustering dendrogram: Clone B9/E2 vs naive/DMSO Ctrls
dgecpm_B9_E2 <- cpm(dge_B9_E2, log = TRUE)
colnames(dgecpm_B9_E2) <- targets_all[c(6,7,9,11),1]
dist_B9_E2 <- dist(t(dgecpm_B9_E2), method = "euclidean")
hclust_B9_E2 <- as.dendrogram(hclust(dist_B9_E2, method = "ward.D2"))
tiff(filename = "hclust_B9_E2.tiff", width = 1920, height = 1080, pointsize = 28)
plot(hclust_B9_E2, ylab = "Distance measure", main="Unsupervised clustering: Clone B9/E2 vs naive/DMSO Ctrls", col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", col.axis = "#F38630", type = "triangle")
axis(side = 2, at = seq(0, 150, 10), col = "#F38630", labels = FALSE, lwd = 2)
dev.off()

## Differentially expressed genes: Clone B9/E2 vs naive/DMSO Ctrls
dge_B9_E2 <- estimateDisp(dge_B9_E2, design = design_B9_E2)
fit_B9_E2 <- glmQLFit(dge_B9_E2, design_B9_E2)
qlf.2vs1_B9_E2 <- glmQLFTest(fit_B9_E2, coef = 2)
test_B9_E2 <- exactTest(dge_B9_E2)
tags_B9_E2 <- topTags(test_B9_E2, n=Inf)
write.table(tags_B9_E2, file = "DE_B9_E2.txt", sep = "\t")

GSEA_B9_E2 <- tags_B9_E2$table
GSEA_B9_E2$stat <- ifelse(GSEA_B9_E2$logFC<0, GSEA_B9_E2$PValue*-1, GSEA_B9_E2$PValue)
GSEA_B9_E2 <- GSEA_B9_E2[,c(4,9)]
GSEA_B9_E2$stat <- 1/GSEA_B9_E2$stat
write.table(GSEA_B9_E2, file = "GSEA_B9_E2.rnk", sep = "\t", row.names = FALSE, quote=FALSE)

logFC_B9_E2 <- tags_B9_E2$table
logFC_B9_E2 <- logFC_B9_E2[,c(4,5)]
write.table(logFC_B9_E2, file = "logFC_B9_E2.rnk", sep = "\t", row.names = FALSE, quote=FALSE)

## Volcano plot: Clone B9/E2 vs naive/DMSO Ctrls
tiff(filename = "volcano_B9_E2.tiff", width = 1920, height = 1080, pointsize = 28)
with(tags_B9_E2$table, plot(logFC, -log10(PValue), pch=20, cex=0.2, main="Volcano plot: Clone B9/E2 vs naive/DMSO Ctrls"))
with(subset(tags_B9_E2$table, FDR<0.005), points(logFC, -log10(PValue), pch=20, cex=0.2, col="red"))
with(subset(tags_B9_E2$table, abs(logFC)>2), points(logFC, -log10(PValue), pch=20,cex=0.2, col="orange"))
with(subset(tags_B9_E2$table, FDR<0.005 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20,cex=0.5, col="green"))
with(subset(tags_B9_E2$table, FDR<0.005 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=Symbol, cex=0.5))
dev.off()
write.table(subset(tags_B9_E2$table, FDR<0.005 & abs(logFC)>2), file = "volcano_B9_E2.txt", sep = "\t")

## Heatmap (top diff): Clone B9/E2 vs naive/DMSO Ctrls
counts_B9_E2 <- as.data.frame(dge_B9_E2$counts)
GeneID <- row.names(counts_B9_E2)
counts_B9_E2 <- cbind(counts_B9_E2, GeneID)
counts_B9_E2 <- merge(tags_B9_E2$table, counts_B9_E2, by="GeneID")
counts_B9_E2 <- counts_B9_E2[order(counts_B9_E2$PValue),]
Symbol <- counts_B9_E2$Symbol
counts_B9_E2$GeneID <- NULL
counts_B9_E2$Length <- NULL
counts_B9_E2$Chr <- NULL
counts_B9_E2$Symbol <- NULL
counts_B9_E2$logFC <- NULL
counts_B9_E2$logCPM <- NULL
counts_B9_E2$PValue <- NULL
counts_B9_E2$FDR <- NULL
counts_B9_E2 <- t(counts_B9_E2)
counts_B9_E2 <- scale(counts_B9_E2)
counts_B9_E2 <- t(counts_B9_E2)
row.names(counts_B9_E2) <- Symbol
tiff(filename = "heatmap_B9_E2.tiff", width = 1920, height = 2160, pointsize = 5)
heatmap.2(counts_B9_E2, col=greenred(75), labCol = dge_all$samples$group[c(6,7,9,11)], trace="none", density="none", srtCol = 20, cexCol = 1.5, cexRow = .7, revC = TRUE, main = "Differential expression: Clone B9/E2 vs naive/DMSO Ctrls")
dev.off()
counts_B9_E2 <- head(counts_B9_E2, n=100L)
tiff(filename = "heatmap100_B9_E2.tiff", width = 1920, height = 2160, pointsize = 30)
heatmap.2(counts_B9_E2, col=greenred(75), labCol = dge_all$samples$group[c(6,7,9,11)], trace="none", density="none", srtCol = 20, cexCol = 1.5, cexRow = .7, revC = TRUE, main = "Differential expression: Clone B9/E2 vs naive/DMSO Ctrls")
dev.off()




## Preamble: Clone 50/200/B9/E2 vs naive/DMSO Ctrls (suffix: _50_200_B9_E2, targets: c(4,5,6,7,9,11), celltype: c(1,1,2,2,1,1))
targets_50_200_B9_E2 <- targets_all[c(4,5,6,7,9,11),]
targets_50_200_B9_E2$CellType <- c(1,1,2,2,1,1)
celltype_50_200_B9_E2 <- factor(targets_50_200_B9_E2$CellType)
design_50_200_B9_E2 <- model.matrix(~celltype_50_200_B9_E2)
dge_50_200_B9_E2 <- dge_all[,c(4,5,6,7,9,11)]
dge_50_200_B9_E2$samples$group <- targets_50_200_B9_E2$CellType

## MDS plot: Clone B9/E2 vs naive/DMSO Ctrls
tiff(filename = "MDS_50_200_B9_E2.tiff", width = 1920, height = 1080, pointsize = 28)
plotMDS(dge_50_200_B9_E2, labels = targets_all$CellType[c(4,5,6,7,9,11)], main = "MDS plot: Clone B9/E2 vs naive/DMSO Ctrls")
dev.off()

## Hierarchical clustering dendrogram: Clone B9/E2 vs naive/DMSO Ctrls
dgecpm_50_200_B9_E2 <- cpm(dge_50_200_B9_E2, log = TRUE)
colnames(dgecpm_50_200_B9_E2) <- targets_all[c(4,5,6,7,9,11),1]
dist_50_200_B9_E2 <- dist(t(dgecpm_50_200_B9_E2), method = "euclidean")
hclust_50_200_B9_E2 <- as.dendrogram(hclust(dist_50_200_B9_E2, method = "ward.D2"))
tiff(filename = "hclust_50_200_B9_E2.tiff", width = 1920, height = 1080, pointsize = 28)
plot(hclust_50_200_B9_E2, ylab = "Distance measure", main="Unsupervised clustering: Clone B9/E2 vs naive/DMSO Ctrls", col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", col.axis = "#F38630", type = "triangle")
axis(side = 2, at = seq(0, 150, 10), col = "#F38630", labels = FALSE, lwd = 2)
dev.off()

## Differentially expressed genes: Clone B9/E2 vs naive/DMSO Ctrls
dge_50_200_B9_E2 <- estimateDisp(dge_50_200_B9_E2, design = design_50_200_B9_E2)
fit_50_200_B9_E2 <- glmQLFit(dge_50_200_B9_E2, design_50_200_B9_E2)
qlf.2vs1_50_200_B9_E2 <- glmQLFTest(fit_50_200_B9_E2, coef = 2)
test_50_200_B9_E2 <- exactTest(dge_50_200_B9_E2)
tags_50_200_B9_E2 <- topTags(test_50_200_B9_E2, n=Inf)
write.table(tags_50_200_B9_E2, file = "DE_50_200_B9_E2.txt", sep = "\t")

GSEA_50_200_B9_E2 <- tags_50_200_B9_E2$table
GSEA_50_200_B9_E2$stat <- ifelse(GSEA_50_200_B9_E2$logFC<0, GSEA_50_200_B9_E2$PValue*-1, GSEA_50_200_B9_E2$PValue)
GSEA_50_200_B9_E2 <- GSEA_50_200_B9_E2[,c(4,9)]
GSEA_50_200_B9_E2$stat <- 1/GSEA_50_200_B9_E2$stat
write.table(GSEA_50_200_B9_E2, file = "GSEA_50_200_B9_E2.rnk", sep = "\t", row.names = FALSE, quote=FALSE)

logFC_50_200_B9_E2 <- tags_50_200_B9_E2$table
logFC_50_200_B9_E2 <- logFC_50_200_B9_E2[,c(4,5)]
write.table(logFC_50_200_B9_E2, file = "logFC_50_200_B9_E2.rnk", sep = "\t", row.names = FALSE, quote=FALSE)

## Volcano plot: Clone B9/E2 vs naive/DMSO Ctrls
tiff(filename = "volcano_50_200_B9_E2.tiff", width = 1920, height = 1080, pointsize = 28)
with(tags_50_200_B9_E2$table, plot(logFC, -log10(PValue), pch=20, cex=0.2, main="Volcano plot: Clone B9/E2 vs naive/DMSO Ctrls"))
with(subset(tags_50_200_B9_E2$table, FDR<0.005), points(logFC, -log10(PValue), pch=20, cex=0.2, col="red"))
with(subset(tags_50_200_B9_E2$table, abs(logFC)>2), points(logFC, -log10(PValue), pch=20,cex=0.2, col="orange"))
with(subset(tags_50_200_B9_E2$table, FDR<0.005 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20,cex=0.5, col="green"))
with(subset(tags_50_200_B9_E2$table, FDR<0.005 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=Symbol, cex=0.5))
dev.off()
write.table(subset(tags_50_200_B9_E2$table, FDR<0.005 & abs(logFC)>2), file = "volcano_50_200_B9_E2.txt", sep = "\t")

## Heatmap (top diff): Clone B9/E2 vs naive/DMSO Ctrls
counts_50_200_B9_E2 <- as.data.frame(dge_50_200_B9_E2$counts)
GeneID <- row.names(counts_50_200_B9_E2)
counts_50_200_B9_E2 <- cbind(counts_50_200_B9_E2, GeneID)
counts_50_200_B9_E2 <- merge(tags_50_200_B9_E2$table, counts_50_200_B9_E2, by="GeneID")
counts_50_200_B9_E2 <- counts_50_200_B9_E2[order(counts_50_200_B9_E2$PValue),]
Symbol <- counts_50_200_B9_E2$Symbol
counts_50_200_B9_E2$GeneID <- NULL
counts_50_200_B9_E2$Length <- NULL
counts_50_200_B9_E2$Chr <- NULL
counts_50_200_B9_E2$Symbol <- NULL
counts_50_200_B9_E2$logFC <- NULL
counts_50_200_B9_E2$logCPM <- NULL
counts_50_200_B9_E2$PValue <- NULL
counts_50_200_B9_E2$FDR <- NULL
counts_50_200_B9_E2 <- t(counts_50_200_B9_E2)
counts_50_200_B9_E2 <- scale(counts_50_200_B9_E2)
counts_50_200_B9_E2 <- t(counts_50_200_B9_E2)
row.names(counts_50_200_B9_E2) <- Symbol
tiff(filename = "heatmap_50_200_B9_E2.tiff", width = 1920, height = 2160, pointsize = 5)
heatmap.2(counts_50_200_B9_E2, col=greenred(75), labCol = dge_all$samples$group[c(6,7,9,11)], trace="none", density="none", srtCol = 20, cexCol = 1.5, cexRow = .7, revC = TRUE, main = "Differential expression: Clone B9/E2 vs naive/DMSO Ctrls")
dev.off()
counts_50_200_B9_E2 <- head(counts_50_200_B9_E2, n=100L)
tiff(filename = "heatmap100_50_200_B9_E2.tiff", width = 1920, height = 2160, pointsize = 30)
heatmap.2(counts_50_200_B9_E2, col=greenred(75), labCol = dge_all$samples$group[c(6,7,9,11)], trace="none", density="none", srtCol = 20, cexCol = 1.5, cexRow = .7, revC = TRUE, main = "Differential expression: Clone B9/E2 vs naive/DMSO Ctrls")
dev.off()




# Common tags between 


save.image(file = "workspace_all_STAR.RData")




