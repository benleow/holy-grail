## Preamble
setwd("~/Desktop/mRNAseq analysis/Workspaces/")
library(limma)
library(edgeR)

targets_all <- readTargets(file = "Targets.txt")
celltype_all <- factor(targets_all$CellType)
design_all <- model.matrix(~celltype_all)

## Perform read counts on all samples
fc <- featureCounts(files = targets_all$OutputFile, annot.inbuilt = "hg19", isPairedEnd = TRUE, nthreads = 4)
dge_all <- DGEList(counts = fc$counts, genes = fc$annotation[,c("GeneID", "Length")], group = celltype_all)

## Filter lowly expressed genes (at least 10cpm, 2 libraries)
isexpr <- rowSums(cpm(dge_all) > 10) >= 2
dge_all <- dge_all[isexpr,]

## Calculate normalisation factors
dge_all <- calcNormFactors(dge_all)

## Create hierarchichal clustering dendrogram
dgecpm_all <- cpm(dge_all, log = TRUE)
colnames(dgecpm_all) <- targets_all[,1]
dist_all <- dist(t(dgecpm_all), method = "manhattan")
hclust_all <- hclust(dist_all, method = "ward.D2")
plot(as.dendrogram(hclust_all))

## Create MDS plot of all samples
plotMDS(dge_all, labels = celltype_all)

## Load ncbi symbol names
ncbi.L1 <- readLines("/Users/benjamin.leow/Desktop/mRNAseq analysis/Homo_sapiens.gene_info", n=1)
ncbi.colname <- unlist(strsplit(substring(ncbi.L1, 10, 234), ' '))
ncbi <- read.delim("/Users/benjamin.leow/Desktop/mRNAseq analysis/Homo_sapiens.gene_info", skip=1, header=FALSE, stringsAsFactors=FALSE)
colnames(ncbi) <- ncbi.colname

## Annotate gene names
match <- match(dge_all$genes$GeneID, ncbi$GeneID)
dge_all$genes$Chr <- ncbi$chromosome[match]
dge_all$genes$Symbol <- ncbi$Symbol[match]
dge_all$genes$Strand <- NULL

## Differential expression: 10/15/25nM Das vs naive/DMSO Ctrls (suffix: _10_15_25)
targets_10_15_25 <- targets_all[c(1,2,3,6,7),]
targets_10_15_25$CellType <- c(2,2,2,1,1)
celltype_10_15_25 <- factor(targets_10_15_25$CellType)
design_10_15_25 <- model.matrix(~celltype_10_15_25)

dge_10_15_25 <- dge_all[,c(1,2,3,6,7)]
dge_10_15_25$samples$group <- targets_10_15_25$CellType

ed_10_15_25 <- estimateDisp(dge_10_15_25, design = design_10_15_25)
fit_10_15_25 <- glmQLFit(ed_10_15_25, design_10_15_25)
qlf.2vs1_10_15_25 <- glmQLFTest(fit_10_15_25, coef = 2)
tags_10_15_25 <- topTags(qlf.2vs1_10_15_25, n=Inf)

write.table(tags_10_15_25, file = "de_10_15_25.txt", sep = "\t")

## Samples only dendrogram: 10/15/25nM Das vs naive/DMSO Ctrls (suffix: _10_15_25)
dgecpm_10_15_25 <- cpm(dge_10_15_25, log = TRUE)
colnames(dgecpm_10_15_25) <- targets_all$CellType[c(1,2,3,6,7)]
dist_10_15_25 <- dist(t(dgecpm_10_15_25), method = "manhattan")
hclust_10_15_25 <- hclust(dist_10_15_25, method = "ward.D2")
plot(as.dendrogram(hclust_10_15_25))

## Samples only MDS: 10/15/25nM Das vs naive/DMSO Ctrls (suffix: _10_15_25)
plotMDS(dge_10_15_25, labels = targets_all$CellType[c(1,2,3,6,7)])

## Differential expression: 50/200nM Das vs Ctrls (suffix: _50_200)
targets_50_200 <- targets_all[c(4,5,6,7),]
targets_50_200$CellType <- c(2,2,1,1)
celltype_50_200 <- factor(targets_50_200$CellType)
design_50_200 <- model.matrix(~celltype_50_200)

dge_50_200 <- dge_all[,c(4,5,6,7)]
dge_50_200$samples$group <- targets_50_200$CellType

ed_50_200 <- estimateDisp(dge_50_200, design = design_50_200)
fit_50_200 <- glmQLFit(ed_50_200, design_50_200)
qlf.2vs1_50_200 <- glmQLFTest(fit_50_200, coef = 2)
tags_50_200 <- topTags(qlf.2vs1_50_200, n=Inf)

write.table(tags_50_200, file = "de_50_200.txt", sep = "\t")

## Samples only dendrogram: 50/200nM Das vs Ctrls (suffix: _50_200)
dgecpm_50_200 <- cpm(dge_50_200, log = TRUE)
colnames(dgecpm_50_200) <- targets_all$CellType[c(4,5,6,7)]
dist_50_200 <- dist(t(dgecpm_50_200), method = "manhattan")
hclust_50_200 <- hclust(dist_50_200, method = "ward.D2")
plot(as.dendrogram(hclust_50_200))

## Samples only MDS: 50/200nM Das vs Ctrls (suffix: _50_200)
plotMDS(dge_50_200, labels = targets_all$CellType[c(4,5,6,7)])

## Differential expression: A8, C10, G5 vs Ctrls (suffix: _A8_C10_G5)
targets_A8_C10_G5 <- targets_all[c(6,7,8,10,12),]
targets_A8_C10_G5$CellType <- c(1,1,2,2,2)
celltype_A8_C10_G5 <- factor(targets_A8_C10_G5$CellType)
design_A8_C10_G5 <- model.matrix(~celltype_A8_C10_G5)

dge_A8_C10_G5 <- dge_all[,c(6,7,8,10,12)]
dge_A8_C10_G5$samples$group <- targets_A8_C10_G5$CellType

ed_A8_C10_G5 <- estimateDisp(dge_A8_C10_G5, design = design_A8_C10_G5)
fit_A8_C10_G5 <- glmQLFit(ed_A8_C10_G5, design_A8_C10_G5)
qlf.2vs1_A8_C10_G5 <- glmQLFTest(fit_A8_C10_G5, coef = 2)
tags_A8_C10_G5 <- topTags(qlf.2vs1_A8_C10_G5, n=Inf)

write.table(tags_A8_C10_G5, file = "de_A8_C10_G5.txt", sep = "\t")

## Samples only dendrogram: A8, C10, G5 vs Ctrls (suffix: _A8_C10_G5)
dgecpm_A8_C10_G5 <- cpm(dge_A8_C10_G5, log = TRUE)
colnames(dgecpm_A8_C10_G5) <- targets_all$CellType[c(6,7,8,10,12)]
dist_A8_C10_G5 <- dist(t(dgecpm_A8_C10_G5), method = "manhattan")
hclust_A8_C10_G5 <- hclust(dist_A8_C10_G5, method = "ward.D2")
plot(as.dendrogram(hclust_A8_C10_G5))

## Samples only MDS: A8, C10, G5 vs Ctrls (suffix: _A8_C10_G5)
plotMDS(dge_A8_C10_G5, labels = targets_all$CellType[c(6,7,8,10,12)])

## Differential expression: B9, E2 vs Ctrls (suffix: _B9_E2)
targets_B9_E2 <- targets_all[c(6,7,9,11),]
targets_B9_E2$CellType <- c(1,1,2,2)
celltype_B9_E2 <- factor(targets_B9_E2$CellType)
design_B9_E2 <- model.matrix(~celltype_B9_E2)

dge_B9_E2 <- dge_all[,c(6,7,9,11)]
dge_B9_E2$samples$group <- targets_B9_E2$CellType

ed_B9_E2 <- estimateDisp(dge_B9_E2, design = design_B9_E2)
fit_B9_E2 <- glmQLFit(ed_B9_E2, design_B9_E2)
qlf.2vs1_B9_E2 <- glmQLFTest(fit_B9_E2, coef = 2)
tags_B9_E2 <- topTags(qlf.2vs1_B9_E2, n=Inf)

write.table(tags_B9_E2, file = "de_B9_E2.txt", sep = "\t")

## Samples only dendrogram: B9, E2 vs Ctrls (suffix: _B9_E2)
dgecpm_B9_E2 <- cpm(dge_B9_E2, log = TRUE)
colnames(dgecpm_B9_E2) <- targets_all$CellType[c(6,7,9,11)]
dist_B9_E2 <- dist(t(dgecpm_B9_E2), method = "manhattan")
hclust_B9_E2 <- hclust(dist_B9_E2, method = "ward.D2")
plot(as.dendrogram(hclust_B9_E2))

## Samples only MDS: B9, E2 vs Ctrls (suffix: _B9_E2)
plotMDS(dge_B9_E2, labels = targets_all$CellType[c(6,7,9,11)])
