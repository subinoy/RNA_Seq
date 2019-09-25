
library(dplyr)

leaf_dat <- read.csv("leaf_389_sample.csv", header=T, row.names=1)

xylem_dat <- read.csv("xylem_385_sample.csv", header=T, row.names=1)

## Experimenting of differential gene expression from FPKM values
## Just testing how it behaves!!

library(DESeq2)
data <- as.data.frame(cbind(leaf_dat, xylem_dat))
total <- as.data.frame(apply(data, 2, as.integer))
rownames(total) <- rownames(leaf_dat)

xlm <- rep("xylem", 385)
x_row <- colnames(xylem_dat)
xylem_row_cond <- as.data.frame(cbind(x_row,xlm))
names(xylem_row_cond) <- c("sample", "condition")

lf <- rep("leaf", 389)
l_row <- colnames(leaf_dat)
leaf_row_cond <- as.data.frame(cbind(l_row,lf))
names(leaf_row_cond) <- c("sample", "condition")

sample_info <- as.data.frame(rbind(leaf_row_cond, xylem_row_cond))
rownames(sample_info) <- sample_info$sample

colData <- sample_info


dds <- DESeqDataSetFromMatrix(countData = total,colData = colData ,design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res
write.table(res, file="leaf_vs_xylem_diff_expr_may_28.tbl")
result_A_B <- results(dds, contrast=c("condition","leaf","xylem")) ## contrast specifies conditions to be tested

resOrdered <- res[order(res$padj),]
head(resOrdered)
summary(resOrdered)
write.csv(resOrdered, file="ordered_results.csv", row.names=T)
#plotMA(res, main="DESeq2", ylim=c(-2,2))

# d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", returnData=TRUE)
# plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

pdf(file="MA plot of Leaf Vs Xylem.pdf")
plotMA( res, ylim = c(-2, 2), main="MA plot of Leaf Vs Xylem" )
dev.off()

pdf("Dispersion plot.pdf")
plotDispEsts( dds, ylim = c(1e-6, 1e1), main="Dispersion plot" )
dev.off()

pdf("Histogram plot of p value")
hist( res$pvalue, breaks=20, col="grey", main="Histogram plot of p value" )
dev.off()

write.csv(as.data.frame(resOrdered), file="condition_treated_results.csv", row.names=T)

resSig <- subset(resOrdered, padj < 0.1)
resSig
write.csv(resSig, file="results_significant_p_l_0p1.csv")

# par( mfrow = c( 1, 2 ) )
# plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
# plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )

pdf("Volcano plot.pdf")
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20,
               main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
dev.off()

save.image(file="full_Deseq2.RData")
