################## package prep
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages(c("Rcpp", "readr", "curl", "data.table","ggplot2"))
library(readr)
library(curl)
library(data.table)
library("ggplot2")


BiocManager::install("DESeq2")
#browseVignettes("DESeq2") documentation 

################## data prep
# import count matrix and column data
# if you need other forms of data: see here 
# http://www.sthda.com/english/wiki/reading-data-from-txt-csv-files-r-base-functions
my_data <- read.delim(file.choose())
cts <- my_data [, -1]
cts = as.matrix(as.data.frame(lapply(cts, as.numeric)))
rownames(cts) <- my_data[,1]

col <- read.delim(file.choose())
col <- subset (col, select = -c(1:4))
coldata <- col [, -1]
rownames(coldata) <- col[,1]

coldata_glo <- coldata[coldata$region == 'glomerulus',]
cts_glo <- cts[,rownames(coldata_glo)]
coldata_tub <- coldata[coldata$region == 'tubule',]
cts_tub <- cts[,rownames(coldata_tub)]

# making sure col names of cts match row name of coldata
all(rownames(coldata_glo) %in% colnames(cts_glo))
all(rownames(coldata_glo) == colnames(cts_glo))
cts_glo <- cts_glo[, rownames(coldata_glo)]
all(rownames(coldata_glo) == colnames(cts_glo))

all(rownames(coldata_tub) %in% colnames(cts_tub))
all(rownames(coldata_tub) == colnames(cts_tub))
cts_tub <- cts_tub[, rownames(coldata_tub)]
all(rownames(coldata_tub) == colnames(cts_tub))

################## generate DeqSeq data input for glomerulus 
library("DESeq2")
dds_glo <- DESeqDataSetFromMatrix(countData = round(cts_glo),
                                  colData = coldata_glo,
                                  design = ~ disease.status)
# can customize the design function
dds_glo

# pre-filtering (eliminate those with fewer reads)
keep <- rowSums(counts(dds_glo)) >= 10
dds_glo <- dds_glo[keep,]

################## generate DeqSeq data input for tubule 
dds_tub <- DESeqDataSetFromMatrix(countData = round(cts_tub),
                                  colData = coldata_tub,
                                  design = ~ disease.status)
# can customize the design function
dds_tub

# pre-filtering (eliminate those with fewer reads)
keep <- rowSums(counts(dds_tub)) >= 10
dds_tub <- dds_tub[keep,]

################## DEQ analysis! 
dds_glo <- DESeq(dds_glo)
res_glo <- results(dds_glo)
resultsNames(dds_glo)

dds_tub <- DESeq(dds_tub)
res_tub <- results(dds_tub)
resultsNames(dds_tub)

# rank the rows using r-values 
resOrdered_glo <- res_glo[order(res_glo$pvalue),]
resOrdered_tub <- res_tub[order(res_tub$pvalue),]

################## DeSeq2 plots
plotMA(res_tub, ylim=c(-2,2))
plotCounts(dds_tub, gene=which.min(res_tub$padj), intgroup="disease.status")
plotMA(res_glo, ylim=c(-2,2))
plotCounts(dds_glo, gene=which.min(res_glo$padj), intgroup="disease.status")

# allow user to interactively click on the plot to identify genes 
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]

################## export results based on statistically difference (p values < 0.1)
resSig_glo <- subset(resOrdered_glo, padj < 0.1)
resSig_glo
write.csv(as.data.frame(resSig_glo), 
          file="disease_status_results_glo.csv")

resSig_tub <- subset(resOrdered_tub, padj < 0.1)
resSig_tub
write.csv(as.data.frame(resSig_tub), 
          file="disease_status_results_tub.csv")

################## 
deq_glo = read.csv("disease_status_results_glo.csv",header= TRUE, row.names=1)
deq_tub = read.csv("disease_status_results_tub.csv",header= TRUE, row.names=1)
par(mar = rep(2, 4))
plot(deq_glo[1:10,]$pvalue,)
axis(1, at=1:10, labels=rownames(deq_glo[1:10,]))
plot(deq_tub[1:10,]$pvalue,)
axis(1, at=1:10, labels=rownames(deq_tub[1:10,]))


