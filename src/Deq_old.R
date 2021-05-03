
################## package prep
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages(c("Rcpp", "readr", "curl", "data.table"))
library(readr)
library(curl)
library(data.table)

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

# making sure col names of cts match row name of coldata
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# http://127.0.0.1:11014/library/DESeq2/doc/DESeq2.html#countmat
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ disease.status)
# can customize the design function
dds

# pre-filtering (eliminate those with fewer reads)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

################## DEQ analysis! 
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)

# rank the rows using r-values 
resOrdered <- res[order(res$pvalue),]

#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

# adjust FDR cutoff, alpha to 0.05
res05 <- results(dds, alpha=0.05)

# RESULT summary
summary(res05)

# log2 fold changes attributable to a given variable over the mean of normalized counts 
#  the fold-change approach, a gene is said to be differentially expressed if 
# the ratio in absolute value of the expression levels between the two conditions 
# exceeds a certain threshold, for example, twofold or threefold change.
plotMA(res, ylim=c(-2,2))

# allow user to interactively click on the plot to identify genes 
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]

# plot counts for each condition 
plotCounts(dds, gene=which.min(res$padj), intgroup="disease.status")

# column description
mcols(res)$description

##################  export results based on statistically difference (p values < 0.1)
resSig <- subset(resOrdered, padj < 0.1)
resSig
write.csv(as.data.frame(resSig), 
          file="disease_status_results.csv")

