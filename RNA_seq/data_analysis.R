library(DESeq2)
library(readr)
# import data
RSEM_raw_counts <- read.table("F:/biostar/test/stddata__2016_07_15/BLCA/20160715/gdac.broadinstitute.org_BLCA.mRNAseq_Preprocess.Level_3.2016071500.0.0/rnaseq_test/BLCA.uncv2.mRNAseq_raw_counts.txt", sep = '\t', header = TRUE)

# find raw_counts matrix of TN samples
x <- colnames(RSEM_raw_counts)
geneID <- RSEM_raw_counts[, 1] # get out the first column
colnames(RSEM_raw_counts)[1] <- "geneID" # then, name it geneID
x <- x[-1] # get rid of 'geneID' field name
x <- substr(x, 14, 15) # subtracting 14th and 15th character which represent the sample type
x <- as.integer(x)
x <- x %in% (10:19) # TN sample belongs to 10-19
y <- colnames(RSEM_raw_counts)[-1]
y <- y[x]
tumor_matched_normal_counts <- cbind(geneID, subset(RSEM_raw_counts, select = y)) # Bingo

# check
colnames(tumor_matched_normal_counts)

x <- colnames(tumor_matched_normal_counts)[-1]
x

# filter out the sample information, 9th to 12th character
y <- substr(x, 9, 12)

# get tumor and TN counts matrix from original counts matrix
z <- subset(RSEM_raw_counts, select = grep(paste(y, collapse = "|"), colnames(RSEM_raw_counts), value = TRUE)) # check details about grep, use help function!
colnames(z)

tumor_counts <- subset(z, select = setdiff(colnames(z), colnames(tumor_matched_normal_counts)))

# add geneID
tumor_counts <- cbind(geneID, tumor_counts)
colnames(tumor_counts)

# At last, merge two matrices
cts_BLCA <- cbind(tumor_matched_normal_counts, tumor_counts[-1])
colnames(cts_BLCA)

# remove some mediate variables and execute garbage collection
rm(x, y, z)
gc()

# Data is tidy now!

rownames(cts_BLCA) <- cts_BLCA$geneID
cts_BLCA <- cts_BLCA[, -1]
cts_BLCA[1:6, 1:3]

cts_BLCA <- round(cts_BLCA)
cts_BLCA[2, 1]

dds_BLCA <- DESeqDataSetFromMatrix(countData = cts_BLCA, colData = coldata_BLCA, design = ~condition)

dds_BLCA

# reset the level of `condition` (factor type), to put `tumor_matched_normal` as first level. Normally we shall put control group in the first level to facilitate the following analysis by DESeq2.
levels(dds_BLCA$condition)

dds_BLCA$condition <- relevel(dds_BLCA$condition, ref = "tumor_matched_normal")
levels(dds_BLCA$condition)

dds_BLCA <- DESeq(dds_BLCA) # core function
results_BLCA <- results(dds_BLCA, alpha = 0.05) # set FDR cutoff as 0.05(5 %)
results_BLCA

topGene <- rownames(results_BLCA)[which.min(results_BLCA$padj)] # find out the gene which has the lowest padj value
plotCounts(dds_BLCA, gene = topGene, intgroup = c("condition"))

bottomGene <- rownames(results_BLCA)[which.max(results_BLCA$padj)] # find out the gene which has the highest padj value
plotCounts(dds_BLCA, gene = bottomGene, intgroup = c("condition"))

# import the first 100 genes which have the lowest padj values
results_BLCA_ordered <- results_BLCA[order(results_BLCA$padj), ] # order
results_ordered_dataframe <- as.data.frame(results_BLCA_ordered)[1:100, ]
write.csv(results_ordered_dataframe, file = "results.csv")
