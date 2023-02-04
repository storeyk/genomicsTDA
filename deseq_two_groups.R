## file: deseq_two_groups
## author: Farzana_Nasrin
## date: 3/14/2022
## This code produces the result table with several values for differential gene expression analysis such as log2 fold
## changes, p-values, and adjusted p-values
## Input: a csv file with the raw count data, and a csv file including gene ids, scores, and conditions. 
## Output: a csv file with the results values for all of the genes ordered 
## according to the p-values from smallest to largest

library(DESeq2) # for differential gene expression analysis 


## input the count data file and convert it to a matrix
cts <-  read.csv('plus_vs_healthy_counts.csv') 
ctsmatrix <- as.matrix(cts[,-1])
rownames(ctsmatrix) <- cts[,1]

## input the gene ids and conditions. We assign scores such as -1, 0, and +1
coldata <- read.csv('plus_vs_healthy_coldata.csv',row.names = 1)
coldata <- coldata[,c("scores", "condition")]
coldata$condition <- factor(coldata$condition, c("healthy", "plus"))
coldata$scores <- factor(coldata$scores)


## convert matrix to deseq analysis data format  
deseq_data <- DESeqDataSetFromMatrix(countData = ctsmatrix,
                              colData = coldata,
                              design = ~ condition)
## generate deseq results
deseq_data <- DESeq(deseq_data)
deseq_res <- results(deseq_data)

deseq_res <- results(deseq_data, contrast=c("condition","healthy","plus"))
deseq_res <- results(deseq_data, name="condition_healthy_vs_plus")
resultsNames(deseq_data)

## order gene ids according to the lowest p-values
resOrdered <- deseq_res[order(deseq_res$pvalue),]
summary(deseq_res)
sum(resOrdered$pvalue < 0.05, na.rm=TRUE)



res05 <- results(deseq_data, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)



write.csv(as.data.frame(resOrdered), 
          file="condition_healthy_vs_minus_results.csv")

