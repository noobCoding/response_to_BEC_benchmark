# BiocManager::install(c("edgeR", "limma", "DESeq2", "ctc", "Biobase", "gplots", "ape", "argparse"))
# BiocManager::install("apeglm")
library(DESeq2)

base_names <- list.dirs('data')
base_names <-base_names[grepl('dgsp', base_names)]

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t")
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t")
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t")
  # count_df<-t(counts)
  count_df<-counts
  # Normalization (Seurat method)
  # myFilteredData <- NormalizeData(count_df, normalization.method = 'LogNormalize', scale.factor = 10000)
  myFilteredData <- count_df
  rv_genes<-which(apply(myFilteredData,1,var)==0)
  rv_genes_names<-rownames(myFilteredData)[rv_genes]
  count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
  geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
  count_df <- round(count_df, 0) + 1 # pseudo count prevents every gene includes at least one 0
  
  # This code chunk assumes that you have a count matrix called cts and a table of sample information called coldata.
  # The design indicates how to model the samples, here, that we want to measure the effect of the condition,
  # controlling for batch differences. The two factor variables batch and condition should be columns of coldata.
  #
  # dds <- DESeqDataSetFromMatrix(countData = cts,
  #                               colData = coldata,
  #                               design= ~ batch + condition)
  # dds <- DESeq(dds)
  # resultsNames(dds) # lists the coefficients
  # res <- results(dds, name="condition_trt_vs_untrt")
  # # or to shrink log fold changes association with condition:
  # res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
  cellinfo$group <- factor(cellinfo$Group)
  cellinfo$batch <- factor(cellinfo$Batch)
  
  dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo, design = ~batch+group)
  dds <- DESeq2::DESeq(dds) #, fitType ='mean')
  # res <- results(dds, name='Group_Group2_vs_Group1')
  # res <- lfcShrink(dds, coef="group_beta_vs_pan", type="apeglm", lfcThreshold=0)
  res <- lfcShrink(dds, coef="group_B_vs_A", type="apeglm", lfcThreshold=0)
  # res05 <- results(dds, pan=0.05)
  # summary(res05)
  #
  
  # summary(decideTests(qlf))
  # (qlf$table)[which(result.table$adjpvalue<0.05&result.table$logFC<0)]
  # length(row.names(qlf$table[which(qlf$table$logFC>0 && qlf$table$adjP<0.05)]))
  result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
  rownames(result.table) <- rownames(dds)
  
  fullDE=rownames(result.table)[which(result.table$adjpvalue<0.05)]
  fullDown=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC>0)]
  fullUp=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC<0)]
  write.table(result.table, file=paste0(base_name,"/all_deseq2_result_table.txt"), sep = "\t")
  write.table(fullDE, file=paste0(base_name,"/all_deseq2_de_genes.txt"), sep = "\t")
  write.table(fullUp, file=paste0(base_name,"/all_deseq2_up_genes.txt"), sep = "\t")
  write.table(fullDown, file=paste0(base_name,"/all_deseq2_down_genes.txt"), sep = "\t")
}

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts_HVG.txt"), sep = "\t")
  counts<-t(counts)
  rownames(counts) <- gsub('.', '-', rownames(counts), fixed = TRUE)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t")
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t")
  count_df<-counts
  
  count_df[is.na(count_df)] = 0.
  
  # Normalization (Seurat method)
  # myFilteredData <- NormalizeData(count_df, normalization.method = 'LogNormalize', scale.factor = 10000)
  myFilteredData <- count_df
  rv_genes<-which(apply(myFilteredData,1,var)==0)
  rv_genes_names<-rownames(myFilteredData)[rv_genes]
  count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
  geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
  count_df <- round(count_df, 0) + 1 # pseudo count prevents every gene includes at least one 0
  
  cellinfo$group <- factor(cellinfo$Group)
  cellinfo$batch <- factor(cellinfo$Batch)
  
  dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo, design = ~batch+group)
  dds <- DESeq(dds) #, fitType = 'mean')
  resultsNames(dds)
  # res <- results(dds, name='Group_Group2_vs_Group1')
  # res <- lfcShrink(dds, coef="Group_Group2_vs_Group1", type="apeglm")
  res <- lfcShrink(dds, coef="group_B_vs_A", type="apeglm", lfcThreshold=0)
  # res <- lfcShrink(dds, coef="group_beta_vs_pan", type="apeglm", lfcThreshold=0)
  
  # res05 <- results(dds, pan=0.05)
  # summary(res05)
  
  result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
  rownames(result.table) <- rownames(dds)
  
  fullDE=rownames(result.table)[which(result.table$adjpvalue<0.05)]
  fullDown=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC>0)]
  fullUp=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC<0)]
  write.table(result.table, file=paste0(base_name,"/HVG_deseq2_result_table.txt"), sep = "\t")
  write.table(fullDE, file=paste0(base_name,"/HVG_deseq2_de_genes.txt"), sep = "\t")
  write.table(fullUp, file=paste0(base_name,"/HVG_deseq2_up_genes.txt"), sep = "\t")
  write.table(fullDown, file=paste0(base_name,"/HVG_deseq2_down_genes.txt"), sep = "\t")
}


# # BiocManager::install(c("edgeR", "limma", "DESeq2", "ctc", "Biobase", "gplots", "ape", "argparse"))
# # BiocManager::install("apeglm")
# library(DESeq2)
# 
# base_names <- list.dirs('data')
# base_names <-base_names[grepl('pan', base_names)]
# 
# for(base_name in base_names){
#   counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t")
#   cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t")
#   geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t")
#   # count_df<-t(counts)
#   count_df<-counts
# 
#   # Normalization (Seurat method)
#   # myFilteredData <- NormalizeData(count_df, normalization.method = 'LogNormalize', scale.factor = 10000)
#   myFilteredData <- count_df
#   rv_genes<-which(apply(myFilteredData,1,var)==0)
#   rv_genes_names<-rownames(myFilteredData)[rv_genes]
#   count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
#   geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
#   count_df <- round(count_df, 0) + 1 # pseudo count prevents every gene includes at least one 0
# 
#   # This code chunk assumes that you have a count matrix called cts and a table of sample information called coldata.
#   # The design indicates how to model the samples, here, that we want to measure the effect of the condition,
#   # controlling for batch differences. The two factor variables batch and condition should be columns of coldata.
#   #
#   # dds <- DESeqDataSetFromMatrix(countData = cts,
#   #                               colData = coldata,
#   #                               design= ~ batch + condition)
#   # dds <- DESeq(dds)
#   # resultsNames(dds) # lists the coefficients
#   # res <- results(dds, name="condition_trt_vs_untrt")
#   # # or to shrink log fold changes association with condition:
#   # res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
#   cellinfo$group <- factor(cellinfo$Group)
#   cellinfo$batch <- factor(cellinfo$Batch)
# 
#   dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo, design = ~batch+group)
#   dds <- DESeq2::DESeq(dds) #, fitType ='mean')
#   # res <- results(dds, name='Group_Group2_vs_Group1')
#   # res <- lfcShrink(dds, coef="group_beta_vs_pan", type="apeglm", lfcThreshold=0)
#   res <- lfcShrink(dds, coef="group_B_vs_A", type="apeglm", lfcThreshold=0)
#   # res05 <- results(dds, pan=0.05)
#   # summary(res05)
#   #
# 
#   # summary(decideTests(qlf))
#   # (qlf$table)[which(result.table$adjpvalue<0.05&result.table$logFC<0)]
#   # length(row.names(qlf$table[which(qlf$table$logFC>0 && qlf$table$adjP<0.05)]))
#   result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
#   rownames(result.table) <- rownames(dds)
# 
#   fullDE=rownames(result.table)[which(result.table$adjpvalue<0.05)]
#   fullDown=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC>0)]
#   fullUp=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC<0)]
#   write.table(result.table, file=paste0(base_name,"/all_deseq2_result_table.txt"), sep = "\t")
#   write.table(fullDE, file=paste0(base_name,"/all_deseq2_de_genes.txt"), sep = "\t")
#   write.table(fullUp, file=paste0(base_name,"/all_deseq2_up_genes.txt"), sep = "\t")
#   write.table(fullDown, file=paste0(base_name,"/all_deseq2_down_genes.txt"), sep = "\t")
# }
# 
# for(base_name in base_names){
#   counts<- read.table(file = paste0(base_name,"/counts_HVG.txt"), sep = "\t")
#   cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t")
#   geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t")
#   count_df<-t(counts)
# 
#   # Normalization (Seurat method)
#   # myFilteredData <- NormalizeData(count_df, normalization.method = 'LogNormalize', scale.factor = 10000)
#   myFilteredData <- count_df
#   rv_genes<-which(apply(myFilteredData,1,var)==0)
#   rv_genes_names<-rownames(myFilteredData)[rv_genes]
#   count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
#   geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
#   count_df <- round(count_df, 0) + 1 # pseudo count prevents every gene includes at least one 0
# 
#   cellinfo$group <- factor(cellinfo$Group)
#   cellinfo$batch <- factor(cellinfo$Batch)
# 
#   dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo, design = ~batch+group)
#   dds <- DESeq(dds) #, fitType = 'mean')
#   resultsNames(dds)
#   # res <- results(dds, name='Group_Group2_vs_Group1')
#   # res <- lfcShrink(dds, coef="Group_Group2_vs_Group1", type="apeglm")
#   res <- lfcShrink(dds, coef="group_B_vs_A", type="apeglm", lfcThreshold=0)
#   # res <- lfcShrink(dds, coef="group_beta_vs_pan", type="apeglm", lfcThreshold=0)
#   
#   # res05 <- results(dds, pan=0.05)
#   # summary(res05)
# 
#   result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
#   rownames(result.table) <- rownames(dds)
# 
#   fullDE=rownames(result.table)[which(result.table$adjpvalue<0.05)]
#   fullDown=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC>0)]
#   fullUp=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC<0)]
#   write.table(result.table, file=paste0(base_name,"/HVG_deseq2_result_table.txt"), sep = "\t")
#   write.table(fullDE, file=paste0(base_name,"/HVG_deseq2_de_genes.txt"), sep = "\t")
#   write.table(fullUp, file=paste0(base_name,"/HVG_deseq2_up_genes.txt"), sep = "\t")
#   write.table(fullDown, file=paste0(base_name,"/HVG_deseq2_down_genes.txt"), sep = "\t")
# }
