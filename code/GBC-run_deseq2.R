# BiocManager::install("apeglm")
# BiocManager::install("DESeq2")
library(DESeq2)

#cutoff = 0.05

base_names <- list.dirs('data')
base_names <-base_names[grepl('dgsp', base_names)]

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
  count_df<-counts
  
  myFilteredData <- count_df
  rv_genes<-which(apply(myFilteredData,1,var)==0)
  rv_genes_names<-rownames(myFilteredData)[rv_genes]
  count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
  geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
  count_df <- round(count_df, 0) + 1 # pseudo count prevents every gene includes at least one 0
  
  cellinfo$group <- factor(cellinfo$Group)
  cellinfo$batch <- factor(cellinfo$Batch)
  
  dds <- DESeqDataSetFromMatrix(countData = count_df, colData = cellinfo, design = ~batch+group)
  dds <- DESeq2::DESeq(dds) #, fitType ='mean')
  res <- lfcShrink(dds, coef="group_B_vs_A", type="apeglm", lfcThreshold=0)
  
  result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
  rownames(result.table) <- rownames(dds)
  
  fullDE=rownames(result.table)[which(result.table$adjpvalue<cutoff)]
  fullDown=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC>0)]
  fullUp=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC<0)]
  write.table(result.table, file=paste0(base_name,"/all_deseq2_result_table.txt"), sep = "\t")
  write.table(fullDE, file=paste0(base_name,"/all_deseq2_de_genes.txt"), sep = "\t")
  write.table(fullUp, file=paste0(base_name,"/all_deseq2_up_genes.txt"), sep = "\t")
  write.table(fullDown, file=paste0(base_name,"/all_deseq2_down_genes.txt"), sep = "\t")
}

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts_HVG.txt"), sep = "\t", fill=T)
  counts<-t(counts)
  rownames(counts) <- gsub('.', '-', rownames(counts), fixed = TRUE)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t", fill=T)
  count_df<-counts
  
  count_df[is.na(count_df)] = 0.
  
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
  res <- lfcShrink(dds, coef="group_B_vs_A", type="apeglm", lfcThreshold=0)

  result.table <- data.frame('pvalue' = res$pvalue, 'adjpvalue' = res$padj, 'logFC' = res$log2FoldChange)
  rownames(result.table) <- rownames(dds)
  
  fullDE=rownames(result.table)[which(result.table$adjpvalue<cutoff)]
  fullDown=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC>0)]
  fullUp=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC<0)]
  write.table(result.table, file=paste0(base_name,"/HVG_deseq2_result_table.txt"), sep = "\t")
  write.table(fullDE, file=paste0(base_name,"/HVG_deseq2_de_genes.txt"), sep = "\t")
  write.table(fullUp, file=paste0(base_name,"/HVG_deseq2_up_genes.txt"), sep = "\t")
  write.table(fullDown, file=paste0(base_name,"/HVG_deseq2_down_genes.txt"), sep = "\t")
}
