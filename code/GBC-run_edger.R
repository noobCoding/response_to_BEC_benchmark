# BiocManager::install('edgeR')
library(edgeR)
library(Seurat)

cutoff = 0.05

base_names <- list.dirs('data')
base_names <-base_names[grepl('dgsp', base_names)]

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t")
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t")
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t")
  count_df<-counts

  # Normalization (Seurat method)
  myFilteredData <- count_df
  rv_genes<-which(apply(myFilteredData,1,var)==0)
  rv_genes_names<-rownames(myFilteredData)[rv_genes]
  count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
  geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
  y <- DGEList(counts=count_df, group=cellinfo$Group)
  y <- calcNormFactors(y)
  
  cellGroup <- factor(cellinfo$Group)
  cellBatch <- factor(cellinfo$Batch)
  design <- model.matrix(~cellBatch+cellGroup)
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design, robust=TRUE)
  y$common.dispersion
  
  design <- model.matrix(~cellBatch+cellGroup)
  rownames(design) <- colnames(y)
  
  fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
  qlf <- glmQLFTest(fit)
  
  FDR<-p.adjust(qlf$table$PValue,method = "BH")
  qlf$table$FDR <- FDR
  
  result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
  rownames(result.table) <- rownames(qlf)
  
  fullDE=rownames(result.table)[which(result.table$adjpvalue<cutoff)]
  fullDown=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC>0)]
  fullUp=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC<0)]
  write.table(result.table, file=paste0(base_name,"/all_edger_result_table.txt"), sep = "\t")
  write.table(fullDE, file=paste0(base_name,"/all_edger_de_genes.txt"), sep = "\t")
  write.table(fullUp, file=paste0(base_name,"/all_edger_up_genes.txt"), sep = "\t")
  write.table(fullDown, file=paste0(base_name,"/all_edger_down_genes.txt"), sep = "\t")
}

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts_HVG.txt"), sep = "\t")
  counts<-t(counts)
  rownames(counts) <- gsub('.', '-', rownames(counts), fixed = TRUE)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t")
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t")
  count_df<-counts
  
  # Normalization (Seurat method)
  myFilteredData <- count_df
  rv_genes<-which(apply(myFilteredData,1,var)==0)
  rv_genes_names<-rownames(myFilteredData)[rv_genes]
  count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
  geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
  
  y <- DGEList(counts=count_df, group=cellinfo$Group)
  y <- calcNormFactors(y)
  
  cellGroup <- factor(cellinfo$Group)
  cellBatch <- factor(cellinfo$Batch)
  design <- model.matrix(~cellBatch+cellGroup)
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design, robust=TRUE, prior.df = 0)  
  y$common.dispersion
  
  design <- model.matrix(~cellBatch+cellGroup)
  rownames(design) <- colnames(y)
  
  fit <- glmQLFit(y, design, robust=TRUE, prior.df = 0)
  qlf <- glmQLFTest(fit)
  
  FDR<-p.adjust(qlf$table$PValue,method = "BH")
  qlf$table$FDR <- FDR
  
  result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
  rownames(result.table) <- rownames(qlf)
  fullDE=rownames(result.table)[which(result.table$adjpvalue<cutoff)]
  fullDown=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC>0)]
  fullUp=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC<0)]
  write.table(result.table, file=paste0(base_name,"/HVG_edger_result_table.txt"), sep = "\t")
  write.table(fullDE, file=paste0(base_name,"/HVG_edger_de_genes.txt"), sep = "\t")
  write.table(fullUp, file=paste0(base_name,"/HVG_edger_up_genes.txt"), sep = "\t")
  write.table(fullDown, file=paste0(base_name,"/HVG_edger_down_genes.txt"), sep = "\t")
}