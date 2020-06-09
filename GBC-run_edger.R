# BiocManager::install(c("edgeR", "limma", "DESeq2", "ctc", "Biobase", "gplots", "ape", "argparse"))
library(edgeR)

base_names <- list.dirs('data')
base_names <-base_names[grepl('simul', base_names)]

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t")
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t")
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t")
  count_df<-t(counts)
  
  # Normalization (Seurat method)
  # myFilteredData <- NormalizeData(count_df, normalization.method = 'LogNormalize', scale.factor = 10000)
  myFilteredData <- count_df
  rv_genes<-which(apply(myFilteredData,1,var)==0)
  rv_genes_names<-rownames(myFilteredData)[rv_genes]
  count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
  geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
  
  y <- DGEList(counts=count_df, group=cellinfo$Group)
  #filtering
  keep <- filterByExpr(y)
  table(keep)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  # y$samples
  
  cellGroup <- factor(cellinfo$Group)
  cellBatch <- factor(cellinfo$Batch)
  design <- model.matrix(~cellBatch+cellGroup)
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design, robust=TRUE, prior.df = 0)
  y$common.dispersion
  design <- model.matrix(~cellBatch+cellGroup)
  rownames(design) <- colnames(y)
  
  fit <- glmQLFit(y, design, robust=TRUE)
  qlf <- glmQLFTest(fit)
  
  FDR<-p.adjust(qlf$table$PValue,method = "BH")
  qlf$table$FDR <- FDR
  
  result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
  rownames(result.table) <- rownames(qlf)
  
  fullDE=rownames(result.table)[which(result.table$adjpvalue<0.05)]
  fullDown=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC>0)]
  fullUp=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC<0)]
  write.table(result.table, file=paste0(base_name,"/all_edger_df_result_table.txt"), sep = "\t")
  write.table(fullDE, file=paste0(base_name,"/all_edger_df_de_genes.txt"), sep = "\t")
  write.table(fullUp, file=paste0(base_name,"/all_edger_df_up_genes.txt"), sep = "\t")
  write.table(fullDown, file=paste0(base_name,"/all_edger_df_down_genes.txt"), sep = "\t")
}

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts_HVG.txt"), sep = "\t")
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t")
  geneinfo<-read.table(file = paste0(base_name,"/geneinfo.txt"), sep = "\t")
  count_df<-t(counts)
  
  # Normalization (Seurat method)
  # myFilteredData <- NormalizeData(count_df, normalization.method = 'LogNormalize', scale.factor = 10000)
  myFilteredData <- count_df
  rv_genes<-which(apply(myFilteredData,1,var)==0)
  rv_genes_names<-rownames(myFilteredData)[rv_genes]
  count_df<-myFilteredData[!(rownames(myFilteredData) %in% rv_genes_names),]
  geneinfo <- geneinfo[!(rownames(geneinfo) %in% rv_genes_names),]
  
  y <- DGEList(counts=count_df, group=cellinfo$Group)
  
  keep <- filterByExpr(y)
  table(keep)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  cellGroup <- factor(cellinfo$Group)
  cellBatch <- factor(cellinfo$Batch)
  design <- model.matrix(~cellBatch+cellGroup)
  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design, robust=TRUE, prior.df = 0)
  y$common.dispersion
  # plotBCV(y)
  
  design <- model.matrix(~cellBatch+cellGroup)
  rownames(design) <- colnames(y)
  
  fit <- glmQLFit(y, design, robust=TRUE)
  qlf <- glmQLFTest(fit)
  
  FDR<-p.adjust(qlf$table$PValue,method = "BH")
  qlf$table$FDR <- FDR
  # summary(decideTests(qlf))
  
  result.table <- data.frame('pvalue' = qlf$table$PValue, 'adjpvalue' = qlf$table$FDR, 'logFC' = qlf$table$logFC)
  rownames(result.table) <- rownames(qlf)
  fullDE=rownames(result.table)[which(result.table$adjpvalue<0.05)]
  fullDown=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC>0)]
  fullUp=rownames(result.table)[which(result.table$adjpvalue<0.05&result.table$logFC<0)]
  write.table(result.table, file=paste0(base_name,"/HVG_edger_df_result_table.txt"), sep = "\t")
  write.table(fullDE, file=paste0(base_name,"/HVG_edger_df_de_genes.txt"), sep = "\t")
  write.table(fullUp, file=paste0(base_name,"/HVG_edger_df_up_genes.txt"), sep = "\t")
  write.table(fullDown, file=paste0(base_name,"/HVG_edger_df_down_genes.txt"), sep = "\t")
}
