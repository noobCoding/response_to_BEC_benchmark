cutoff = 0.05

base_names <- list.dirs('data')
base_names <-base_names[grepl('dgsp', base_names)]

for(base_name in base_names){
  counts<- read.table(file = paste0(base_name,"/counts.txt"), sep = "\t", fill=T)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  count_df<-counts
  
  count_df[is.na(count_df)] = 0.
  
  nf <- edgeR::calcNormFactors(count_df, method = 'TMM')
  
  mod <- model.matrix(~Group+Batch, data=cellinfo)
  
  voom.data <- limma::voom(count_df, design = mod, lib.size = colSums(count_df) * nf)
  voom.data$genes <- rownames(count_df)
  voom.fitlimma <- limma::lmFit(voom.data, design = mod)
  voom.fitbayes <- limma::eBayes(voom.fitlimma)
  voom.pvalues <- voom.fitbayes$p.value[, 2]
  voom.adjpvalues <- p.adjust(voom.pvalues, method = 'BH')
  voom.logFC <- voom.fitbayes$coefficients[, 2]
  voom.score <- 1 - voom.pvalues
  result.table <- data.frame('pvalue' = voom.pvalues, 'adjpvalue' = voom.adjpvalues, 'logFC' = voom.logFC, 'score' = voom.score)
  rownames(result.table) <- rownames(count_df)
  
  voomDE=rownames(result.table)[which(result.table$adjpvalue<cutoff)]
  voomDown=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC>0)]
  voomUp=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC<0)]
  write.table(result.table, file=paste0(base_name,"/voom_result_table.txt"), sep = "\t")
  write.table(voomDE, file=paste0(base_name,"/voom_de_genes.txt"), sep = "\t")
  write.table(voomUp, file=paste0(base_name,"/voom_up_genes.txt"), sep = "\t")
  write.table(voomDown, file=paste0(base_name,"/voom_down_genes.txt"), sep = "\t")
}

for(base_name in base_names){
  
  counts<- read.table(file = paste0(base_name,"/counts_HVG.txt"), sep = "\t", fill=T)
  counts<-t(counts)
  rownames(counts) <- gsub('.', '-', rownames(counts), fixed = TRUE)
  cellinfo<-read.table(file = paste0(base_name,"/cellinfo.txt"), sep = "\t", fill=T)
  count_df<-counts

  count_df[is.na(count_df)] = 0.
  
  nf <- edgeR::calcNormFactors(count_df, method = 'TMM')
  
  mod <- model.matrix(~Group+Batch, data=cellinfo)
  
  voom.data <- limma::voom(count_df, design = mod, lib.size = colSums(count_df) * nf)
  voom.data$genes <- rownames(count_df)
  voom.fitlimma <- limma::lmFit(voom.data, design = mod)
  voom.fitbayes <- limma::eBayes(voom.fitlimma)
  voom.pvalues <- voom.fitbayes$p.value[, 2]
  voom.adjpvalues <- p.adjust(voom.pvalues, method = 'BH')
  voom.logFC <- voom.fitbayes$coefficients[, 2]
  voom.score <- 1 - voom.pvalues
  result.table <- data.frame('pvalue' = voom.pvalues, 'adjpvalue' = voom.adjpvalues, 'logFC' = voom.logFC, 'score' = voom.score)
  rownames(result.table) <- rownames(count_df)
 
  voomDE=rownames(result.table)[which(result.table$adjpvalue<cutoff)]
  voomDown=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC>0)]
  voomUp=rownames(result.table)[which(result.table$adjpvalue<cutoff&result.table$logFC<0)]
  write.table(result.table, file=paste0(base_name,"/HVG_voom_result_table.txt"), sep = "\t")
  write.table(voomDE, file=paste0(base_name,"/HVG_voom_de_genes.txt"), sep = "\t")
  write.table(voomUp, file=paste0(base_name,"/HVG_voom_up_genes.txt"), sep = "\t")
  write.table(voomDown, file=paste0(base_name,"/HVG_voom_down_genes.txt"), sep = "\t")
}
