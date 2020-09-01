
rm(list=ls())

library(SingleCellExperiment)
library(scater)
#BiocManager::install("scMerge")
library(scMerge)
# 

lsdir <- list.dirs('data', recursive=FALSE) 

sapply(lsdir,function(x){
  
  x2 <- gsub('data/','',x)
  dir.create('demo_scmerge')
  dir.create(paste0('demo_scmerge/',x2), showWarnings = FALSE)
  
  selection <- c('all', 'HVG')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_scmerge/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    if(s=='HVG'){
      counts <- read.table(paste0(x,'/counts_HVG.txt'), head=T, sep='\t')
      counts<-t(counts)
      rownames(counts) <- gsub('.', '-', rownames(counts), fixed = TRUE)
    } else {
      counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t')
      counts_HVG <- read.table(paste0(x,'/counts_HVG.txt'), head=T, sep='\t')
      counts_HVG <- t(counts_HVG)
      rownames(counts_HVG) <- gsub('.', '-', rownames(counts_HVG), fixed = TRUE)
    }
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t')
    cellinfo <- subset(cellinfo,select=c('Batch','Group'))
    colnames(cellinfo) <- c('batch','cell_type')
    rownames(cellinfo) <- factor(colnames(counts))
    
    counts <- counts[rowSums(counts)>1,]
    
    # remove genes with variance equals to 0
    rv_genes <- which(apply(counts, 1, var)==0) 
    rv_genes_names <- rownames(counts)[rv_genes]
    count_df <- counts[!(rownames(counts) %in% rv_genes_names),]
    
    count_df <- count_df[,rownames(cellinfo)]
    
    sce <- SingleCellExperiment(assays = list(counts = count_df),colData = cellinfo)
    sce <- logNormCounts(sce)
    
    counts_b1 <- count_df[,colnames(count_df) %in% rownames(cellinfo[cellinfo$batch==unique(cellinfo$batch)[1],])]
    counts_b2 <- count_df[,colnames(count_df) %in% rownames(cellinfo[cellinfo$batch==unique(cellinfo$batch)[2],])]
    data <- list(counts_b1,counts_b2)
    
    geneinfo <- read.table(paste0(x,'/geneinfo.txt'), head=T, sep='\t')
    up_genes <- read.table(paste0(x,'/true_up_genes.txt'), head=T, sep='\t')
    down_genes <- read.table(paste0(x,'/true_down_genes.txt'), head=T, sep='\t')
    gi_wobatch <- rbind(up_genes, down_genes)
    geneinfo_wobatch <- geneinfo$x[gi_wobatch$x]
    counts2 <- count_df[rownames(count_df) %in% as.character(geneinfo_wobatch),]
    gene_lowvar <- names(sort(apply(counts2,1,var), decreasing=F)[1:20])

    if(s=='HVG'){
      # Run ScMerge
      t1 = Sys.time()
      scMerge_res <- scMerge(
        sce_combine = sce,
        ctl = gene_lowvar,
        kmeansK = c(2,2),
        assay_name = "scMerge_res",
        replicate_prop = 0.5, verbose=T,
        marker = rownames(counts))
      t2 = Sys.time()
      print(t2-t1)
    } else{
    # Run ScMerge
      t1 = Sys.time()
      scMerge_res <- scMerge(
        sce_combine = sce, 
        ctl = gene_lowvar,
        kmeansK = c(2,2),
        assay_name = "scMerge_res",
        replicate_prop = 0.5, verbose=T,
        marker = rownames(counts_HVG))
      t2 = Sys.time()
      print(t2-t1)
    }
    
    # save the output
    save(scMerge_res, file = paste0('demo_scmerge/',x2,'/',s,"/output.rda"))
    scmerge_norm <- scMerge_res@assays@data$scMerge_res       
    write.table(scmerge_norm, file = paste0('demo_scmerge/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
    
    
  })
  
})


