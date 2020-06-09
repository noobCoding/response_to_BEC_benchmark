
library(Seurat)
library(sva)

dir.create('demo_Combat_seq')
lsdir <- list.dirs('data', recursive=FALSE) 

sapply(lsdir,function(x){

  x2 <- gsub('data/','',x)

  dir.create(paste0('demo_Combat_seq/',x2), showWarnings = FALSE)
  
  selection <- c('HVG', 'all')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_Combat_seq/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    if(s=='HVG'){
      counts <- read.table(paste0(x,'/counts_HVG.txt'), head=T, sep='\t')
    } else {
      counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t')
    }
    counts <- t(counts)
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t')
    cellinfo <- cellinfo[colnames(counts),]
    
    # Normalize each library to the median of the transcript counts across all cells 
    # Then, log transform expression values   
    print("Median normalizing counts and log-transforming")
    col_sums = apply(counts,2, sum)
    med_trans = median(col_sums)

    # remove genes with variance equals to 0
    rv_genes <- which(apply(counts, 1, var)==0) # apply on normalized data
    rv_genes_names <- rownames(counts)[rv_genes]
    count_df <- counts[!(rownames(counts) %in% rv_genes_names),]
    
    # Run COMBAT_seq
    t1 = Sys.time()
    
    cellGroup <- factor(cellinfo$Group)
    cellBatch <- factor(cellinfo$Batch)
    combat_seq <- ComBat_seq(count_df, batch=cellBatch, group=cellGroup, full_mod = TRUE)
    
    t2 = Sys.time()
    
    # save the output
    save(combat_seq,file=paste0('demo_Combat_seq/',x2,'/',s,"/output.rda"))
    write.table(combat_seq, file = paste0('demo_Combat_seq/',x2,'/',s,"/output.txt"), quote=FALSE, row.names = T, col.names = T, sep="\t")
    
  })
})
