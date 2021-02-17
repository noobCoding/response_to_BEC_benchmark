
rm(list=ls())
# BiocManager::install("scran")
# BiocManager::install("scales")
# BiocManager::install("Rtsne")
# BiocManager::install(version='devel')
# BiocManager::install("batchelor")
library(scran)
library(scales)
require(Rtsne)
library(Seurat)
library(ggplot2)
library(cowplot)
library(batchelor)
# 

dir.create('demo_MNN')
lsdir <- list.dirs('data', recursive=FALSE) 

sapply(lsdir,function(x){
  
  x2 <- gsub('data/','',x)
  dir.create(paste0('demo_MNN/',x2), showWarnings = FALSE)
  selection <- c('HVG','all')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_MNN/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    if(s=='HVG'){
      counts <- read.table(paste0(x,'/counts_HVG.txt'), head=T, sep='\t', fill=T)
      counts<-t(counts)
      rownames(counts) <- gsub('.', '-', rownames(counts), fixed = TRUE)
    } else {
      counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t', fill=T)
    }
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t', fill=T)
    rownames(cellinfo) <- factor(colnames(counts))
    
    pbmc <- CreateSeuratObject(counts = counts, meta.data = cellinfo)
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")
    
    pbmc.list <- SplitObject(pbmc, split.by = "Batch")
    myData1 <- pbmc.list[[1]]@assays$RNA@data
    myData2 <- pbmc.list[[2]]@assays$RNA@data
    
    # Run MNN
    t1 = Sys.time()
    out.mnn.total <- batchelor::mnnCorrect(myData1, myData2, k=20, sigma=1, cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE)
    t2 = Sys.time()
    print(t2-t1)
    
    # save the output
    save(out.mnn.total, file = paste0('demo_MNN/',x2,'/',s,"/output.rda"))
    corre.mnn <- out.mnn.total@assays@data$corrected # @assays[['corrected']]
    write.table(corre.mnn, file = paste0('demo_MNN/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
   
  })
  
})


