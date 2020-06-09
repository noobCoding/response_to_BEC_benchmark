
rm(list=ls())
# BiocManager::install('Seurat')
library(Seurat)  # Seurat v3 version
library(cowplot)
library(ggplot2)
# 
lsdir <- list.dirs('data', recursive=FALSE) 
dir.create('demo_seurat3')
sapply(lsdir,function(x){
  
  x2 <- gsub('data/','',x)

  dir.create(paste0('demo_seurat3/',x2), showWarnings = FALSE)
  
  selection <- c('all', 'HVG')
  sapply(selection, function(s){
    
    dir.create(paste0('demo_seurat3/',x2,'/',s), showWarnings = FALSE)
    
    # read data counts and cellinfo
    if(s=='HVG'){
      counts <- read.table(paste0(x,'/counts_HVG.txt'), head=T, sep='\t')
    } else {
      counts <- read.table(paste0(x,'/counts.txt'), head=T, sep='\t')
    }
    counts <- t(counts)
    cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t')
    cellinfo <- cellinfo[colnames(counts),]
    
    pbmc <- CreateSeuratObject(counts = counts, project = '', min.cells = 0, min.features = 0, meta.data = subset(cellinfo,select=c('Batch','Group')))
    
    pbmc.list <- SplitObject(pbmc, split.by = "Batch")
    
    # Run Seurat V3 integration
    t1 = Sys.time()
    for (i in names(pbmc.list)) {
      pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = FALSE)
    }
    pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = pbmc@assays$RNA@counts@Dim[1])
    pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features)
    pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", anchor.features = pbmc.features, dims = 1:70)
    immune.combined <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT", dims = 1:70)
    t2 = Sys.time()
    print(t2-t1)
    
    # save the output
    save(immune.combined, file = paste0('demo_seurat3/',x2,'/',s,"/output.rda"))
    seuratv3_integrated <- immune.combined@assays$integrated@data
    write.table(seuratv3_integrated, file = paste0('demo_seurat3/',x2,'/',s,"/output.txt"), row.names = T, col.names = T, sep="\t")
    
  })
  
})


