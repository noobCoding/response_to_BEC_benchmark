
# library(Seurat)
# install.packages('remotes')
# remotes::install_version("SDMTools", "1.1-221")
# remotes::install_version("multtest")
# BiocManager::install('devtools')
library(ggplot2)
# install.packages('Seurat', lib = '~/R/x86_64-pc-linux-gnu-library/4.0/Seurat3' )
# devtools::install_version(package = 'Seurat', version = package_version('2.3.0'), )
# library(Seurat, lib.loc = '/home/admin/R/x86_64-pc-linux-gnu-library/4.0/Seurat2/libs/')
library(Seurat)

lsdir <- list.dirs('data', recursive=FALSE)

sapply(lsdir,function(x){
  # read data counts and cellinfo
  counts <- read.table(paste0(x,'/counts.txt'), header = TRUE, sep='\t', fill = TRUE)
  # counts <- t(counts)
  cellinfo <- read.table(paste0(x,'/cellinfo.txt'), head=T, sep='\t', fill = TRUE)
  # cellinfo <- cellinfo[colnames(counts),]
  
  pbmc <- CreateSeuratObject(raw.data = counts, min.cells = 3)
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
  
  png(paste0(x,'/variableGenes.png'),height = 480, width = 480, res = 72)
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  dev.off()
  length(x = pbmc@var.genes)
  
  counts_HVG <- counts[rownames(counts) %in% pbmc@var.genes,]
  write.table(t(counts_HVG), file = paste0(x,'/counts_HVG.txt'), quote=FALSE, row.names = T, col.names = T, sep="\t")
})
