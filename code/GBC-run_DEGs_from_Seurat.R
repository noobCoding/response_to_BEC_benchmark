
rm(list=ls())
# this.dir <- '/home/marion/Marion/Project/Hoa_batch_normalization/simulation_dataset_V3/'
# setwd(this.dir)

source('GBC-Seurat_DEG_analysis.R')
# source('GBC-Seurat_DEG_analysis_pval.R')

######## S3 batch12 after normalization 

# METHODS
# vect_method <- c('Combat','limma','MNN')
vect_method <- c('Combat', 'limma', 'MNN', 'scmerge', 'seurat3')
# vect_method <- c('scmerge')

# HVG (all/as Seurat)
vect_HVG <- c('all','HVG')
# vect_HVG <- c('all')
dir.create('Seurat_DEGs')

# SIMULATIONS
vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))

sapply(vect_simu,function(simu){
  sapply(vect_method,function(method){
    sapply(vect_HVG,function(HVG){

      base_name <- paste0('Seurat_DEGs/',simu,'/')
      dir.create(base_name, showWarnings = FALSE)
      base_name <- paste0(base_name,'S3_batch12/')
      dir.create(base_name, showWarnings = FALSE)
      base_name <- paste0(base_name,'after_',method,'_',HVG,'/')
      dir.create(base_name, showWarnings = FALSE)

      # import data, sample
      inside <- list.files(paste0('demo_',method,'/',simu,'/',HVG), recursive=FALSE)
      if(grepl('.txt',inside[grep('output\\.((txt)|(csv))',inside)])){
        output_batch12 <- read.table(file = paste0('demo_',method,'/',simu,'/',HVG,'/output.txt'),sep="\t",header=T,row.names=1,check.names = F)
      } else {
        output_batch12 <- read.csv(file = paste0('demo_',method,'/',simu,'/',HVG,'/output.csv'),sep=",",header=T,row.names=1,check.names = F)
        output_batch12 <- t(output_batch12)
      }

      sample_batch12 <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1)
      rownames(sample_batch12) = colnames(output_batch12)
      # dim(sample_batch12)

      #### DEGs Group1 vs Group2 (batch 1 + batch 2)
      seurat_analysis_deg(TPM=output_batch12,
                          sample=sample_batch12,
                          group_col="Group",
                          base_name=paste0(base_name,'degs_batch12'),
                          group1='A',
                          group2='B',
                          test.use="wilcox",
                          logfc.threshold = 0)

      # #### DEGs Group1 vs Group2 (batch 1)
      # seurat_analysis_deg(TPM=output_batch12[,colnames(output_batch12) %in% rownames(sample_batch12[sample_batch12$Batch=='Batch1',])],
      #                     sample=sample_batch12[sample_batch12$Batch=='Batch1',],
      #                     group_col="Group",
      #                     base_name=paste0(base_name,'degs_batch1'),
      #                     group1='Group1',
      #                     group2='Group2',
      #                     test.use="bimod",
      #                     logfc.threshold = 0)
      #
      # #### DEGs Group1 vs Group2 (batch 2)
      # seurat_analysis_deg(TPM=output_batch12[,colnames(output_batch12) %in% rownames(sample_batch12[sample_batch12$Batch=='Batch2',])],
      #                     sample=sample_batch12[sample_batch12$Batch=='Batch2',],
      #                     group_col="Group",
      #                     base_name=paste0(base_name,'degs_batch2'),
      #                     group1='Group1',
      #                     group2='Group2',
      #                     test.use="bimod",
      #                     logfc.threshold = 0)

    })
  })
})


###### S3 batch1+batch2 without correction

# HVG (all/as Seurat)
vect_HVG <- c('all','HVG')
# vect_HVG <- c('HVG')

# SIMULATIONS
vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))

sapply(vect_simu,function(simu){
  sapply(vect_HVG,function(HVG){

    dir.create('Seurat_DEGs')
    base_name <- paste0('Seurat_DEGs/',simu,'/')
    dir.create(base_name, showWarnings = FALSE)
    base_name <- paste0(base_name,'S3_batch12/')
    dir.create(base_name, showWarnings = FALSE)
    base_name <- paste0(base_name,'raw_data_',HVG,'/')
    dir.create(base_name, showWarnings = FALSE)

    sample <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1)

    # import data, sample
    if(HVG=='HVG'){
      data <- read.table(file = paste0('data/',simu,'/counts_HVG.txt'),sep="\t",header=T,row.names=1,check.names = F)
      data <- t(data)
      rownames(data) <- gsub('.', '-', rownames(data), fixed = TRUE)
    } else {
      data <- read.table(file = paste0('data/',simu,'/counts.txt'),sep="\t",header=T,row.names=1,check.names = F)
    }
    rownames(sample) = colnames(data)

    #### DEGs Group1 vs Group2 (batches 1+2)
    seurat_analysis_deg2(TPM=data,
                         sample=sample,
                         group_col="Group",
                         base_name=paste0(base_name,'degs_batch12_'),
                         group1='A',
                         group2='B',
                         test.use="wilcox",
                         logfc.threshold = 0)
  })
})

# 
# ######## S1 batch1 without correction
# 
# # HVG (all/as Seurat)
# vect_HVG <- c('all','HVG')
# 
# # SIMULATIONS
# vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
# 
# sapply(vect_simu,function(simu){
#     sapply(vect_HVG,function(HVG){
# 
#       base_name <- paste0('Seurat_DEGs/',simu,'/')
#       dir.create(base_name, showWarnings = FALSE)
#       base_name <- paste0(base_name,'S1_batch1/')
#       dir.create(base_name, showWarnings = FALSE)
#       base_name <- paste0(base_name,'raw_data_',HVG,'/')
#       dir.create(base_name, showWarnings = FALSE)
#       
#       sample <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1)
# 
#       # import data, sample
#       if(HVG=='HVG'){
#         data <- read.table(file = paste0('data/',simu,'/counts.txt'),sep="\t",header=T,row.names=1,check.names = F)
#       } else {
#         data <- read.table(file = paste0('data/',simu,'/counts_HVG.txt'),sep="\t",header=T,row.names=1,check.names = F)
#       }
#       data <- t(data)
# 
#       #### DEGs Group1 vs Group2 (batch 1)
#       seurat_analysis_deg2(TPM=data[,colnames(data) %in% rownames(sample[sample$Batch=='Batch1',])],
#                           sample=sample[sample$Batch=='Batch1',],
#                           group_col="Group",
#                           base_name=paste0(base_name,'degs_batch1_'),
#                           group1='Group1',
#                           group2='Group2',
#                           test.use="bimod",
#                           logfc.threshold = 0)
#     })
# })
# 
# 
# ######## S2 batch2 without correction
# 
# # HVG (all/as Seurat)
# vect_HVG <- c('all','HVG')
# 
# # SIMULATIONS
# vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
# 
# sapply(vect_simu,function(simu){
#   sapply(vect_HVG,function(HVG){
# 
#     base_name <- paste0('Seurat_DEGs/',simu,'/')
#     dir.create(base_name, showWarnings = FALSE)
#     base_name <- paste0(base_name,'S2_batch2/')
#     dir.create(base_name, showWarnings = FALSE)
#     base_name <- paste0(base_name,'raw_data_',HVG,'/')
#     dir.create(base_name, showWarnings = FALSE)
#     
#     sample <- read.table(paste0("data/",simu,"/cellinfo.txt"),sep="\t",header=T,row.names=1)
# 
#     # import data, sample
#     if(HVG=='HVG'){
#       data <- read.table(file = paste0('data/',simu,'/counts.txt'),sep="\t",header=T,row.names=1,check.names = F)
#     } else {
#       data <- read.table(file = paste0('data/',simu,'/counts_HVG.txt'),sep="\t",header=T,row.names=1,check.names = F)
#     }
#     data <- t(data)
# 
#     #### DEGs Group1 vs Group2 (batch 2)
#     seurat_analysis_deg2(TPM=data[,colnames(data) %in% rownames(sample[sample$Batch=='Batch2',])],
#                         sample=sample[sample$Batch=='Batch2',],
#                         group_col="Group",
#                         base_name=paste0(base_name,'degs_batch2_'),
#                         group1='Group1',
#                         group2='Group2',
#                         test.use="bimod",
#                         logfc.threshold = 0)
#   })
# })

