########## summary confusion matrix

rm(list=ls())

# load libraries
library(gridExtra)
library(openxlsx)
dir.create('Venn_diagram')

main_Fscore <- function(select){
  # METHODS
  vect_method <- c('raw_data','MNN','Combat', 'limma',"limma_voom",
                   'edger', 'deseq2', 'seurat3' , 'scmerge', 'zinbwave'
                   )
  # HVG (all/as Seurat)
  vect_HVG <- c('HVG', 'all')
  
  # SIMULATIONS
  vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
  
  Fscore_list <- lapply(vect_HVG,function(HVG){
    base_name <- paste0('Venn_diagram/')
    
    df_all <- lapply(vect_simu,function(simu){
      
      df <- sapply(vect_method,function(method){
        
        if(method=='raw_data'){
          S3 <- read.table(paste0('Seurat_DEGs/',simu,'/','S3_batch12/',method,'_',HVG,'/degs_batch12__seurat_bimod_DEG.txt'), head=T, sep='\t', fill=T)
          S3 <- data.frame(genes = rownames(S3), S3)
        } else {
          if (method=='deseq2'){
            if (HVG=='HVG'){
              S3 <-read.table(paste0('data/',simu,'/','HVG_deseq2_result_table.txt'), head=T, sep='\t', fill=T)
            } else {
              S3 <-read.table(paste0('data/',simu,'/','all_deseq2_result_table.txt'), head=T,  sep='\t', fill=T)
            } 
            S3 <- data.frame(genes = row.names(S3), S3)
          } else if (method=='edger'){
            if (HVG=='HVG'){
              S3 <-read.table(paste0('data/',simu,'/','HVG_edger_result_table.txt'), head=T, sep='\t', fill=T)
            } else {
              S3 <-read.table(paste0('data/',simu,'/','all_edger_result_table.txt'), head=T,  sep='\t', fill=T)
            } 
            S3 <- data.frame(genes = row.names(S3), S3)
          } else if (method=='limma_voom'){
            if (HVG=='HVG'){
              S3 <-read.table(paste0('data/',simu,'/','HVG_voom_result_table.txt'), head=T, sep='\t', fill=T)
            } else {
              S3 <-read.table(paste0('data/',simu,'/','voom_result_table.txt'), head=T,  sep='\t', fill=T)
            } 
            S3 <- data.frame(genes = row.names(S3), S3)
          } 
          else {
            S3 <- read.table(paste0('Seurat_DEGs/',simu,'/','S3_batch12/','after_',method,'_',HVG,'/degs_batch12_seurat_bimod_DEG.txt'), head=T, sep='\t', fill=T)
          }
        }
        geneinfo <- read.table(paste0('data/',simu,'/geneinfo.txt'), head=T, fill=T)
        real_de_genes_ls <- read.table(paste0('data/',simu,'/de_genes.txt'), head=T, fill=T)
        S7 <- as.data.frame(geneinfo[which(geneinfo$x %in% real_de_genes_ls$x),])
        rownames(S7) <- S7[[1]]
        colnames(S7) <- c('Gene')
        
        cutoff = 0.05
        if(select=='UP'){
          if (method=='deseq2'){
            S3 <- S3[S3$logFC<0,]
            S3 <- S3[S3$adjpvalue<cutoff,]
          }
          else
            if (method=='edger'){
              S3 <- S3[S3$logFC<0,]
              S3 <- S3[S3$adjpvalue<cutoff,]
            }
          else
            if (method=='limma_voom'){
              S3 <- S3[S3$logFC<0,]
              S3 <- S3[S3$adjpvalue<cutoff,]
            } else {
              S3 <- S3[S3$avg_logFC>0,]
            }
          
          up_genes <- read.table(paste0('data/',simu,'/true_up_genes.txt'), head=T, fill = T)
          S7 <- as.data.frame(geneinfo[which(geneinfo$x %in% up_genes$x),])
          rownames(S7) <- S7[[1]]
          colnames(S7) <- c('Gene')
          if(HVG=='HVG')
          {
            hvgenes <- read.table(paste0('data/',simu,'/counts_HVG.txt'), head=T, row.names=NULL, fill = T)
            hvgenes <- colnames(hvgenes)
            hvgenes <- gsub('.', '-', hvgenes, fixed = TRUE)
            S7 <- as.data.frame(S7[rownames(S7) %in% hvgenes,])
            colnames(S7) <- c('Gene')
          } 
          
        } else if(select=='DOWN'){
          if (method=='deseq2'){
            S3 <- S3[S3$logFC>0,]
            S3 <- S3[S3$adjpvalue<cutoff,]
          } else
            if (method=='edger'){
              S3 <- S3[S3$logFC>0,]
              S3 <- S3[S3$adjpvalue<cutoff,]
            } else
              if (method=='limma_voom'){
                S3 <- S3[S3$logFC>0,]
                S3 <- S3[S3$adjpvalue<cutoff,]
              } 
          else {
            S3 <- S3[S3$avg_logFC<0,]
          }
          # S7 <- S7[S7$DEFacGroup1<S7$DEFacGroup2,]
          down_genes <- read.table(paste0('data/',simu,'/true_down_genes.txt'), head=T, fill = T)
          S7 <- as.data.frame(geneinfo[which(geneinfo$x %in% down_genes$x),])
          rownames(S7) <- S7[[1]]
          colnames(S7) <- c('Gene')
          if(HVG=='HVG'){
            # read HVG table 
            hvgenes <- read.table(paste0('data/',simu,'/counts_HVG.txt'), head=T, row.names=NULL, fill = T)
            hvgenes <- colnames(hvgenes)
            hvgenes <- gsub('.', '-', hvgenes, fixed = TRUE)
            S7 <- as.data.frame(S7[rownames(S7) %in% hvgenes,])
            colnames(S7) <- c('Gene')
          } 
        } else {
          stop('select UP or DOWN DEGs')
        }
        
        GT <- S7$Gene
        if (method=='deseq2'){
          norm <- S3$genes
        } else
          if (method=='edger'){
            norm <- S3$genes
          } else if (method=='limma_voom'){
            norm <- S3$genes
          }
        else {
          norm <- S3$X
        }
        
        # table
        TP <- sum(GT %in% norm)
        FN <- length(GT)-TP
        FP <- length(norm)-TP
        
        # recall (sensitivity)
        TPR <- round(TP/(TP + FN),3)
        
        # precision (positive predictive value)
        PPV <- round(TP/(TP + FP),3)
        
        # F-score
        Fscore <- 2*((PPV*TPR)/(PPV+TPR))
        Fscore <- round(Fscore,3)
        
        # table matrix confusion 
        data <- data.frame(TP=TP,FN=FN,FP=FP,recall=TPR,precision=PPV,Fscore=Fscore, row.names = method)
        return(data)
      })
      
      mef <- t(df)
      rownames(mef) <- c('Raw','MNN', 'Combat','limma', 'Seurat 3', 'scmerge', 'zinbwave',
                         'limma_voom' , 'edger', 'deseq2')
      
      mef <- rbind(rep('',dim(mef)[2]),mef)
      rownames(mef)[rownames(mef)==""] <- simu
      # print (mef)
      return(list(mef,df['Fscore',]))
      
    })
    
    df_all_all <- do.call(rbind,lapply(df_all,function(l){l[[1]]}))
    
    # add Average block
    df_fscore <- do.call(rbind,lapply(df_all,function(l){unlist(l[[2]])}))
    
    write.table(df_all_all,paste0(base_name,'summary_confusion_matrix_',HVG,'_DEGs_',select,'.txt'),sep='\t', quote=F)
    return(df_fscore)
    
  })
  names(Fscore_list) <- vect_HVG
  return(Fscore_list)
}


######################### BOXPLOT OF FSCORE
library(ggplot2)
library(cowplot)

# Prepare the data
Fscore_up <- main_Fscore(select='UP') # comment write.table(df_all_all)
Fscore_down <- main_Fscore(select='DOWN') # comment write.table(df_all_all)
plotdata <- reshape2::melt(list(Fscore_up,Fscore_down), value.name = "F.score")
plotdata <- plotdata[,-1]
plotdata$F.score[is.nan(plotdata$F.score)] = 0.

plotdata$L1 <- factor(plotdata$L1,c(1,2),c('Up-regulated in Group 1','Down-regulated in Group 1'))

plotdata$Var2 <- relevel(plotdata$Var2, ref = "Combat")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "limma")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "MNN")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "scmerge")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "seurat3")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "zinbwave")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "raw_data")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "deseq2")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "edger")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "limma_voom")

meanplot = aggregate(plotdata[2],list(rep(1:(nrow(plotdata[2])%/%6+1),each=6,len=nrow(plotdata))),median)[-1]

t = length(unique(plotdata$Var2))
rawmedian = rep(c(meanplot$F.score[1], meanplot$F.score[1 + t], 
                  meanplot$F.score[1 + 2*t], meanplot$F.score[1 + 3*t]), 
                each=nrow(plotdata[2])%/%4)
plotdata$RawMedian = as.numeric(rawmedian)

# df2 = as.numeric(rawmedian)
# df2 <- as.data.frame(df2)
# colnames(df2)<-c('rawmedian')

library(plyr)
df <- ddply(plotdata,.(L2, L1),summarise,median=median(RawMedian, na.rm = TRUE))
# df2 <- ddply(plotdata,.(RawMedian),summarise,median=median(RawMedian, na.rm = TRUE))
# df$median = df2$median

p4 <- ggplot(plotdata,aes(x=Var2,y=F.score, color=Var2)) + geom_boxplot(outlier.size = 1) + coord_flip() +
  geom_hline(data=df, aes(yintercept=median),linetype="dashed", color='red') + 
  labs(y = 'F-score') #+  ylim(0.00, 0.50)

  # geom_hline(plotdata$RawMedian, linetype="dashed", color = "red") + labs(y = 'F-score') +

L2.labs <- c('All genes', 'HVG')
names(L2.labs) <- c('all','HVG')
# p4 <- p4 + facet_grid(L2~L1, labeller = labeller(L2 = L2.labs))
library(lemon)
p4 <- p4 + facet_rep_grid(L1~L2, labeller = labeller(L2 = L2.labs),repeat.tick.labels = 'x', scales='fixed')
p4 <- p4 + theme(axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                 axis.title.x = element_text(size=16),
                 axis.title.y = element_blank(),
                 #axis.text.x = element_text(size=11,colour = 'black',angle = 45, hjust = 1),
                 axis.text = element_text(size=14,colour = 'black'),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_blank(),
                 panel.spacing.x = unit(0.5, "lines"),
                 panel.spacing.y = unit(1, "lines"),
                 strip.text = element_text(size=17, color="black"),
                 strip.background.x = element_rect(fill="#CDE8DF"),
                 # strip.background.x = element_blank(),
                 # strip.background.y = element_blank(),
                 legend.position="none") 


# p4 <- p4 + scale_x_discrete(breaks=c("raw_data","zinb_wave","seurat3","scmerge","scGen",'scanorama2',"MNN","limma","Combat"),
#                             labels=c("Raw","ZINB-WaVE","Seurat 3","scMerge","scGen",'Scanorama',"MNN Correct","limma","Combat"))
# p4 <- p4 + scale_x_discrete(breaks=c("raw_data","seurat3","limma"),
#                             labels=c("Raw","Seurat 3","limma"))
p4 <- p4 + scale_x_discrete(breaks=c("raw_data","seurat3","MNN","Combat","Combat_seq","limma","limma_voom",
                                     "scmerge", 'edger', 'deseq2', "zinbwave"),
                                     # "scmerge", 'edger', 'deseq2'),
                            labels=c("Raw","Seurat3","MNNCorrect","Combat", "Combat_seq", "limma-bec","limma", 
                                     "scMerge", 'edgeR', 'DESeq2', "ZINB-WaVE"))
                                     # "scMerge", 'edgeR', 'DESeq2'))
 
# #save
# # ggsave(filename=paste0('Venn_diagram/Fscore_boxplot_flip.png'), plot=p4, width = 18, height = 16)
# tiff("Venn_diagram/Fscore_boxplot_flip.tiff", units="px", width=1000, height=750, res=300)
p4
dev.off()
# 
p5 <- p4 + geom_point()#geom_jitter(shape=16, position=position_jitter(0.2),size=1)
p5
# ggsave(filename=paste0('Venn_diagram/Fscore_boxplot_flip2.png'), plot=p4, width = 18, height = 16)
tiff(filename = "Venn_diagram/Fscore_boxplot.tiff", units="px", width=900, height=800)
p5
dev.off()
