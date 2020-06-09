########## summary confusion matrix

rm(list=ls())

# load libraries
library(gridExtra)
library(openxlsx)
dir.create('Venn_diagram')

main_Fscore <- function(select){
  
  # METHODS
  vect_method <- c('raw_data','seurat3','MNN','Combat', 'limma',"limma_voom", 'zinb_wave',
                    'scmerge', 'edger', 'deseq2')
  
  vect_HVG <- c('all','HVG')
  
  # SIMULATIONS
  vect_simu <- gsub('data/','',list.dirs('data', recursive=FALSE))
  
  Fscore_list <- lapply(vect_HVG,function(HVG){
    
    base_name <- paste0('Venn_diagram/')
    
    df_all <- lapply(vect_simu,function(simu){
      
      df <- sapply(vect_method,function(method){
        
        if(method=='raw_data'){
          S3 <- read.table(paste0('Seurat_DEGs/',simu,'/','S3_batch12/',method,'_',HVG,'/degs_batch12__seurat_bimod_DEG.txt'), head=T, sep='\t')
        } else {
          if (method=='deseq2'){
            if (HVG=='HVG'){
              S3 <-read.table(paste0('data/',simu,'/','HVG_deseq2_result_table.txt'), head=T, sep='\t')
            } else {
              S3 <-read.table(paste0('data/',simu,'/','all_deseq2_result_table.txt'), head=T,  sep='\t')
            } 
            S3 <- data.frame(genes = row.names(S3), S3)
          } else if (method=='edger'){
            if (HVG=='HVG'){
              S3 <-read.table(paste0('data/',simu,'/','HVG_edger_result_table.txt'), head=T, sep='\t')
            } else {
              S3 <-read.table(paste0('data/',simu,'/','all_edger_result_table.txt'), head=T,  sep='\t')
            } 
            S3 <- data.frame(genes = row.names(S3), S3)
          } else if (method=='limma_voom'){
            if (HVG=='HVG'){
              S3 <-read.table(paste0('data/',simu,'/','HVG_voom_cov_result_table.txt'), head=T, sep='\t')
            } else {
              S3 <-read.table(paste0('data/',simu,'/','voom_cov_result_table.txt'), head=T,  sep='\t')
            } 
            S3 <- data.frame(genes = row.names(S3), S3)
          } 
          else {
            S3 <- read.table(paste0('Seurat_DEGs/',simu,'/','S3_batch12/','after_',method,'_',HVG,'/degs_batch12_seurat_bimod_DEG.txt'), head=T, sep='\t')
          }
        }
        geneinfo <- read.table(paste0('data/',simu,'/geneinfo.txt'), head=T)
        real_de_genes_ls <- rownames(geneinfo[(geneinfo$DEFacGroup1+geneinfo$DEFacGroup2)!=2,])
        S7 <- geneinfo[real_de_genes_ls,]
        
        if(HVG=='HVG'){
          # read HVG table 
          HVG <- read.table(paste0('data/',simu,'/counts_HVG.txt'), head=T, row.names=1)
          HVG <- colnames(HVG)
          S7 <- S7[rownames(S7) %in% HVG,] 
        } 
        
        if(select=='UP'){
          if (method=='deseq2'){
            S3 <- S3[S3$logFC<0,]
            S3 <- S3[S3$adjpvalue<0.05,]
          }
          else
          if (method=='edger'){
            S3 <- S3[S3$logFC<0,]
            S3 <- S3[S3$adjpvalue<0.05,]
          }
          else
          if (method=='limma_voom'){
            S3 <- S3[S3$logFC<0,]
            S3 <- S3[S3$adjpvalue<0.05,]
          } else {
            S3 <- S3[S3$avg_logFC>0,]
          }
          S7 <- S7[S7$DEFacGroup1>S7$DEFacGroup2,]
        } else if(select=='DOWN'){
          if (method=='deseq2'){
            S3 <- S3[S3$logFC>0,]
            S3 <- S3[S3$adjpvalue<0.05,]
          } else
          if (method=='edger'){
            S3 <- S3[S3$logFC>0,]
            S3 <- S3[S3$adjpvalue<0.05,]
          } else
          if (method=='limma_voom'){
            S3 <- S3[S3$logFC>0,]
            S3 <- S3[S3$adjpvalue<0.05,]
          } 
          else {
            S3 <- S3[S3$avg_logFC<0,]
          }
          S7 <- S7[S7$DEFacGroup1<S7$DEFacGroup2,]
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
        data <- data.frame(TP=TP,FN=FN,FP=FP,recall=TPR,precision=PPV,Fscore=Fscore,row.names = method)
        return(data)
        
      })
      
      mef <- t(df)
      rownames(mef) <- c('Raw','Seurat 3','MNN correct','Combat','limma','limma_voom', 'zinb_wave',
                         'scMerge', 'edger', 'deseq2')
     
      mef <- rbind(rep('',dim(mef)[2]),mef)
      rownames(mef)[rownames(mef)==""] <- simu
      return(list(mef,df['Fscore',]))
      
    })
    
    df_all_all <- do.call(rbind,lapply(df_all,function(l){l[[1]]}))
    
    # add Average block
    df_fscore <- do.call(rbind,lapply(df_all,function(l){unlist(l[[2]])}))
    
    write.table(df_all_all,paste0(base_name,'summary_confusion_matrix_',HVG,'_DEGs_',select,'.txt'),sep='\t', quote=F, row.names=T, col.names=NA)
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

plotdata$L1 <- factor(plotdata$L1,c(1,2),c('Up-regulated in Group 1','Down-regulated in Group 1'))

plotdata$Var2 <- relevel(plotdata$Var2, ref = "Combat")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "limma")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "MNN")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "scmerge")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "seurat3")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "zinb_wave")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "raw_data")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "deseq2")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "edger")
plotdata$Var2 <- relevel(plotdata$Var2, ref = "limma_voom")

meanplot = aggregate(plotdata[2],list(rep(1:(nrow(plotdata[2])%/%6+1),each=6,len=nrow(plotdata))),median)[-1]
rawmedian = rep(c(meanplot$F.score[1], meanplot$F.score[11], meanplot$F.score[21], meanplot$F.score[31]), 
                each=nrow(plotdata[2])%/%4)
plotdata$RawMedian = as.numeric(rawmedian)

# df2 = as.numeric(rawmedian)
# df2 <- as.data.frame(df2)
# colnames(df2)<-c('rawmedian')

library(plyr)
df <- ddply(plotdata,.(L2, L1),summarise,median=median(RawMedian, na.rm = TRUE))

p4 <- ggplot(plotdata,aes(x=Var2,y=F.score, color=Var2)) + geom_boxplot(outlier.size = 1) + coord_flip() +
  geom_hline(data=df, aes(yintercept=median),linetype="dashed", color='red') +  ylim(0.00, 1.00)+ labs(y = 'F-score')

L2.labs <- c('All genes', 'HVG')
names(L2.labs) <- c('all','HVG')

library(lemon)
p4 <- p4 + facet_rep_grid(L1~L2, labeller = labeller(L2 = L2.labs),repeat.tick.labels = 'x', scales='fixed')
p4 <- p4 + theme(axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                 axis.title.x = element_text(size=20),
                 axis.title.y = element_blank(),
                 #axis.text.x = element_text(size=11,colour = 'black',angle = 45, hjust = 1),
                 axis.text = element_text(size=17,colour = 'black'),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_blank(),
                 panel.spacing.x = unit(0.5, "lines"),
                 panel.spacing.y = unit(1, "lines"),
                 strip.text = element_text(size=20, color="black"),
                 strip.background.x = element_rect(fill="#CDE8DF"),
                 legend.position="none") 

p4 <- p4 + scale_x_discrete(breaks=c("raw_data","seurat3","MNN","Combat","Combat_seq","limma","limma_voom",
                                     "scmerge", 'edger', 'deseq2', "zinb_wave"),
                            labels=c("Raw","Seurat3","MNNCorrect","Combat", "Combat_seq", "limma-bec","limma", 
                                     "scMerge", 'edgeR', 'DESeq2', "ZINB-WaVE"))

p5 <- p4 + geom_point()#geom_jitter(shape=16, position=position_jitter(0.2),size=1)
p5

tiff(filename = "Venn_diagram/Fscore_boxplot.tiff", units="px", width=1200, height=900)
p5
dev.off()