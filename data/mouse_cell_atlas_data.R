library(scater)
rm(list=ls())

base_name='mouse_cell_atlas'
counts <- read.table("ds2/filtered_total_batch1_seqwell_batch2_10x.txt", head=T, sep='\t')
metadata <- read.table("ds2/filtered_total_sample_ext_organ_celltype_batch.txt", head=T, sep='\t')
rownames(metadata) <- gsub("-", "\\.", rownames(metadata))

K_intervals = c(2, 3, 4, 6, 7, 8)
cutoff = 0.95
a = b = 5
pp = 1

for (K in K_intervals)
{
  batch <- metadata$batch
  tcell = which(metadata$orig.ident=='T-cell')
  tcount= counts[, tcell]
  tbatch= batch[tcell]
  
  batch1 = tcount[, which(tbatch==1)]
  batch2 = tcount[, which(tbatch==2)]
  
  # Shuffling replicates    
  ori_cell_order <- colnames(batch1)
  sp_cells <- colnames(batch1)
  sp_cells <- sample(sp_cells, length(sp_cells), replace = FALSE)
  batch1 <- batch1[, sp_cells]
  colnames(batch1) <- ori_cell_order
  
  ori_cell_order <- colnames(batch2)
  sp_cells <- colnames(batch2)
  sp_cells <- sample(sp_cells, length(sp_cells), replace = FALSE)
  batch2 <- batch2[, sp_cells]
  colnames(batch2) <- ori_cell_order
  
  group <- tbatch
  for (i in unique(group)){
    gid = which(group==i)
    group[gid[1:(length(gid)%/%10*K)]] = paste0(i, "_A")
    group[gid[(length(gid)%/%10*K + 1):length(gid) ]] = paste0(i, "_B")
  }
  
  rowsum = rowSums(as.matrix(batch1)==0)/length(colnames(batch1))
  rowsum2 = rowSums(as.matrix(batch2)==0)/length(colnames(batch2))
  gb1 = rownames(batch1)[which(rowsum < cutoff)]
  print(length(gb1))
  gb2 = rownames(batch2)[which(rowsum2 < cutoff)]
  print(length(gb2))
  common_genes = intersect(gb1, gb2)
  print(length(common_genes))
  
  batch1 = batch1[common_genes, ]   
  batch2 = batch2[common_genes, ]   
  colnames(batch1)<-group[which(tbatch==1)]  
  colnames(batch2)<-group[which(tbatch==2)]
  
  N = round(length(common_genes) *pp /50)
  if (N%%2==1) N = N + 1
  print (N)
  DEG_groundthruth = sample(common_genes, N, replace=FALSE)
  second_N_genes  = sample(DEG_groundthruth, N%/%2, replace = FALSE)
  first_N_genes= DEG_groundthruth[!DEG_groundthruth %in% second_N_genes]
  
  #Down sampling
  batch1_ori = batch1
  batch2_ori = batch2
  
  tmp <- batch1[first_N_genes, which(colnames(batch1)==unique(colnames(batch1))[1])]
  for (g in rownames(tmp)){
    for (c in colnames(tmp)){
      prob = rbeta(1, a, b)
      tmp[g, c] = rbinom(1, tmp[g, c], prob=prob)
    }
  }
  batch1[first_N_genes, which(colnames(batch1)==unique(colnames(batch1))[1])] <- tmp
  
  tmp <- batch1[second_N_genes, which(colnames(batch1)==unique(colnames(batch1))[2])]
  for (g in rownames(tmp)){
    for (c in colnames(tmp)){
      prob = rbeta(1, a, b)
      tmp[g, c] = rbinom(1, tmp[g, c], prob=prob)
    }
  }
  batch1[second_N_genes, which(colnames(batch1)==unique(colnames(batch1))[2])] <- tmp
  
  ############################# batch 2
  tmp <- batch2[first_N_genes, which(colnames(batch2)==unique(colnames(batch2))[1])]
  for (g in rownames(tmp)){
    for (c in colnames(tmp)){
      prob = rbeta(1, a, b)
      tmp[g, c] = rbinom(1, tmp[g, c], prob=prob)
    }
  }
  batch2[first_N_genes, which(colnames(batch2)==unique(colnames(batch2))[1])] <- tmp
  
  tmp <- batch2[second_N_genes, which(colnames(batch2)==unique(colnames(batch2))[2])]
  for (g in rownames(tmp)){
    for (c in colnames(tmp)){
      prob = rbeta(1, a, b)
      tmp[g, c] = rbinom(1, tmp[g, c], prob=prob)
    }
  }
  batch2[second_N_genes, which(colnames(batch2)==unique(colnames(batch2))[2])] <- tmp
  
  ##################################################
  newmat <- cbind.DataFrame(batch1, batch2)
  
  if (anyNA(newmat)) {
    cat(paste0(K, "-NAs exist!"))
    exit
  }
  # colsum = colSums(as.matrix(newmat) == 0)
  # simsparsity = sum(colsum)/(length(rownames(newmat))*length(colnames(newmat)))
  # print(simsparsity)
  
  newgroup<-colnames(newmat)
  newgroup[newgroup=='1_A'] = 'A'
  newgroup[newgroup=='1_B'] = 'B'
  newgroup[newgroup=='2_A'] = 'A'
  newgroup[newgroup=='2_B'] = 'B'
  
  newlabel<-c('label',newgroup)
  
  newbatch<-tbatch
  geneinfo <- common_genes
  cellinfo <- as.data.frame(newgroup)
  cellinfo$batch <- factor(newbatch)
  
  colnames(cellinfo) <- c("Group", "Batch")
  up_genes <- second_N_genes
  down_genes <- first_N_genes
  de_genes <- DEG_groundthruth
  
  ds = paste0("mca_dgsp95_b155_",pp,"_", K)
  dir.create(ds)
  write.table(newmat, file = paste0(ds,"/counts.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(geneinfo, file = paste0(ds,"/geneinfo.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(cellinfo, file = paste0(ds,"/cellinfo.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(down_genes, file = paste0(ds,"/true_down_genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(up_genes, file = paste0(ds,"/true_up_genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(de_genes, file = paste0(ds,"/de_genes.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  
}  
