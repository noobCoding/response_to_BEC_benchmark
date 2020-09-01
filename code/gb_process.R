#library('Seurat') 

library('rstudioapi')
jobRunScript("GBC-run_limma.R")
jobRunScript("GBC-run_combat.R")
jobRunScript("GBC-run_edger.R")
jobRunScript("GBC-run_deseq2.R")
jobRunScript("GBC-run_MNN.R")
jobRunScript("GBC-run_limma_voom.R")
jobRunScript("GBC-run_scmerge.R")
jobRunScript("GBC-run_seurat3.R")

jobRunScript("GBC-run_DEGs_from_Seurat.R")
