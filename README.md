# On using batch-effect-corrected data for differential expressionanalysis of single-cell RNA sequencing data

### Group vs. Batch data distribution - Mouse Cell Atlas
 * We have plotted the distribution of 6 datasets we have simulated from Mouse Cell Atlas dataset. Each subplot clearly shows the large batch difference and the tiny difference between groups. Moreover, plots also show the various ratios of cells in each batch and group that we have implemented on purpose.
 <img src="data/group_vs_batch_distribution.png" width="900"> 

### Performance (F-score, Precision-recall curve) of selected methods on specific datasets
| **Splatter simulation dataset** | **Mouse Cell Atlas dataset** | **Human Pancreas dataset** |
| --- | --- | --- |
| <img src="data/splatter_80_gf_pca.png" width="400"> | <img src="data/mca_10pp_pca.png" width="400"> | <img src="data/pan_98_pca_o.png" width="400"> |
| <img src="data/splatter_80_gf_after.png" width="400"> | <img src="data/mca_gene_sparsity.png" width="400"> | <img src="data/pan_gene_sparsity.png" width="400"> |
| <img src="data/splatter_80gf_fscore.png" width="400"> | <img src="data/mca_10pp_fscore.png" width="400"> | <img src="data/pan_98_10pp_fscore.png" width="400"> |
| <img src="data/splatter_80gf_PR_curve.png" width="400"> | <img src="data/mca_10pp_PR_curve.png" width="400"> | <img src="data/pan_98_10pp_PR_curve.png" width="400"> |

### Comparison of sparsity level distribution over cells  
<img src="data/sparsity_distribution_over_cells.png" width="600"> 

### data 
  * Including code and the results of processing 2 real datasets: **Mouse Cell Atlas** and **Human Pancreas**

### code
  * This github includes **modified and additional code** for supporting the response to the article at: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9

  * The original code can be found at: https://github.com/JinmiaoChenLab/Batch-effect-removal-benchmarking
