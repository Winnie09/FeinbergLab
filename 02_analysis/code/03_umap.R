# ---- 
# umap
# ----
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
mat <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/result/pca/pr.rds')
UMAP(samplebyfeature_mat = mat, save.umap = TRUE, result.dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/result')

