# ---
# PCA
# ---
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
mat <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/matrix/scran.rds')
tmp <- PCA(genebycell_mat = mat, save.pca = TRUE, plot.statistics=TRUE, plot.dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/plot/', result.dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/result/pca/', PC_for_columns = TRUE, findVariableGenes = TRUE, maxVariableGenes = 3000, numPC = NULL, smoothFittingMethod = 'loess')

