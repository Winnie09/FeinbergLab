library(Seurat)
setwd('/home/whou10/scratch4/whou10/FeinbergLab/mousePancreas/visium/deconvolute/RCTD')
load('adm_spatial_RCTD_manual_and_pseudospot_normalized_weights_matrix.rda')
load('/home/whou10/scratch4/whou10/FeinbergLab/mousePancreas/visium/data/seurat_adm_spatial_for_ji_lab.rda')

adm.combined@meta.data <- data.frame(adm.combined@meta.data,norm_weights[match(rownames(adm.combined@meta.data),rownames(norm_weights)),])
for (ct in colnames(norm_weights)) {
  pdf(paste0('plot/',ct,'.pdf'),width=12,height=12)
  print(SpatialFeaturePlot(adm.combined,features=ct,ncol=5))
  dev.off()
}

