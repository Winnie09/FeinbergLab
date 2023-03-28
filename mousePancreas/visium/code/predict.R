library(Seurat)

source('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/software/predict.R')
mod <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/spatial/model/model.rds')

af <- list.files('/work-zfs/hji7/whou10/feinberg/pancreasVisium/visium_slides/',pattern = '-outs')

gl <- c('Klf4','Ep300','Kras','Cdkn1a','Notch1')

for (f in af) {
  print(f)
  d <- Load10X_Spatial(paste0('/work-zfs/hji7/whou10/feinberg/pancreasVisium/visium_slides/',f))
  d <- NormalizeData(d)
  
  pred <- predict(d@assays$Spatial@data,mod)
  
  for (gene in gl) {
    pdf(paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/spatial/plot/',sub('-outs','',f),'_rna.pdf'),width=15,height=3)
    print(SpatialFeaturePlot(d, features = gl,ncol=length(gl)))
    dev.off()

    d@meta.data <- data.frame(d@meta.data,t(pred[gl,rownames(d@meta.data)]))
    colnames(d@meta.data)[4:ncol(d@meta.data)] <- paste0(colnames(d@meta.data)[4:ncol(d@meta.data)],'_pred')
    pdf(paste0('/home-4/whou10@jhu.edu/work-zfs/whou10/spatial/plot/',sub('-outs','',f),'_metpred.pdf'),width=15,height=3)
    print(SpatialFeaturePlot(d, features = paste0(gl,'_pred'),ncol=length(gl)))
    dev.off()    
  }
}

