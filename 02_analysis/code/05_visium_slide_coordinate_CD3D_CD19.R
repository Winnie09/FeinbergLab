# -------------
# spatial spots
# -------------
loc <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/spatial/spatial/tissue_positions_list.csv', as.is = TRUE)

colnames(loc) <- c('barcode', 'in_tissue', 'array_row','array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres')
rownames(loc) <- loc[,1]
loc <- loc[,-1]
dnam <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/03_pred_DNAm/result/predicted_DNAm.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/saver/saver.rds')
colnames(expr) <- sub('.*_', '', colnames(expr))
colnames(pred) <- sub('.*_', '', colnames(pred))
gene = c('CD3D', 'CD19')
pd = data.frame(x = loc[colnames(expr),4], y = loc[colnames(expr),3], g1 = expr[gene[1],], g2 = expr[gene[2],], stringsAsFactors = FALSE)

library(ggplot2)
library(scattermore)
library(viridis)
library(gridExtra)
p1 <- ggplot(pd, aes(x = x, y = y, color = g1)) + 
  geom_scattermore(pointsize = 3)+
  theme_classic() +
  scale_colour_viridis() +
  # scale_color_gradientn(colors=rev(rainbow(15)[-c(12:15)])) + 
  ggtitle(gene[1]) +
  xlab('10x Visium Slide Coordinate 1') +
  ylab('10x Visium Slide Coordinate 2')
p2 <- ggplot(pd, aes(x = x, y = y, color = g2)) + 
  geom_scattermore(pointsize = 3)+
  theme_classic()+
  scale_colour_viridis() +
  # scale_color_gradientn(colors=rev(rainbow(15)[-c(12:15)])) + 
  ggtitle(gene[2])+
  xlab('10x Visium Slide Coordinate 1') +
  ylab('10x Visium Slide Coordinate 2')
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/plot/visium_slide_coordinate_CD3D_CD19.pdf', width = 7, height = 3)
grid.arrange(p1,p2, nrow =1)
dev.off()
