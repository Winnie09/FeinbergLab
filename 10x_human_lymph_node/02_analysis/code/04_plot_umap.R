# ------------
# plot my umap
# ------------
plotdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/plot/'
umap <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/result/umap/umap.rds')
cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/matrix/scran.rds')
gene = c('CD3D', 'CD19')
pd <- data.frame(umap1 = umap[,1], umap2 = umap[,2], g1 = cnt[gene[1],], g2 = cnt[gene[2],], stringsAsFactors = FALSE)

library(ggplot2)
library(scattermore)
library(viridis)
library(gridExtra)
p1 <- ggplot(pd, aes(x = umap1, y = umap2, color = g1)) + 
  geom_scattermore(pointsize = 1.5)+
  theme_classic() +
  scale_colour_viridis() +
  # scale_color_gradientn(colors=rev(rainbow(15)[-c(12:15)])) + 
  ggtitle(gene[1])
p2 <- ggplot(pd, aes(x = umap1, y = umap2, color = g2)) + 
  geom_scattermore(pointsize = 1.5)+
  theme_classic()+
  scale_colour_viridis() +
  # scale_color_gradientn(colors=rev(rainbow(15)[-c(12:15)])) + 
  ggtitle(gene[2])
pdf(paste0(plotdir, 'umap_my.pdf'), width = 7, height = 3)
grid.arrange(p1,p2, nrow =1)
dev.off()


