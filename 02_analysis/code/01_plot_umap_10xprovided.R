# -------------------------
# plot umap provided by 10x
# -------------------------
plotdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/plot/'
umap.tb <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/spatial/analysis/umap/2_components/projection.csv', as.is = TRUE)
cnt <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/matrix/rawcount.rds')
colnames(cnt) <- sub('.*_', '', colnames(cnt))
cnt <- cnt[, umap.tb[,1]]
gene = c('CD3D', 'CD19')
pd <- data.frame(umap1 = umap.tb[,2], umap2 = umap.tb[,3], g1 = cnt[gene[1],], g2 = cnt[gene[2],], stringsAsFactors = FALSE)

library(ggplot2)
library(scattermore)
library(viridis)
library(gridExtra)
p1 <- ggplot(pd, aes(x = umap1, y = umap2, color = g1)) + 
  geom_scattermore(pointsize = 1.5)+
  theme_classic() +
  scale_colour_viridis() +
  ggtitle(gene[1])
p2 <- ggplot(pd, aes(x = umap1, y = umap2, color = g2)) + 
  geom_scattermore(pointsize = 1.5)+
  theme_classic()+
  scale_colour_viridis() +
  ggtitle(gene[2])
pdf(paste0(plotdir, 'umap_10xprovided.pdf'), width = 7, height = 3)
grid.arrange(p1,p2, nrow =1)
dev.off()

