# ----------------------------------
# read raw files and save matrix, qc
# ----------------------------------
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/proc/read10x.R')
read10x(list.files('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/raw',full.names = T),'/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc', 'HumanLymphNode',verbose = T)

rm(list=ls())

# -------
# plot qc
# -------
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/proc/qcplot.R')
qcplot(datadir = '/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/qc', plotfn = '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/01_process/plot/qcplot.pdf', width = 7, height = 2.5, num.expressed.gene.cutoff1 = 500, num.expressed.gene.cutoff2 = NULL, mito.prop.cutoff = 0.2, totalcount.cutoff =NULL)
  
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

# ----
# QC
# ----
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/proc/filtercount.R')
filtercount(savedir = '/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc', sep='_',mingn=500,maxgn=NULL,minrc=NULL,maxrc=NULL,maxmito=0.2,mingeneprop=0.01) 
rm(list=ls())

# -----
# scran
# ------
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/proc/runscran.R')
runscran(savedir = '/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc', normmat=TRUE) 

# -----
# saver
# -----
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/proc/runsaver.R')
runsaver(savedir = '/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc', count=NULL,sep='_',pid=1,ncores=20,multiple.patient=FALSE)

# ----
# PCA
# ----
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
mat <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/matrix/scran.rds')
tmp <- PCA(genebycell_mat = mat, save.pca = TRUE, plot.statistics=TRUE, plot.dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/plot/', result.dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/result/pca/', PC_for_columns = TRUE, findVariableGenes = TRUE, maxVariableGenes = 3000, numPC = NULL, smoothFittingMethod = 'loess')
rm(tmp)

# ---- 
# umap
# ----
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
mat <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/result/pca/pr.rds')
UMAP(samplebyfeature_mat = mat, save.umap = TRUE, result.dir = '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/02_analysis/result')

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

# ---------------------------
# train DNAm prediction model
# ---------------------------
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/predict.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/trainmodel.R')
meth <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/combine/wgbs/sampleme.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/combine/wgbs/filterge.rds')
genenametb <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg38_geneid_genename_genelength.rds')
genenametb[,1] <- gsub('\\..*','', genenametb[,1])
rownames(expr) <- paste0(genenametb[match(rownames(expr), genenametb[,1]) , 2], ';', rownames(expr))
m <- trainmodel(expr, meth)
saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs/model.rds')

# ------------
# predict DNAm
# ------------
m <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs/model.rds')
mat <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/saver/saver.rds')
genenametb <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg38_geneid_genename_genelength.rds')
rownames(mat) <- paste0(rownames(mat), ';', genenametb[match(rownames(mat), genenametb[,2]),1])
pred <- predict(mat, m)
saveRDS(pred, '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/03_pred_DNAm/result/predicted_DNAm.rds')

# -----------------------
# evaluate predicted DNAm: to do
# -----------------------

scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- sqrt((rowMeans(data*data) - cm^2) / (ncol(data) - 1) * ncol(data))
  (data - cm) / csd
}

corfunc <- function(m1,m2,type='concordant') {
  if (type=='concordant') {
    rowSums(scalematrix(m1) * scalematrix(m2))/(ncol(m1)-1)
  } else {
    scalematrix(t(m1)) %*% t(scalematrix(t(m2)))/(nrow(m1)-1)            
  }
}
cv <- corfunc(pred, meth[,ds %in% testds])
print(summary(cv))
cv <- corfunc(t(pred),t(meth[,ds %in% testds]))
print(summary(cv))
}

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
# ---------------------------------------------------------------------------------------
# visium RNAseq: identify cell types (use reference data: bulk RNA-seq of different cell types from ENCODE. software: singleR)
# ---------------------------------------------------------------------------------------
genenametb <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg38_geneid_genename_genelength.rds')
ref.expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/combine/wgbs/filterge.rds')
genenametb[,1] <- gsub('\\..*','', genenametb[,1])
rownames(ref.expr) <- paste0(genenametb[match(rownames(ref.expr), genenametb[,1]) , 2], ';', rownames(ref.expr))

que.expr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/saver/saver.rds')
rownames(que.expr) <- paste0(rownames(que.expr), ';', genenametb[match(rownames(que.expr), genenametb[,2]), 1])

int <- intersect(rownames(que.expr), rownames(ref.expr))
ref.expr <- ref.expr[int, ]
que.expr <- que.expr[int,]


source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
vargene.ref <- findVariableGene(genebycell_mat = ref.expr, num.gene = NULL ,plot.statistics=TRUE, plot.fn = '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/04_identifycelltype/plot/vargene_ref.pdf')
vargene.que <- findVariableGene(genebycell_mat = que.expr, num.gene = NULL ,plot.statistics=TRUE, plot.fn = '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/04_identifycelltyp/plot/vargene_que.pdf')
int <- intersect(vargene.ref, vargene.que)
ref.expr <- ref.expr[int, ]
que.expr <- que.expr[int,]
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
cormat <- corfunc(ref.expr, que.expr)
saveRDS(cormat, '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/04_identifycelltype/result/bulk_sc_expr_pearson_correlation_matrix.rds')

library(data.table)
ref.expr.rank <- apply(ref.expr,2,frank)
que.expr.rank <- apply(que.expr,2,frank)
cormat <- corfunc(ref.expr.rank, que.expr.rank)
saveRDS(cormat, '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/04_identifycelltype/result/bulk_sc_expr_spearman_correlation_matrix.rds')

ct=apply(cormat,2,which.max)
table(rownames(cormat)[ct])
# encode:spleen_male_child_(3_years) 
#                               4021 


###
library(Seurat)
library(ggplot2)
setwd('/home-4/whou10@jhu.edu/scratch/Wenpin/alsf_filbin/identifyCellType/')
# setwd('/users/whou/alsf_filbin/identifyCellType/')
d = readRDS(file='./data/filter.rds')
expr = d@assays$RNA@counts
df = readRDS('./data/allen_tumor_meta.rds')

ref.data = as.matrix(expr[, df$sample == 'allen'])
label = df[df$sample=='allen', 'type']
test = as.matrix(expr[, df$group == 'tumor_NonMalignant'])

saveRDS(ref.data, './data/refdata.rds')
saveRDS(label, './data/label.rds')
saveRDS(test, './data/test.rds')

# [jhpce01 /users/whou/alsf_filbin/identifyCellType/code]$ cat 02_singleR.R 
setwd('/users/whou/alsf_filbin/identifyCellType/')
ref.data = readRDS('./data/refdata.rds')
label =  readRDS('./data/label.rds')
test = readRDS('./data/test.rds')

library(SingleR)
library(ggplot2)
pred = SingleR(test = test, ref = ref.data, labels = colnames(ref.data))
pred2 = do.call(cbind,as.matrix(pred)@listData)
pred2 = pred@listData$scores 
rownames(pred2) = colnames(test)
donor = sub('-.*','',colnames(test))

labels = names(sort(table(pred$labels),decreasing = T)[1:20])
pdf('/users/whou/alsf_filbin/identifyCellType/plot/allen_hm.pdf',width=12,height=12.5)
plotScoreHeatmap(pred, clusters=donor)
dev.off()
# ---------------------------------------------------------------------------
# predicted DNAm: identify differential CpG using T cell and B cell bulk wgbs
# ---------------------------------------------------------------------------

# --------------------------------------------------------------------------------
# plot the differential CpG signals on visium slide, to show both predicted DNAm and bulk WGBS side by side
# --------------------------------------------------------------------------------
