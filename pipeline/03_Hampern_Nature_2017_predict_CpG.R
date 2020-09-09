# -----------------------
# transfer to data matrix
# -----------------------
tb <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/raw/Hampem_Nature_2017/NIHMS70855-supplement-Supplementary_Table_3.csv', as.is = TRUE)
rownames(tb) <- tb[,1]
tb <- tb[,-1]
dn <- list(rownames = rownames(tb), colnames = colnames(tb))
tb <- as.matrix(tb)
dimnames(tb) <- dn
saveRDS(tb, '/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/Hampem_Nature_2017/matrix/average_expr.rds')

# check if values are same as in paper Fig.3
g <- c('Axin2','Hamp', 'Arg1')
par(mfrow=c(1,3))
for (i in g){
  print(plot(tb[i,] ~ seq(1,9), pch = 20))
}

# --------------------------------  
# prepare DNAm training model data
# --------------------------------
meth <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/data/data/mm10/mm10filtermat/mat.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/data/data/mm10/TPM.rds')
colnames(expr) <- sub('.tsv', '', colnames(expr))
mm10experiment <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/data/data/mm10/meta/mm10_experiment.rds')
mm10file <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/data/data/mm10/meta/mm10_file.rds')
mm10experiment[, 'Biosample summary'] <- gsub(' ', '_', mm10experiment[, 'Biosample summary'])
acc <- mm10experiment[match(colnames(meth), mm10experiment[, 'Biosample summary']), 'Accession']
facc <- mm10file[match(acc, mm10file[,'Experiment accession']), 'File accession']
expr <- expr[, facc]
colnames(expr) <- colnames(meth)
saveRDS(meth, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/meth.rds')
saveRDS(expr, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/expr.rds')

 
# ---------------------------
# train DNAm prediction model
# ---------------------------
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/predict.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/trainmodel.R')
meth <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/sampleme.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/expr.rds')

genenametb <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/grcm38_geneid_genename_genelength.rds')
rownames(expr) <- genenametb[match(rownames(expr), genenametb[,1]), 2]
expr <- expr[!is.na(rownames(expr)),]
m <- trainmodel(expr, meth)
saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs_mm10/model_sub.rds')

meth <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/meth.rds')
m <- trainmodel(expr, meth)
saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs_mm10/model.rds')

# ------------
# predict DNAm
# ------------
# use sub model
m <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs_mm10/model_sub.rds')
mat <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/Hampem_Nature_2017/matrix/average_expr.rds')
pred <- predict(mat, m)
saveRDS(pred, '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/Hampem_Nature_2017/01_pred_DNAm/result/predicted_DNAm_using_submodel.rds')
pred <- pred[order(apply(pred,1,cor,seq(1,9))),]
library(pheatmap)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/Hampem_Nature_2017/01_pred_DNAm/plot/CpG_by_layer_hm.pdf', width = 3.5, height = 7)
pheatmap(pred, cluster_cols = FALSE, scale = 'row', cluster_rows = FALSE)
dev.off()

# use foll model
m <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/model/wgbs_mm10_model.rds')
mat <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/Hampem_Nature_2017/matrix/average_expr.rds')
pred <- predict(mat, m)
saveRDS(pred, '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/01_pred_DNAm/result/predicted_DNAm_using_fullmodel.rds')

pred <- pred[complete.cases(pred), ]
corvec <- apply(pred,1,cor,seq(1,9))
set.seed(12345)
pred.sub <- pred[sample(1:nrow(pred), 1e4), ]
pred.sub <- pred.sub[order(apply(pred.sub,1,cor,seq(1,9))),]
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/01_pred_DNAm/plot/CpG_by_layer_hm_using_fullmodel_sub1e4CpG.pdf', width = 4, height = 7)
pheatmap(pred.sub, cluster_cols = FALSE, scale = 'row', cluster_rows = FALSE,
         main = 'Predicted CpG DNAm')
dev.off()



