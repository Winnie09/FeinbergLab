# ------------
# predict DNAm
# ------------
m <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs/model.rds')
mat <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/saver/saver.rds')
genenametb <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/hg38_geneid_genename_genelength.rds')
rownames(mat) <- paste0(rownames(mat), ';', genenametb[match(rownames(mat), genenametb[,2]),1])
pred <- predict(mat, m)
saveRDS(pred, '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/03_pred_DNAm/result/predicted_DNAm.rds')

