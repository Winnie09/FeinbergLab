source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/predict.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/trainmodel.R')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/expr.rds')
genenametb <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/grcm38_geneid_genename_genelength.rds')
rownames(expr) <- genenametb[match(rownames(expr), genenametb[,1]), 2]
expr <- expr[!is.na(rownames(expr)),]


pred = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/promoter/01_pred_DNAm/result/predicted_DNAm_using_zonageneagg_filtergene_up250down0.rds')
int = intersect(rownames(expr), rownames(pred))

predexpr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/Hampem_Nature_2017/matrix/average_expr.rds')
predexpr <- predexpr[rowSums(predexpr) > 0,]
predexpr <- predexpr[int, ]
predexpr <- predexpr[order(apply(predexpr,1,cor,seq(1,9))),]

library(viridis)
library(pheatmap)
mycol = viridis(100)
mycol1 = c(rep(mycol[1],10),mycol[1:95], rep(mycol[91:100], each = 5))
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/promoter/01_pred_DNAm/plot/1058_zonageneagg_by_layer_expr_hm.pdf', width = 5, height = 8)
pheatmap(predexpr, cluster_cols = FALSE, scale = 'row', cluster_rows = F,
         color = mycol2,
         show_rownames = F, main = 'zonated genes expression (1058 genes)')
dev.off()

mycol2 = c(rep(mycol[1],10),mycol[1:95], rep(mycol[96:100], each = 10))
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/promoter/01_pred_DNAm/plot/1058_zonageneagg_by_layer_predDNAm_hm_up250down0.pdf', width = 5, height = 8)
pheatmap(pred[rownames(predexpr), ], cluster_cols = FALSE, scale = 'row', cluster_rows = F,
         show_rownames = F, main = '1058 zonated genes body aggregated predDNAm',
         color = mycol2)
dev.off()

