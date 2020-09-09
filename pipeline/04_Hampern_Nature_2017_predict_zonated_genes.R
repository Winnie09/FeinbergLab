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
# meth <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/sampleme.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/expr.rds')

genenametb <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/grcm38_geneid_genename_genelength.rds')
rownames(expr) <- genenametb[match(rownames(expr), genenametb[,1]), 2]
expr <- expr[!is.na(rownames(expr)),]
# m <- trainmodel(expr, meth)
# saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs_mm10/model_sub.rds')
meth <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/meth.rds')

expr <- expr[rowSums(expr >= 1) >= 1,]

predexpr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/Hampem_Nature_2017/matrix/average_expr.rds')
predexpr <- predexpr[rowSums(predexpr) > 0,]
int <- intersect(rownames(predexpr),rownames(expr))

library(GenomicRanges)
gr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/data/data/mm10/mm10filtermat/gr.rds')
df <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/grcm38_geneid_genename_genelength.rds')
gene_gr <- GRanges(seqnames = df[,'chr'],
        ranges = IRanges(start = df[,'start'], end = df[,'end']),
        strand = df[,'strand'])
names(gene_gr) <- df[,'genename']
gene_gr <- promoters(gene_gr,upstream=200,downstream=100)

zonated.genes <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/raw/Hampem_Nature_2017/NIHMS70855-supplement-Supplementary_Table_3._qvalues.csv', as.is = TRUE)
zonated.genes <- zonated.genes[,1]
zonated.genes <- intersect(names(gene_gr), zonated.genes)
zonated.gene.gr <- gene_gr[zonated.genes]
ove <- as.matrix(findOverlaps(gr, zonated.gene.gr))

ove <- data.frame(ove[,1],names(zonated.gene.gr)[ove[,2]], stringsAsFactors = FALSE)

summeth <- rowsum(meth[ove[,1],],ove[,2])
tab <- table(ove[,2])
v <- as.vector(tab)
names(v) <- names(tab)
summeth <- summeth/v[rownames(summeth)]
summeth <- summeth[rowSums(summeth) > 0,]

m <- trainmodel(expr[int,], summeth)
saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs_mm10/model_zonegeneagg_filtergene_up200down100.rds')

pred <- predict(predexpr[int,], m)
sum(apply(pred,1,sd) > 0.01)
filterpred <- pred[apply(pred,1,sd) > 0.01,]
saveRDS(pred, '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/promoter/01_pred_DNAm/result/predicted_DNAm_using_zonageneagg_filtergene_up200down100.rds')

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/promoter/01_pred_DNAm/plot/zonageneagg_by_layer_hm_up200down100.pdf', width = 6, height = 8)
pheatmap(pred, cluster_cols = FALSE, scale = 'row',
         show_rownames = F, main = 'zonated  genes aggregated predDNAm (700 sd > 0.01)')
dev.off()

corvec <- sapply(1:nrow(filterpred),function(i) {
  cor(predexpr[rownames(filterpred),][i,],filterpred[i,],method='spearman')
})

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/promoter/01_pred_DNAm/plot/zonageneagg_by_layer_cv_up200down100.pdf', width = 6, height = 4)
hist(corvec, col='grey', breaks = 100, xlab = 'SCC(filterpred.DNAm, filterpred.expr)', ylab = 'Number of zonated genes', main = 'train&pred aggreated DNAm of zonated genes promoters')
dev.off()


# ------------
# predict DNAm
# ------------
# use sub model
m <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs_mm10/model_sub.rds')
mat <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/Hampem_Nature_2017/matrix/average_expr.rds')
pred <- predict(mat, m)
saveRDS(pred, '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/promoter/01_pred_DNAm/result/predicted_DNAm_using_submodel.rds')
pred <- pred[order(apply(pred,1,cor,seq(1,9))),]
library(pheatmap)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/promoter/01_pred_DNAm/plot/CpG_by_layer_hm.pdf', width = 3.5, height = 7)
pheatmap(pred, cluster_cols = FALSE, scale = 'row', cluster_rows = FALSE)
dev.off()

# use foll model
m <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/model/wgbs_mm10_model.rds')
mat <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/Hampem_Nature_2017/matrix/average_expr.rds')
pred <- predict(mat, m)
saveRDS(pred, '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/promoter/01_pred_DNAm/result/predicted_DNAm_using_fullmodel.rds')


## findoverlap of predicted CpG sites and zonated genes promoters 
cm <- rowMeans(pred)
csv <- sqrt((rowMeans(pred*pred) - cm^2) / (ncol(pred) - 1) * ncol(pred))



