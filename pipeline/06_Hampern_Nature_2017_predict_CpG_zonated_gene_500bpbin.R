# -----------
# train model
# -----------
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/predict.R')
source('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/software/trainmodel.R')
# meth <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/sampleme.rds')
expr <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/expr.rds')
genenametb <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/grcm38_geneid_genename_genelength.rds')
rownames(expr) <- genenametb[match(rownames(expr), genenametb[,1]), 2]
expr <- expr[!is.na(rownames(expr)),]
expr <- expr[rowSums(expr >= 1) >= 1,  ]
meth <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/meth.rds')

predexpr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/Hampem_Nature_2017/matrix/average_expr.rds')
predexpr <- predexpr[rowSums(predexpr)  > 0,]
int <- intersect(rownames(predexpr),rownames(expr))
predexpr = predexpr[int, ]
predexpr <- predexpr[order(apply(predexpr,1,cor,seq(1,9))),]
m <- trainmodel(expr[int,], summeth)
saveRDS(m, '/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs_mm10/model_zonegeneagg_filtergene_pred_CpG.rds')
# -----------------
# predict CpG DNAm
# -----------------
pred <- predict(predexpr[int,], m)
saveRDS(pred, '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/CpG/01_pred_DNAm/result/predicted_CpG_DNAm_using_zonated_filtergene.rds')
# --------------
# filter and plot
# --------------
mysd <- function(mat){
   ## calculate the sd for mat rows.
   m <- rowMeans(mat)
   s <- sqrt((rowMeans(mat*mat) - m^2) * ncol(mat)/(ncol(mat)-1))
}
pred.sd <- mysd(pred)
mark.sd <- (pred.sd > 0.01 & complete.cases(pred)) ###
filterpred <- pred[mark.sd, ]

library(pheatmap)
library(viridis)
pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/CpG/01_pred_DNAm/plot/CpG_predDNAm_by_layer_hm_sample1e4CpG.pdf', width = 6, height = 8)
set.seed(12345)
id = sample(1:nrow(filterpred), 1e4)
pred.sp <- filterpred[id, ]
pred.sp <- pred.sp[oder(apply(pred.sp, 1, cor, seq(1,9))), ]
pheatmap(pred.sp, cluster_cols = FALSE, scale = 'row',
       show_rownames = F, main = 'predDNAm using zonated genes',
       colors = viridis(100))
dev.off()

library(GenomicRanges)
gr <- readRDS('/home-4/whou10@jhu.edu/data2/whou10/metpred/data/data/mm10/mm10filtermat/gr.rds')
gr <- gr[mark.sd]

## DNAm on gene upstream bin
df <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/grcm38_geneid_genename_genelength.rds')
gene_gr <- GRanges(seqnames = df[,'chr'],
        ranges = IRanges(start = df[,'start'], end = df[,'end']),
        strand = df[,'strand'])
names(gene_gr) <- df[,'genename']
# gene_gr <- promoters(gene_gr,upstream=200,downstream=100)
zonated.genes <- read.csv('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/raw/Hampem_Nature_2017/NIHMS70855-supplement-Supplementary_Table_3._qvalues.csv', as.is = TRUE)
zonated.genes <- zonated.genes[,1]
zonated.genes <- intersect(names(gene_gr), zonated.genes)
zonated.gene.gr <- gene_gr[zonated.genes]

ove <- as.matrix(findOverlaps(gr, zonated.gene.gr))
ove <- data.frame(ove[,1],names(zonated.gene.gr)[ove[,2]], stringsAsFactors = FALSE)
summeth <- rowsum(filterpred[ove[,1],],ove[,2])
tab <- table(ove[,2])
v <- as.vector(tab)
names(v) <- names(tab)
summeth <- summeth/v[rownames(summeth)]
summeth <- summeth[rowSums(summeth) > 0,]
int = intersect(rownames(predexpr), rownames(summeth))
print(length(int))
plotlist = list()
x <- pheatmap(predexpr[int, ],
        cluster_rows = FALSE,cluster_cols = FALSE, 
        show_rownames = F, main = paste0(length(int), ' zonated genes expression'),
        colors = viridis(100))
plotlist[[1]] <- x[[4]]
x <- pheatmap(summeth[int, ], scale = 'row',
        cluster_rows = FALSE,cluster_cols = FALSE, 
        show_rownames = F, main = paste0('average genebody predDNAm using zonated genes'),
        colors = viridis(100))
plotlist[[2]] <- x[[4]]


g<-do.call(grid.arrange,plotlist)
ggsave(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/CpG/01_pred_DNAm/plot/', 'genebody.pdf'),g, width=6, height = 10)

library(pheatmap)
library(gridExtra)
library(viridis)

for (up in c(1:20)*500){
  print(up)
  gene_gr.start = ifelse(as.vector(strand(gene_gr))=='+', start(gene_gr)-up, end(gene_gr)+up-500)
  gene_gr.end = ifelse(as.vector(strand(gene_gr))=='+', start(gene_gr)-up+500,end(gene_gr)+up)
  gene_gr.pro = GRanges(seqnames=as.character(seqnames(gene_gr)),IRanges(start=gene_gr.start,end=gene_gr.end),strand=as.vector(strand(gene_gr)))
  names(gene_gr.pro) <- names(gene_gr)
  
  zonated.gene.gr.pro = gene_gr.pro[zonated.genes]
  ove <- as.matrix(findOverlaps(gr, zonated.gene.gr.pro))
  ove <- data.frame(ove[,1],names(zonated.gene.gr)[ove[,2]], stringsAsFactors = FALSE)
  summeth <- rowsum(filterpred[ove[,1],],ove[,2])
  tab <- table(ove[,2])
  v <- as.vector(tab)
  names(v) <- names(tab)
  summeth <- summeth/v[rownames(summeth)]  ## average CpG intensity
  summeth <- summeth[rowSums(summeth) > 0,]
  int = intersect(rownames(predexpr), rownames(summeth))
  print(length(int))
  plotlist = list()
  x <- pheatmap(predexpr[int, ],
          cluster_rows = FALSE,cluster_cols = FALSE, 
          show_rownames = F, main = paste0(length(int), ' zonated genes expression'),
          colors = viridis(100))
  plotlist[[1]] <- x[[4]]
  x <- pheatmap(summeth[int, ], scale = 'row',
          cluster_rows = FALSE,cluster_cols = FALSE, 
          show_rownames = F, main = paste0('up', as.character(up+500), 'up', as.character(up),' bin of predDNAm using zonated genes'),
          colors = viridis(100))
  plotlist[[2]] <- x[[4]]
  
  
  g<-do.call(grid.arrange,plotlist)
  ggsave(paste0('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/CpG/01_pred_DNAm/plot/', 'up', as.character(up+500), 'up', as.character(up), '.pdf'),g, width=6, height = 10)
}
  
