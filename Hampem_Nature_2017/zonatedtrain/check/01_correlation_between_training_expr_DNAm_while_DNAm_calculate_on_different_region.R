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
expr <- expr[rowSums(expr >= 1) >= 1,]

meth <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/data/mm10/meth.rds')
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

m = readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/metpred/model/wgbs_mm10/model_zonegeneagg_filtergene_genebody.rds')

### in training data, expr and meth cor
corlist = list()
int2 <- intersect(rownames(summeth), rownames(expr))
corvec2 <- sapply(int2, function(i){
  cor(summeth[i, ], expr[i, ], method = 'spearman')
})
corvec.genebody = corvec2
corlist[['genebody']] <- corvec.genebody
# ----------------------
### use promoter regions
# ----------------------

for (up in c(200, 300, 400, 500, 1000, 2000)){
  print(up)
  for (down in c(0, 100, 200, 500, 1000)){
    print(down)
    if (up > down){
      gene_gr.pro = promoters(gene_gr, upstream = up, downstream = down)
        zonated.gene.gr.pro = gene_gr.pro[zonated.genes]
        ove <- as.matrix(findOverlaps(gr, zonated.gene.gr.pro))
        ove <- data.frame(ove[,1],names(zonated.gene.gr)[ove[,2]], stringsAsFactors = FALSE)
        summeth <- rowsum(meth[ove[,1],],ove[,2])
        
        tab <- table(ove[,2])
        v <- as.vector(tab)
        names(v) <- names(tab)
        summeth <- summeth/v[rownames(summeth)]
        summeth <- summeth[rowSums(summeth) > 0,]
        
        int2 <- intersect(rownames(summeth), rownames(expr))
        corvec2 <- sapply(int2, function(i){
          cor(summeth[i, ], expr[i, ], method = 'spearman')
        })
        corlist[[paste0('up', up, 'down', down)]] <- corvec2
    }
  }
}

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/spatial/Hampem_Nature_2017/zonatedtrain/check/correlation_between_training_expr_DNAm_while_DNAm_calculate_on_different_region.pdf', width = 12, height = 7)
par(mfrow=c(4,6))
for (i in names(corlist)){
  hist(corlist[[i]], col = 'grey', breaks = 100, 
       xlab = 'SCC(training expr and DNAm)', ylab = 'num.genes',
       main = i)
}
dev.off()
