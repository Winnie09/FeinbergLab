## set up location: local or cluster
# setwd('/Users/wenpinhou/Dropbox/FeinbergLab/mousePancreas/wgbs/')
setwd('/home/whou10/scratch4/whou10/FeinbergLab/mousePancreas/wgbs/')
# source('/Users/wenpinhou/Dropbox/resource/ggplot_theme.R')
source('/home/whou10/scratch16/whou10/resource/ggplot_theme.R')

#####################################
## calculate promoter DNA methylation
#####################################
suppressMessages(library(bsseq))
library(data.table)
load('data/adm_bsseq_smoothed_cov_filtered.rda')
m <- assays(bsseq_smoothed_sub)@listData$M
cov <- assays(bsseq_smoothed_sub)@listData$Cov

gtf <- fread('/home/whou10/scratch16/whou10/resource/gtf/grcm38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gene <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
names(gene) <- sub('".*','',sub('.*gene_name "','',gtf[,9]))
pro <- promoters(gene,upstream=1000,downstream=0)
o <- as.matrix(findOverlaps(bsseq_smoothed_sub@rowRanges,pro))
sm <- rowsum(m[o[,1],],names(gene)[o[,2]])
scov <- rowsum(cov[o[,1],],names(gene)[o[,2]])
prop <- sm/scov
saveRDS(prop, 'promoter_methylation/promoter_tssup1kb_averaged_dnam.rds')




## box plot of promoter DNA methylation

### prepare box plot data 
library(reshape2)
library(ggplot2)
pd <- melt(prop)
meta <- read.csv('doc/all_adm_sample_metadata.csv')
pd$experiment <- meta[match(sub('_unrecovered$','',sub('_[a-zA-Z]*$','',pd$Var2)),meta[,1]),2]

pd$ct <- sub('.*_','',pd[,2])
pd$ct <- pd[,2]
saveRDS(pd, 'promoter_methylation/promoter_tssup1kb_averaged_dnam_boxplot_pd.rds')

### plot box plot
pd <- readRDS('promoter_methylation/promoter_tssup1kb_averaged_dnam_boxplot_pd.rds')
str(pd)
pd[,1] <- as.character(pd[,1])
pd[,2] <- as.character(pd[,2])
pd <- rbind(pd[pd[,5]=='acinar',],
            pd[pd[,5]=='ADM',],
            pd[pd[,5]=='duct',])
str(pd)
pd[,'ct'] <- factor(pd[,'ct'], levels = c('acinar', 'ADM', 'duct'))
pd[,2] <- factor(pd[,2], levels = unique(pd[,2]))
str(pd)
pdf('promoter_methylation/promoter_dnam_vs_sample_boxplot.pdf', width = 7, height = 2.5) ## this is strange, can't repeat the boxplot
ggplot(data = pd, aes(x = Var2, y = value, fill = ct)) +
  geom_boxplot() +
  facet_wrap(~experiment, scale='free_x') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab('WGBS experiment') + 
  ylab('Promoter DNA methylation') + 
  ggtitle('promoter: TSS upstream 1kb') +
  scale_fill_brewer(palette = 'Set2')
dev.off()
  

ggplot(pd,aes(x=Var2,fill=ct,y=value)) + geom_boxplot() + theme_classic() + facet_wrap(~experiment,scale='free_x')
dev.off()

## volcano plot of differential promoter DNA methylation 
## between Caerulein_ADM and Klf4_ADM samples
## Differential is done for each ct: adinar, ADM, duct
prop <- prop[,!grepl('unrecovered_ADM',colnames(prop))] ## remove a abnormal sample
library(limma)
pdmeta <- unique(pd[,c('Var2','ct','experiment')])
for (sct in unique(pdmeta$ct)) {
  print(sct)
  tmpmeta <- pdmeta[pdmeta$ct==sct,]
  res <- topTable(eBayes(lmFit(prop[,tmpmeta[,1]],cbind(1,tmpmeta$experiment=='Caerulein_ADM'))),n=nrow(prop),coef=2)
  res <- res[,c('logFC','P.Value','adj.P.Val')]
  colnames(res)[1] <- 'logFC_Caer-Klf4'
  ## prepare vocano plot data
  pd.diff <- data.frame(logfc = res[,'logFC_Caer-Klf4'],
                        logfdr = -log10(res$adj.P.Val),
                        sig=ifelse(res$adj.P.Val < 0.05 & abs(res$`logFC_Caer-Klf4`)>0.1,'significant','not_significant'),
                        gene = rownames(res),
                        stringsAsFactors = FALSE)
  pd.diff$gene[pd.diff[,3]=='not_significant'] <- ''
  ## plot 
  pdf(paste0('promoter_methylation/volcano_',sct,'.pdf'), width = 6, height = 6)  
  print(ggplot()  +
    geom_point(data = pd.diff,
               aes(x = logfc, 
                   y = logfdr,
                   color = sig),
               size = 1, stroke = 0)+
    ggrepel::geom_text_repel(data = pd.diff[1:100,],
                              aes(x = logfc,
                                  y = logfdr,
                                  label=gene),
                             size = 2) +
    scale_color_manual(values = c('black','red')) + 
    xlab('log2 fold change (Caer-Klf4)') + 
    ylab('log10 FDR')+
    ggtitle(paste0(sum(pd.diff=='significant'), ' out of ', nrow(res), ' are significant')))
   dev.off()
}



################
## GO analysis
################
## enrichR 
## load enrichR
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr")
## dbs <- listEnrichrDbs() ## all databses

library(limma)
pdmeta <- unique(pd[,c('Var2','ct','experiment')])
allres <- NULL
res.diff <- NULL
for (sct in unique(pdmeta$ct)) {
  print(sct)
  tmpmeta <- pdmeta[pdmeta$ct==sct,]
  res <- topTable(eBayes(lmFit(prop[,tmpmeta[,1]],cbind(1,tmpmeta$experiment=='Caerulein_ADM'))),n=nrow(prop),coef=2)
  res <- res[,c('logFC','P.Value','adj.P.Val')]
  colnames(res)[1] <- 'logFC_Caer-Klf4'
  res.diff <- rbind(res.diff, data.frame(res[res$adj.P.Val < 0.05,], celltype = sct))
    
  gl <- rownames(res)[res$adj.P.Val < 0.05 & res[,1] > 0.1] ## dhceck cutoffs
  dbs <- c("GO_Biological_Process_2015") ## omit other "GO_Molecular_Function_2015", "GO_Cellular_Component_2015"
  enriched <- enrichr(gl, dbs)
  sigres <- enriched[["GO_Biological_Process_2015"]]
  if (nrow(sigres) > 0){
    allres <- rbind(allres, data.frame(sigres,type=paste0(sct,'_Caer>Klf4')))
  }
  
  gl <- rownames(res)[res$adj.P.Val < 0.05 & res[,1] < -0.1]
  enriched <- enrichr(gl, dbs)
  sigres <- enriched[["GO_Biological_Process_2015"]]
  if (nrow(sigres) > 0){
    allres <- rbind(allres, data.frame(sigres,type=paste0(sct,'_Caer<Klf4')))
  }
}
saveRDS(res.diff, file='promoter_methylation/differential_promoter_dnam_caer_minus_klf4.rds')
write.csv(res.diff, file='promoter_methylation/differential_promoter_dnam_caer_minus_klf4.csv')
allres <- allres[allres$Adjusted.P.value < 0.1,]
saveRDS(allres, file='promoter_methylation/significantGO_enrichR.rds')
write.csv(allres, file='promoter_methylation/significantGO_enrichR.csv')


## top GO
suppressMessages(library(topGO))
library(limma)
pdmeta <- unique(pd[,c('Var2','ct','experiment')])
allres <- NULL
for (sct in unique(pdmeta$ct)) {
  print(sct)
  tmpmeta <- pdmeta[pdmeta$ct==sct,]
  res <- topTable(eBayes(lmFit(prop[,tmpmeta[,1]],cbind(1,tmpmeta$experiment=='Caerulein_ADM'))),n=nrow(prop),coef=2)
  res <- res[,c('logFC','P.Value','adj.P.Val')]
  colnames(res)[1] <- 'logFC_Caer-Klf4'

    
  back <- rownames(res)
  gl <- rownames(res)[res$adj.P.Val < 0.05 & res[,1] > 0.1] ## dhceck cutoffs
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)})
  sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
  sigres$FC <- sigres$Significant/sigres$Expected
  allres <- rbind(allres, data.frame(sigres,type=paste0(sct,'_Caer>Klf4')))

  back <- rownames(res)
  gl <- rownames(res)[res$adj.P.Val < 0.05 & res[,1] < -0.1]
  geneList <- factor(as.integer(back %in% gl))
  names(geneList) <- back
  suppressMessages({GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,geneSel=function(a) {a},annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "Symbol")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  sigres <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(resultFisher@score),orderBy="classicFisher",numChar=1000)})
  sigres$classicFisher[sigres$classicFisher=="< 1e-30"] <- 0
  sigres <- sigres[sigres$Annotated >= 10,]
  sigres$FDR <- p.adjust(sigres$classicFisher,method="fdr")
  sigres$FC <- sigres$Significant/sigres$Expected
  allres <- rbind(allres, data.frame(sigres,type=paste0(sct,'_Caer<Klf4')))

}

allres <- allres[allres$FDR < 0.05,]
saveRDS(allres, file='promoter_methylation/significantGO_topGO.rds')
write.csv(allres, file='promoter_methylation/significantGO_topGO.csv')
