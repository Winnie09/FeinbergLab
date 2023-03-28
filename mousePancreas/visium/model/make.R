library(GenomicRanges)
e <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred20210520/data/mm10/expr.rds')
m <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred20210520/data/mm10/meth.rds')
gr <- readRDS('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred20210520/data/mm10/gr.rds')

library(data.table)
gtf <- fread('/home-4/zji4@jhu.edu/scratch/resource/gtf/grcm38.gtf',data.table=F)
gtf <- gtf[gtf[,3]=='gene',]
gn <- sub('\".*','',sub('.*gene_name "','',gtf[,9]))
gid <- sub('\\..*','',sub('\".*','',sub('.*gene_id "','',gtf[,9])))
did <- duplicated(gn)
gn <- gn[!did]
gtf <- gtf[!did,]
gener <- GRanges(seqnames=gtf[,1],IRanges(start=gtf[,4],end=gtf[,5]),strand=gtf[,7])
pro <- promoters(gener,upstream=200,downstream=0)
over <- as.matrix(findOverlaps(gr,pro))
gnover <- gn[over[,2]]
tarm <- m[over[,1],]

am <- rowsum(tarm,gnover)
am <- am/as.vector(table(gnover)[rownames(am)])

rownames(e) <- sub('\\..*','',rownames(e))
e <- e[rownames(e)%in%gid,]
gne <- gn[match(rownames(e),gid)]
id <- !duplicated(gne)
e <- e[id,]
rownames(e) <- gne[id]

source('/home-4/whou10@jhu.edu/work-zfs/whou10/metpred/software/trainmodel.R')
mod <- trainmodel(e,am)
saveRDS(mod,file='/home-4/whou10@jhu.edu/work-zfs/whou10/spatial/model/model.rds')

