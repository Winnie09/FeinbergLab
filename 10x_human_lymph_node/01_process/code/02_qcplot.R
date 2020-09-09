# -------
# plot qc
# -------
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/proc/qcplot.R')
qcplot(datadir = '/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc/qc', plotfn = '/home-4/whou10@jhu.edu/scratch/Wenpin/spatial10x/01_process/plot/qcplot.pdf', width = 7, height = 2.5, num.expressed.gene.cutoff1 = 500, num.expressed.gene.cutoff2 = NULL, mito.prop.cutoff = 0.2, totalcount.cutoff =NULL)

