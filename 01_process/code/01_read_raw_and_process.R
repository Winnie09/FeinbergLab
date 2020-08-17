# ----------------------------------
# read raw files and save matrix, qc
# ----------------------------------
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/proc/read10x.R')
read10x(list.files('/home-4/whou10@jhu.edu/data2/whou10/spatial/data/raw',full.names = T),'/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc', 'HumanLymphNode',verbose = T)

rm(list=ls())

