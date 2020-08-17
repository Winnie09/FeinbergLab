# --
# QC
# --
source('/home-4/whou10@jhu.edu/scratch/Wenpin/raisin/proc/filtercount.R')
filtercount(savedir = '/home-4/whou10@jhu.edu/data2/whou10/spatial/data/proc', sep='_',mingn=500,maxgn=NULL,minrc=NULL,maxrc=NULL,maxmito=0.2,mingeneprop=0.01) 
rm(list=ls())

