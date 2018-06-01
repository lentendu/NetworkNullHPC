library(plyr)
library(tidyr)
library(dplyr)
library(vegan)
library(foreach)
library(doParallel)
library(rhdf5)
library(Hmisc)
library(igraph)

# read options
config<-read.table("config",h=T,as.is=2)
mat<-read.table(config$mat,h=T)
write(nrow(mat),"nbsamp_ori")
write(ncol(mat),"nbotu_ori")
for (i in 3:5){assign(names(config)[i],config[1,i])}

# Drop samples with too few reads
if(mincount>1) {
  minc<-mincount
} else {
  minc<-median(rowSums(mat))*mincount
}
mat_tmp<-mat[which(rowSums(mat)>=minc),]

# Drop OTUs with low occurrence
if(minocc>1) {
  mino<-minocc
} else {
  mino<-ceiling(nrow(mat_tmp)*minocc)
}
mat_ab<-mat_tmp[,which(colSums(mat_tmp>0*1)>=mino)]

# Normalize read counts by multiplying each sample read counts by the ratio of half the read count median of all samples divide by total sample read counts
# better than rarefaction (but see McMurdie & Holmes, 2014)
if(depth>1) {
  mat_norm<-round(mat_ab*(depth/rowSums(mat_ab)))
} else {
  mat_norm<-round(mat_ab*(median(rowSums(mat_ab))*depth/rowSums(mat_ab)))
}

# save normalized OTU table and its size
saveRDS(as.matrix(mat_norm),"mat")
write(nrow(mat_norm),"nbsamp")
write(ncol(mat_norm),"nbotu")
otus<-colnames(mat_norm)
write(otus,"otus",ncolumns=1)
write(mino,"minocc")
write(minc,"mincount")
