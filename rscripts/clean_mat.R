library(rhdf5)
library(vegan)
library(Hmisc)
library(igraph)
library(plyr)
library(tidyr)
library(dplyr)
library(foreach)
library(doParallel)

# read options
config<-read.table("config",h=T,as.is=2)
mat<-read.table(config$mat,h=T)
# transpose matrix if necessary
if(nrow(mat)>ncol(mat)) {
	mat<-t(mat)
}
write(nrow(mat),"nbsamp_ori")
write(ncol(mat),"nbotu_ori")
for (i in 4:ncol(config)){assign(names(config)[i],config[1,i])}

# Drop samples with too few reads
if(mincount>1) {
  minc<-mincount
} else {
  minc<-ceiling(median(rowSums(mat))*mincount)
}
mat_tmp<-mat[which(rowSums(mat)>=minc),]
mat_tmp2<-mat_tmp[,which(colSums(mat_tmp)>0)]

# Drop OTUs with low occurrence
if(minocc>1) {
  mino<-minocc
} else {
  mino<-ceiling(nrow(mat_tmp2)*minocc)
}
mat_ab<-mat_tmp2[,which(colSums((mat_tmp2>0)*1)>=mino)]

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

# Environmental table check and export
if ( ! is.na(config$env)) {
	env<-read.table(file.path(config$env),h=T)
	if ( nrow(env) != nrow(mat) ) {
		quit(save="no",status=2,runLast=F)
	}
	if ( ! all(apply(env,2,is.numeric)) ) {
		quit(save="no",status=3,runLast=F)
	}
	env_ab<-env[rownames(mat_ab),]
	if( nrow(env_ab) != nrow(mat_ab) ) {
		quit(save="no",status=4,runLast=F)
	}
	saveRDS(as.matrix(env_ab),"env")
	write(ncol(env_ab),"nbenv")
}

