suppressMessages(library(vegan))
suppressMessages(library(foreach))
suppressMessages(library(Hmisc))
suppressMessages(library(igraph,quietly=T))
suppressMessages(library(rhdf5))

# read options
config<-read.table("config",h=T,colClasses=c(cksum="character",mat="character"))
for (i in 4:ncol(config)){assign(names(config)[i],config[1,i])}

# load normalize matrix 
mat<-readRDS("mat")
net_name<-t(combn(colnames(mat),2))

# randomize the normalize matrix
seed<-1
if(is.numeric(nullm)) {
  if(nullm==0) {
    set.seed(seed+12345); mat_rand<-matrix(sample(c(mat)),nrow=nrow(mat))
  } else if(nullm==1) {
    set.seed(seed+12345); mat_rand<-t(apply(mat,nullm,sample))
  } else {
    set.seed(seed+12345); mat_rand<-apply(mat,nullm,sample)
  }
} else{
  set.seed(seed+12345); mat_rand<-permatfull(mat,fixedmar=nullm,times=1)$perm[[1]]
}

# re-normalize per sample if necessary
if(sd(rowSums(mat))/ncol(mat_rand)>1) {
  if(depth>1) {
    mat_rand_norm<-round(mat_rand*(depth/rowSums(mat_rand)))
  } else {
    mat_rand_norm<-round(mat_rand*(median(rowSums(mat_rand))*depth/rowSums(mat_rand)))
  }
} else {
  mat_rand_norm<-mat_rand
}

# add noise
b<-1e-4
set.seed(seed+12345)
mat_rand_norm_noise<-mat_rand_norm+(-b+(2*b)*matrix(runif(length(c(mat))),nrow=nrow(mat)))

# Spearman's rho
cor_rand<-rcorr(mat_rand_norm_noise,type="spearman")
cor_rand_r<-c(as.dist(cor_rand$r))

# positive correlation
pos_tresh<-NULL
for(i in seq(0.9,0.01,-0.01)) {
  tmp_edges<-which(cor_rand_r>=i)
  if(length(tmp_edges)>1) {
    tmp_c<-max(components(graph.data.frame(net_name[tmp_edges,],directed=F))$csize)
    if(tmp_c/ncol(mat)>largecp/100) {
      pos_tresh<-i
      break
    }
  }
}

neg_tresh<-NULL
for(i in seq(-0.9,-0.01,0.01)) {
  tmp_edges<-which(cor_rand_r<=i)
  if(length(tmp_edges)>1) {
    tmp_c<-max(components(graph.data.frame(net_name[tmp_edges,],directed=F))$csize)
    if(tmp_c/ncol(mat)>largecp/100) {
      neg_tresh<-i
      break
    }
  }
}

write(pos_tresh,"pos_tresh_range")
write(neg_tresh,"neg_tresh_range")

# append environmental parameters to normalize matrix if necessary
if ( ! is.na(config$env)) {
	env<-readRDS("env")
	mat<-cbind(mat,env)
	net_name<-t(combn(colnames(mat),2))
}

# prepare spearman's rho storage
size<-ncol(mat)*(ncol(mat)-1)/2
blocks<-ceiling(size/1e5)
otusnum<-seq(1,ncol(mat))
name_length<-nchar(tail(otusnum,1))
pair_names<-combn(sprintf(paste0("%0",name_length,"d"),otusnum),2)
full_net_name<-matrix(c(paste(pair_names[1,],pair_names[2,],sep="-"),rep(NA,ceiling(size/1e5)*1e5-size)),nrow=1e5)
h5createFile("net_name.h5")
h5createDataset("net_name.h5","name",c(1e5,blocks),storage.mode="character",size=name_length*2+2,chunk=c(1e5,1),level=6)
h5write(full_net_name,file="net_name.h5",name="name")
