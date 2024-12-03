suppressMessages(library(vegan))
suppressMessages(library(plyr))
suppressMessages(library(Hmisc))
suppressMessages(library(igraph,quietly=T))

# read input
seed<-as.numeric(commandArgs()[7])

# Read options
config<-read.table("config",h=T,colClasses=c(cksum="character",mat="character"))
for (i in 4:ncol(config)){assign(names(config)[i],config[1,i])}

# load normalize matrix 
mat<-readRDS("mat")
if ( ncol(config)>10 ) {
  sec=T
  second<-readRDS("second")
  net_name<-expand.grid(colnames(mat),colnames(second))
} else {
  sec=F
  net_name<-as.data.frame(t(combn(colnames(mat),2)))
}
pos_tresh<-scan("pos_tresh_range",quiet=T)
neg_tresh<-scan("neg_tresh_range",quiet=T)

# randomize the normalize matrix
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

if ( sec ) {
  if(is.numeric(nullm)) {
    if(nullm==0) {
      set.seed(seed+12345); second_rand<-matrix(sample(c(second)),nrow=nrow(second))
    } else if(nullm==1) {
      set.seed(seed+12345); second_rand<-t(apply(second,nullm,sample))
    } else {
      set.seed(seed+12345); second_rand<-apply(second,nullm,sample)
    }
  } else{
    set.seed(seed+12345); second_rand<-permatfull(second,fixedmar=nullm,times=1)$perm[[1]]
  }
}

# re-normalize per sample if necessary
if( sd(rowSums(mat_rand))/ncol(mat_rand)>1 & norm!="no" ) {
  if(depth>1) {
    mat_rand_norm<-round(mat_rand*(depth/rowSums(mat_rand)))
  } else {
    mat_rand_norm<-round(mat_rand*(median(rowSums(mat_rand))*depth/rowSums(mat_rand)))
  }
} else {
  mat_rand_norm<-mat_rand
}

if ( sec ) {
  if ( same && sd(rowSums(second))/ncol(second_rand)>1 ) {
    if(depth>1) {
      second_rand_norm<-round(second_rand*(depth/rowSums(second_rand)))
    } else {
      second_rand_norm<-round(second_rand*(median(rowSums(second_rand))*depth/rowSums(second_rand)))
    }
  } else {
    second_rand_norm<-second_rand
  }
}

# add noise
b<-1e-4
set.seed(seed+12345)
mat_rand_norm_noise<-mat_rand_norm+(-b+(2*b)*matrix(runif(length(c(mat))),nrow=nrow(mat)))

if ( sec ) {
  set.seed(seed+12345)
  second_rand_norm_noise<-second_rand_norm+(-b+(2*b)*matrix(runif(length(c(second))),nrow=nrow(second)))
}

# Spearman's rho
if ( sec ) {
  cor_rand<-rcorr(cbind(mat_rand_norm_noise,second_rand_norm_noise),type="spearman")
  cor_rand_r<-c(cor_rand$r[1:ncol(mat),(ncol(mat)+1):(ncol(mat)+ncol(second))])
} else {
  cor_rand<-rcorr(mat_rand_norm_noise,type="spearman")
  cor_rand_r<-c(as.dist(cor_rand$r))
}

# construct network and get size of the largest connected component at different thresholds for positive correlations and negative correlations
rand_cc<-laply(seq(pos_tresh-0.05,pos_tresh+0.15,0.01),function(x){
  tmp_edges<-which(cor_rand_r>=x)
  if(length(tmp_edges)>0) {
    tmp_c<-max(components(graph_from_data_frame(net_name[tmp_edges,],directed=F))$csize)
  } else {tmp_c<-0}
})
rand_cc_ex<-laply(seq(neg_tresh+0.05,neg_tresh-0.15,-0.01),function(x){
  tmp_edges<-which(cor_rand_r<=x)
  if(length(tmp_edges)>0) {
    tmp_c<-max(components(graph_from_data_frame(net_name[tmp_edges,],directed=F))$csize)
  } else {tmp_c<-0}
})

# save
write(rand_cc,file.path("spearman_rand_r",paste0("cc_",seed)),ncolumns=1)
write(rand_cc_ex,file.path("spearman_rand_r",paste0("cc_ex_",seed)),ncolumns=1)
