suppressMessages(library(Hmisc))
suppressMessages(library(rhdf5))

# read input
seed<-as.numeric(commandArgs()[7])

# read options
config<-read.table("config",h=T,colClasses=c(cksum="character",mat="character"))

# load normalize matrix
mat<-readRDS("mat")
if ( ncol(config) > 10 ) {
  sec=T
  second<-readRDS("second")
  net_name<-expand.grid(colnames(mat),colnames(second))
} else {
  sec=F
}

# append environmental parameters if necessary
if ( ! is.na(config$env)) {
	env<-readRDS("env")
	mat<-cbind(mat,env)
}

# add noise to matrix
b<-1e-4
set.seed(seed+12345)
mat_noise<-mat+(-b+(2*b)*matrix(runif(length(c(mat))),nrow=nrow(mat)))

if ( sec ) {
  set.seed(seed+12345)
  second_noise<-second+(-b+(2*b)*matrix(runif(length(c(second))),nrow=nrow(second)))
  size<-ncol(mat)*ncol(second)
  blocks<-ceiling(size/1e5)
} else {
  size<-ncol(mat)*(ncol(mat)-1)/2
  blocks<-ceiling(size/1e5)
}

# Spearman's rho
if ( sec ) {
  cor_noise<-rcorr(cbind(mat_noise,second_noise),type="spearman")
  cor_noise_r<-matrix(c(round(c(cor_noise$r[1:ncol(mat),(ncol(mat)+1):(ncol(mat)+ncol(second))]),digits=4),rep(NA,blocks*1e5-size)),nrow=1e5)
  cor_noise_p<-matrix(c(round(c(cor_noise$P[1:ncol(mat),(ncol(mat)+1):(ncol(mat)+ncol(second))]),digits=4),rep(NA,blocks*1e5-size)),nrow=1e5)
} else {
  cor_noise<-rcorr(mat_noise,type="spearman")
  cor_noise_r<-matrix(c(round(c(as.dist(cor_noise$r)),digits=4),rep(NA,blocks*1e5-size)),nrow=1e5)
  cor_noise_p<-matrix(c(round(c(as.dist(cor_noise$P)),digits=4),rep(NA,blocks*1e5-size)),nrow=1e5)
}

# save noise matrix results into blocks with the hdf5 format (1e5 values per block)
for(i in c("r","p")) {
  tmp_path<-file.path(paste0("spearman_noise_",i),paste0(seed,".h5"))
  h5createFile(tmp_path)
  h5createDataset(tmp_path,as.character(seed),c(1e5,blocks),storage.mode="double",chunk=c(1e5,1),level=7)
  h5write(get(paste0("cor_noise_",i)),file=tmp_path,name=as.character(seed))
}
