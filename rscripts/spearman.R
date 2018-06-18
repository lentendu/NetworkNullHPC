suppressMessages(library(Hmisc))
suppressMessages(library(rhdf5))

# read input
seed<-as.numeric(commandArgs()[7])

# read options
config<-read.table("config",h=T,as.is=2)

# load normalize matrix
mat<-readRDS("mat")

# append environmental parameters if necessary
if ( ! is.na(config$env)) {
	env<-readRDS("env")
	mat<-cbind(mat,env)
}

# add noise to matrix
size<-ncol(mat)*(ncol(mat)-1)/2
blocks<-ceiling(size/1e5)
b<-1e-4
set.seed(seed+12345)
mat_noise<-mat+(-b+(2*b)*matrix(runif(length(c(mat))),nrow=nrow(mat)))

# Spearman's rho
cor_noise<-rcorr(mat_noise,type="spearman")
cor_noise_r<-matrix(c(round(c(as.dist(cor_noise$r)),digits=4),rep(NA,blocks*1e5-size)),nrow=1e5)
cor_noise_p<-matrix(c(round(c(as.dist(cor_noise$P)),digits=4),rep(NA,blocks*1e5-size)),nrow=1e5)

# save noise matrix results in splitted blocks with the hdf5 format (1e5 values per block)
for(i in c("r","p")) {
  tmp_path<-file.path(paste0("spearman_noise_",i),paste0(seed,".h5"))
  h5createFile(tmp_path)
  h5createDataset(tmp_path,as.character(seed),c(1e5,blocks),storage.mode="double",chunk=c(1e5,1),level=7)
  h5write(get(paste0("cor_noise_",i)),file=tmp_path,name=as.character(seed))
}
