suppressMessages(library(plyr))
suppressMessages(library(Hmisc))
suppressMessages(library(igraph,quietly=T))

# read input
seed<-as.numeric(commandArgs()[7])

# Read options
config<-read.table("config",h=T,as.is=1)
for (i in 3:ncol(config)){assign(names(config)[i],config[1,i])}

# load normalize matrix 
mat<-readRDS("mat")
size<-ncol(mat)*(ncol(mat)-1)/2
net_name<-as.data.frame(t(combn(colnames(mat),2)))
pos_tresh<-scan("pos_tresh_range",quiet=T)
neg_tresh<-scan("neg_tresh_range",quiet=T)

# add noise to matrix
b<-1e-4
set.seed(seed+12345)
mat_noise<-mat+(-b+(2*b)*matrix(runif(length(c(mat))),nrow=nrow(mat)))

# randomize the normalize matrix with re-introduced compositionnality, following Faust et al. (2012):
# first randomize the matrix for each OTU separetedly, re-normalize per sample and add noise
set.seed(seed+12345)
mat_rand<-apply(mat,2,sample)
if(depth>1) {
  mat_rand_norm<-round(mat_rand*(depth/rowSums(mat_rand)))
} else {
  mat_rand_norm<-round(mat_rand*(median(rowSums(mat_rand))*depth/rowSums(mat_rand)))
}
set.seed(seed+12345)
mat_rand_norm_noise<-mat_rand_norm+(-b+(2*b)*matrix(runif(length(c(mat))),nrow=nrow(mat)))

# Spearman's rho
cor_rand<-rcorr(mat_rand_norm_noise,type="spearman")
cor_rand_r<-c(as.dist(cor_rand$r))

# construct network and get size of the largest connected component at different thresholds for positive correlations and negative correlations
rand_cc<-laply(seq(pos_tresh-0.05,pos_tresh+0.15,0.01),function(x){
  tmp_edges<-which(cor_rand_r>=x)
  if(length(tmp_edges)>0) {
    tmp_c<-max(components(graph.data.frame(net_name[tmp_edges,],directed=F))$csize)
  } else {tmp_c<-0}
})
rand_cc_ex<-laply(seq(neg_tresh+0.05,neg_tresh-0.15,-0.01),function(x){
  tmp_edges<-which(cor_rand_r<=x)
  if(length(tmp_edges)>0) {
    tmp_c<-max(components(graph.data.frame(net_name[tmp_edges,],directed=F))$csize)
  } else {tmp_c<-0}
})

# save
write(rand_cc,file.path("spearman_rand_r",paste0("cc_",seed)),ncolumns=1)
write(rand_cc_ex,file.path("spearman_rand_r",paste0("cc_ex_",seed)),ncolumns=1)
