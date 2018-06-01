suppressMessages(library(plyr))
suppressMessages(library(doParallel))
suppressMessages(library(rhdf5))

# read input
ncores<-as.numeric(commandArgs()[7])
blocks<-as.numeric(commandArgs()[8])

# read options
config<-read.table("config",h=T,as.is=2)
for (i in 3:5){assign(names(config)[i],config[1,i])}

# load edges
files<-paste("edges",1:blocks,sep="_")
files_ex<-paste("edges_ex",1:blocks,sep="_")

cl<-makeCluster(ncores)
edges_ok_signif<-parLapply(cl,files,function(i) {readRDS(file.path("spearman_noise_p",i))})
edges_ex_ok_signif<-parLapply(cl,files_ex,function(i) {readRDS(file.path("spearman_noise_p",i))})
stopCluster(cl)

# network
otus<-scan("otus",quiet=T)
net_name<-foreach(i=1:blocks,.combine=c) %do% {
  h5read("net_name.h5","name",index=list(edges_ok_signif[[i]],i))
}
net_edges<-ldply(net_name,function(x) otus[as.numeric(unlist(sub("^0*","",unlist(strsplit(x,"-")))))])

net_ex_name<-foreach(i=1:blocks,.combine=c) %do% {
  h5read(paste0(project,"net_name.h5"),"name",index=list(edges_ex_ok_signif[[i]],i))
}
net_ex_edges<-ldply(net_ex_name,function(x) otus[as.numeric(unlist(sub("^0*","",unlist(strsplit(x,"-")))))])

# save
write.table(net_edges,file.path("..",paste("cooccurrence",cksum,sub("\\.[^\\.]*$","",basename(config$mat)),"txt",sep=".")),col.names=F,row.names=F,sep=" ",quote=F)
write.table(net_ex_edges,file.path("..",paste("coexclusion",cksum,sub("\\.[^\\.]*$","",basename(config$mat)),"txt",sep=".")),col.names=F,row.names=F,sep=" ",quote=F)
