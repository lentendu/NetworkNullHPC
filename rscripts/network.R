suppressMessages(library(plyr))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(rhdf5))

# read input
ncores<-as.numeric(commandArgs()[7])
blocks<-as.numeric(commandArgs()[8])

# read options
config<-read.table("config",h=T,colClasses=c(cksum="character",mat="character"))
for (i in 4:ncol(config)){assign(names(config)[i],config[1,i])}

# load edges
files<-paste("edges",1:blocks,sep="_")
files_ex<-paste("edges_ex",1:blocks,sep="_")

cl<-makeCluster(ncores)
registerDoParallel(cl)
edges_ok_signif<-foreach(i=files) %dopar% {readRDS(file.path("spearman_noise_p",i))}
edges_ex_ok_signif<-foreach(i=files_ex) %dopar% {readRDS(file.path("spearman_noise_p",i))}
stopCluster(cl)

# network
otus<-scan("otus",what="character",quiet=T)
if(length(unlist(edges_ok_signif))>0) {
  net_name<-foreach(i=1:blocks,.combine=c) %do% {
	if (length(edges_ok_signif[[i]])>0) {
	  h5read("net_name.h5","name",index=list(edges_ok_signif[[i]]$edges,i))
	}
  }
  net_edges<-cbind(ldply(net_name,function(x) otus[as.numeric(unlist(sub("^0*","",unlist(strsplit(x,"-")))))]),ldply(edges_ok_signif)$cor)
  write.table(net_edges,file.path("..",paste("cooccurrence",config$cksum,sub("\\.[^\\.]*$","",basename(config$mat)),"txt",sep=".")),col.names=F,row.names=F,sep=" ",quote=F)
  outcooc<-0
} else {
  outcooc<-2
}

if(length(unlist(edges_ex_ok_signif))>0) {
  net_ex_name<-foreach(i=1:blocks,.combine=c) %do% {
    if (length(edges_ex_ok_signif[[i]])>0) {
      h5read("net_name.h5","name",index=list(edges_ex_ok_signif[[i]]$edges,i))
    }
  }
  net_ex_edges<-cbind(ldply(net_ex_name,function(x) otus[as.numeric(unlist(sub("^0*","",unlist(strsplit(x,"-")))))]),ldply(edges_ex_ok_signif)$cor)
  write.table(net_ex_edges,file.path("..",paste("coexclusion",config$cksum,sub("\\.[^\\.]*$","",basename(config$mat)),"txt",sep=".")),col.names=F,row.names=F,sep=" ",quote=F)
  outcoex<-0
} else {
  outcoex<-3
}

quit(save="no",status=(outcooc+outcoex))
