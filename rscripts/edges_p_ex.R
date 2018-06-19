suppressMessages(library(foreach))
suppressMessages(library(rhdf5))

# read input
block<-as.numeric(commandArgs()[7])

# read options
config<-read.table("config",h=T,as.is=2)
for (i in 4:ncol(config)){assign(names(config)[i],config[1,i])}

# get p.values of edges above threshold
edges_ex_ok<-readRDS(file.path("spearman_noise_r",paste0("edges_ex_",block)))

if (length(edges_ex_ok)>0) {
  # read p.values
  pvals<-foreach(i=1:nboot,.combine=cbind) %do% {
    tmp_path<-file.path("spearman_noise_p",paste0(i,".h5"))
    h5read(tmp_path,as.character(i),index=list(edges_ex_ok$edges,block))
  }
  # Select edges for which 90% of corrected pvalues are significant
  pvals_ex_ok<-which(apply(pvals,1,function(x) {sum((p.adjust(x,method="BH")<=0.01)*1)})>=nboot*0.9)
  # intersection of Markov sampling and threshold validated edges with significant edges
  edges_ex_ok_signif<-edges_ex_ok[pvals_ex_ok]
} else {
  edges_ex_ok_signif<-NULL
}

# save edges
saveRDS(edges_ex_ok_signif,file.path("spearman_noise_p",paste0("edges_ex_",block)))
