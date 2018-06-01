suppressMessages(library(foreach))
suppressMessages(library(rhdf5))

# read input
block<-as.numeric(commandArgs()[7])
threshold<-scan("threshold",quiet=T)
threshold_ex<-scan("ex_threshold",quiet=T)

# Read Spearman's rho from noise added matrix (one block at a time)
edges<-foreach(i=1:1000,.combine=cbind) %do% {
  tmp_path<-file.path("spearman_noise_r",paste0(i,".h5"))
  h5read(tmp_path,as.character(i),index=list(NULL,block))
}
# Select edges above threshold and present in at least 90 % of the Markov samplings
edges_ok<-which(rowSums((edges>threshold)*1)>=1000*0.9)
saveRDS(edges_ok,file.path("spearman_noise_r",paste0("edges_",block)))

# Exclusion edges below the exclusion threshold
edges_ex_ok<-which(rowSums((edges<threshold_ex)*1)>=1000*0.9)
saveRDS(edges_ex_ok,file.path("spearman_noise_r",paste0("edges_ex_",block)))