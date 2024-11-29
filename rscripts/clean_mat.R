suppressMessages(library(rhdf5))
suppressMessages(library(vegan))
suppressMessages(library(Hmisc))
suppressMessages(library(igraph,quietly=T))
suppressMessages(library(plyr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(compositions))

# read options
config<-read.table("config",h=T,colClasses=c(cksum="character",mat="character"))
if ( grepl("\\.rds$",config$mat) ) {
  mat<-readRDS(config$mat)
} else {
  mat<-read.table(config$mat,h=T)
}
write(nrow(mat),"nbsamp_ori")
write(ncol(mat),"nbotu_ori")
for (i in 4:ncol(config)){assign(names(config)[i],config[1,i])}
if ( ncol(config) > 10 ) {
  sec=T
  if ( grepl("\\.rds$",config$second) ) {
    second<-readRDS(config$second)
  } else {
    second<-read.table(config$second,h=T)
  }
  if ( nrow(mat) != nrow(second) ) {
    quit(save="no",status=5,runLast=F)
  }
  write(ncol(second),"nbsec_ori")
} else {
  sec=F
}

# Drop OTUs with low occurrence
if(minocc>1) {
  mino<-minocc
} else {
  mino<-ceiling(nrow(mat)*minocc)
}
mat_tmp<-mat[,which(colSums((mat>0)*1)>=mino)]
# same filter on second matrix
if ( sec & same ) {
  if(minocc>1) {
    minos<-minocc
  } else {
    minos<-ceiling(nrow(second)*minocc)
  }
  second_tmp<-second[,which(colSums((second>0)*1)>=minos)]
}

# Drop samples with too few reads
if(mincount>1) {
  minc<-mincount
} else {
  minc<-ceiling(median(rowSums(mat_tmp))*mincount)
}
mat_ab<-mat_tmp[which(rowSums(mat_tmp)>=minc),]
# keep same samples on second matrix
if ( sec ) {
  if (same ) {
    second_ab<-second_tmp[rownames(mat_ab),]
  } else {
    second_ab<-second[rownames(mat_ab),]
  }
}

# Normalize read counts by multiplying each sample read counts by the ratio of half the read count median of all samples divide by total sample read counts
# better than rarefaction (but see McMurdie & Holmes, 2014)
# additional log, centered-log or sqrt transformation of ratio is to reduce the exponential abundance bias due to PCR
if(norm=="no") {
  mat_norm<-mat_ab
} else {
  if(norm=="ratio") {
    mat_tmp<-mat_ab
  } else if(norm=="ratio_log") {
    mat_tmp<-decostand(mat_ab,"log")
  } else if(norm=="ratio_sqrt") {
    mat_tmp<-decostand(mat_ab,"hellinger")
  } else if(norm=="clr") {
    mat_clr<-clr(mat_ab)
    mat_tmp<-t(apply(mat_clr,1,function(x){(x-min(x)*1.01)*((x!=0)*1)})) # shift all transformed values to positive, except zeroes
  }
  if(depth>1) {
    mat_norm<-round(mat_tmp*(depth/rowSums(mat_tmp)))
  } else {
    mat_norm<-round(mat_tmp*(median(rowSums(mat_ab))*depth/rowSums(mat_tmp)))
  }
}
# same for second
if ( sec & same ) {
  if(norm=="no") {
    second_norm<-second_ab
  } else {
    if(norm=="ratio") {
      second_tmp<-second_ab
    } else if(norm=="ratio_log") {
      second_tmp<-decostand(second_ab,"log")
    } else if(norm=="ratio_sqrt") {
      second_tmp<-decostand(second_ab,"hellinger")
    } else if(norm=="clr") {
      second_clr<-clr(second_ab)
      second_tmp<-t(apply(second_clr,1,function(x){(x-min(x)*1.01)*((x!=0)*1)})) # shift all transformed values to positive, except zeroes
    }
    if(depth>1) {
      second_norm<-round(second_tmp*(depth/rowSums(second_tmp)))
    } else {
      second_norm<-round(second_tmp*(median(rowSums(second_ab))*depth/rowSums(second_tmp)))
    }
  }
}

# save normalized OTU table and its size
saveRDS(as.matrix(mat_norm),"mat")
write(nrow(mat_norm),"nbsamp")
write(ncol(mat_norm),"nbotu")
otus<-colnames(mat_norm)
write(mino,"minocc")
write(minc,"mincount")

# save second matrix
if ( sec ) {
  if( same ) {
    saveRDS(as.matrix(second_norm),"second")
    write(ncol(second_norm),"nbsec")
    write(minos,"minoccsec")
  } else {
    saveRDS(as.matrix(second_ab),"second")
    write(ncol(second_ab),"nbsec")
  }
}

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
	write(c(otus,colnames(env_ab)),"otus",ncolumns=1)
} else {
	write(otus,"otus",ncolumns=1)
}

