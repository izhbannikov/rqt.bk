#### 

source("/Users/ilya/Projects/rqt/rqt.R")

rQTest <- function(pheno, geno, perm=0, STT=0.2,weight=FALSE, cumvar.threshold=90, reg.family="binomial") {
  if(perm==0){
    rslt<-try(QTest.one(y=pheno,covadat=NULL,newgeno=geno, STT=STT, cumvar.threshold=cumvar.threshold, reg.family=reg.family)$p.value,TRUE)
  } else {
    rsltPP <- do.call(rbind,mclapply(1:perm, function(k){
      yP <- pheno[sample(1:length(pheno),length(pheno),replace=F)]
      QTest.one(y=yP,covadat=NULL, newgeno=geno,STT=STT,weight=weight, cumvar.threshold=cumvar.threshold, reg.family=reg.family)$p.value
    },mc.cores=4))
      rslt0<-try(as.numeric(QTest.one(y=pheno, covadat=NULL, newgeno=geno, STT=STT, weight=weight, cumvar.threshold=cumvar.threshold, reg.family=reg.family)$p.value),TRUE)
      rslt <- data.frame(p.Q1 = (length(which(rsltPP[,1] < rslt0[1]))+1)/(perm+1),
                         p.Q2 = (length(which(rsltPP[,2] < rslt0[2]))+1)/(perm+1),
                         p.Q3 = (length(which(rsltPP[,3] < rslt0[3]))+1)/(perm+1))
  }
  
  rslt
}

data <- read.table("/Users/ilya/Projects/rqt/data/phengen2.dat", header=F, skip=1)
pheno <- data[,1]
geno <- data[,-1]

rQTest(pheno=pheno, geno=geno, perm=10000)
