# A basic function for gene-level association test #

library(plyr)
library(CCP)
library(homals)
library(Matrix)
library(CompQuadForm)
library(parallel)

## Empirical null distribution for Q3 test
null.dist.Q3 <- read.table("~/Projects/rqt/n.log10.minp.1e09.txt",header=T)

## Get a given STT
get.a<-function(L,STT){
  aa<-diff<-seq(0,1,length=200)
  for(i in 1:length(aa)){
    diff[i]<-abs(get.stt(L,aa[i],STT)-STT)
  }
  return(aa[which.min(diff)])
}
get.stt<-function(L,a,STT){
  1-pgamma(L*qgamma(1-STT,a,1),L*a,1)
}
## get minor allele frequency
get.maf <-function(vec){
  (length(which(vec==1))+2*length(which(vec==2)))/(2*length(vec))
}




###############################################################################################
## QTest.one(phe.cova,geno,yname, r2)  #QTest for single gene                                
##                                                                                           
## y: trait, newgeno: genotype matrix                                           
## cut.r2: r^2 cutoff, n.perm: # of pumutations                                              
## n.perm: number of permutation for GM method                                              
## a: shape parameter for GM method                                                          
###############################################################################################

QTest.one<-function(y,covadat=NULL,newgeno,STT=0.2,weight=FALSE, cumvar.threshold=90, eig.threshold=0.0005){
  ### DEBUG ###
  y <- matrix(ncol=1,data[,1])
  newgeno <- data[,2:dim(data)[2]]
  covadat <- NULL
  STT <- 0.2
  weight <- FALSE
  
  ####
  ## If covariates exist
  #if(length(covadat)!=0){
  #  resid<-try(resid(glm(y~.,data=data.frame(covadat),na.action=na.exclude)),TRUE)
  #} else {
  #  resid<-try(resid(glm(y~1,na.action=na.exclude, family = binomial(link = log))),TRUE)
  #}
  
  #### Calculating PCA ####
  pcadata <- cbind(y, newgeno)
  res.pca <- prcomp(newgeno)
  # Eigenvalues
  eig <- (res.pca$sdev)^2
  # Variances in percentage
  variance <- eig*100/sum(eig)
  # Cumulative variances
  cumvar <- cumsum(variance)
  eig.decathlon2.active <- data.frame(eig = eig, variance = variance, cumvariance = cumvar)
  #head(eig.decathlon2.active)
  #### End of calculating PCA ####
  
  #### Regression after PCA ####
  S <- try(as.matrix(res.pca$x[,which(eig.decathlon2.active$cumvar <= cumvar.threshold)]),TRUE)
  t.fit <- try(glm(y ~ .,data=data.frame(S), family = poisson(link=log)),TRUE)
  fit <- try(glm(y ~ .,data=data.frame(S), family = binomial(link = log), start=coef(t.fit)),TRUE)
  na.S <- try(which(is.na(coef(fit)[-1]) == TRUE),TRUE)
  
  #### Calculating statistics and p-values ####
  if(length(na.S)>0){
    S <- try(as.matrix(S[,-na.S]),TRUE)
    fit <- try(glm(resid~.,data=data.frame(S)),TRUE)
  }
  coef <- try(coef(summary(fit))[-1,1:2],TRUE)
  
  
  if(mode(fit)=="character"){length(coef)<-0}
  
  if(length(coef)!=0){
    if(length(coef)!=2){beta1<-coef[,1];se1<-coef[,2]}
    if(length(coef)==2){beta1<-coef[1];se1<-coef[2]}
    SS<-cbind(1,S); n<-length(y)
    vv<-vcov(fit)[-1,-1]
    alpha<-(1/(se1^2))/sum(1/(se1^2))
    #vv<-solve(t(SS)%*%SS)[-1,-1]*var(y)*(n-1)/n;alpha<-(1/(se1^2))
    
    ##QTest1##
    if(weight==FALSE){
      var.pool0<-t(alpha)%*%vv%*%alpha
      beta.pool0<-t(alpha)%*%beta1
      z.score0<-beta.pool0/sqrt(var.pool0)
      Q1<-z.score0^2
      p.Q1<-pchisq(Q1,df=1,lower.tail=FALSE)
    }
    if(weight==TRUE){
      maf.S<-apply(S,2,function(v)mean(na.omit(v))/2)
      #w.S0<-1/sqrt(maf.S*(1-maf.S));
      w.S0<-qbeta(maf.S,1,25,lower.tail=F)
      WS<-diag(w.S0);if(length(beta1)==1){WS<-w.S0}
      var.pool<-t(alpha)%*%WS%*%vv%*%WS%*%alpha
      beta.pool<-t(alpha)%*%WS%*%beta1
      z.score<-beta.pool/sqrt(var.pool)
      Q1<-z.score^2
      p.Q1<-pchisq(Q1,df=1,lower.tail=FALSE)
    }
    
    ## QTest2 ##
    Q2.eigen<-eigen(vv);U2<-Q2.eigen$vectors; l2<-Q2.eigen$values;na.l2<-which(l2/mean(l2)<0.01)
    if(length(na.l2)>0){l2<-l2[-na.l2];U2<-U2[,-na.l2]}
    p2<-pchisq((t(U2)%*%beta1)^2/l2,df=1,lower.tail=FALSE)
    a<-get.a(length(p2),STT)
    q2<-2*(qgamma(p2,a,1,lower.tail=FALSE))
    Q2<-sum(q2)
    p.Q2<-pchisq(Q2,df=2*a*length(q2),lower.tail=FALSE)
    
    
    ##QTest3 ##
    
    if(weight==FALSE){
      c(t(alpha)%*%vv)->cov.beta
      b.star<-beta1-beta.pool0*cov.beta/var.pool0[1]
      vv.star<-vv-(cov.beta%*%t(cov.beta))/var.pool0[1]
    }
    if(weight==TRUE){
      w.vv<-WS%*%vv%*%WS
      c(t(alpha)%*%w.vv)->cov.beta
      b.star<-WS%*%beta1-beta.pool*cov.beta/var.pool[1]
      vv.star<-w.vv-(cov.beta%*%t(cov.beta))/var.pool[1]
    }
    
    if(length(beta1)!=1){
      Q3.eigen<-eigen(vv.star); U3<-Q3.eigen$vectors; l3<-Q3.eigen$values;na.l3<-which(l3/mean(l3)<0.001)
      if(length(na.l3)>0){l3<-l3[-na.l3];U3<-U3[,-na.l3]}
      L3<-try(diag(l3),TRUE);if(length(l3)==1){L3<-l3}
      q2.proj<-(t(U3)%*%b.star)^2/l3; 
      p2.1<-pchisq(q2.proj,df=1,lower.tail=FALSE)
      a<-get.a(length(p2.1),STT)
      Q2.proj<-sum(2*qgamma(p2.1,a,1,lower.tail=FALSE))
      p.Q2.proj<-pchisq(Q2.proj,df=2*a*length(l3),lower.tail=FALSE)
      Q2.1<-qchisq(p.Q2.proj,df=1,lower.tail=FALSE)
      
      pi0<-seq(0,1,by=0.1);p.Q3.can<-Q3<-rep(1,11)
      p.Q3.can[1]<-p.Q2.proj; Q3[1] <- Q2.1
      p.Q3.can[11]<-p.Q1; Q3[11]<-Q1
      for(h in 2:10){
        Q3[h]<-pi0[h]*Q1+(1-pi0[h])*Q2.1
        p.Q3.can[h]<-davies(Q3[h],c(pi0[h],(1-pi0[h])),c(1,1))$Qq
        if(p.Q3.can[h]<=0|p.Q3.can[h]>1){p.Q3.can[h]<-imhof(Q3[h],c(pi0[h],(1-pi0[h])),c(1,1))$Qq}
        if(p.Q3.can[h]<=0|p.Q3.can[h]>1){p.Q3.can[h]<-liu(Q3[h],c(pi0[h],(1-pi0[h])),c(1,1))[1]}
      }
      
      
      Q3final<-Q3[which.min(p.Q3.can)]
      p.Q3<-(sum(null.dist.Q3[null.dist.Q3[,1]>-log10(min(p.Q3.can)),2])+1)/(sum(null.dist.Q3[,2])+1)
      
    }
    
    if(length(beta1)==1){p.Q3<-p.Q1; Q3final<-Q1}
    rslt<-list( data.frame(Q1, Q2, Q3=Q3final), data.frame(p.Q1,p.Q2,p.Q3))
    names(rslt)<-c("Qstatistic", "p.value")	
  }
  
  if(length(coef)==0){rslt<-NA}
  options (warn=-1)
  return(rslt)	
}

data <- read.table("~/Projects/rqt/data/phengen2.dat", skip=1)
pheno <- matrix(ncol=1,data[,1])
geno <- data[,2:dim(data)[2]]

QTest.one(y=pheno,newgeno=geno)
