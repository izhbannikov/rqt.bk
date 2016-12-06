# A basic function for gene-level association test #


#' Get a given STT
#' @param L TODO
#' @param STT Numeric indicating soft truncation threshold (STT) to convert 
#' to gamma parameter (must be <= 0.4).
#' @return a TODO
get.a<-function(L,STT) {
    aa <- diff <- seq(0,1,length=200)
    for(i in 1:length(aa)){
      diff[i] <- abs(get.stt(L,aa[i], STT) - STT)
    }
    return(aa[which.min(diff)])
}

get.stt<-function(L,a,STT){
    1-pgamma(L*qgamma(1-STT,a,1),L*a,1)
}


get.reg.family <- function(out.type) {
    if(out.type == "D") {
      reg.family="binomial"
    } else if(out.type == "C") {
      reg.family="gaussian"
    } else {
      stop(paste("Unknown out.type:", out.type))
    }
    return(reg.family)
}

QTest.one <- function(phenotype, genotype, covariates, 
                      STT=0.2,weight=FALSE, 
                      cumvar.threshold=90, 
                      method="pca", out.type="D", 
                      scale=FALSE) {
    ### Data preprocessing, (scaling if needed) ###
    #### Binding predictors (genotype and covariates) ####
    preddata <- genotype
    
    if(length(covariates) != 0) {
      tryCatch({
        preddata <- cbind(genotype, covariates)
      }, error=function(e) {
        stop(print(e))
      })
    }
  
    if(scale) {
      preddata <- as.data.frame(scale(preddata))
    }
  
    reg.family <- get.reg.family(out.type)
  
    rslt <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                data.frame(p.Q1=NA, p.Q2=NA, p.Q3=NA) )
    names(rslt)<-c("Qstatistic", "p.value")
    res <- list()
    tryCatch({
      if(dim(preddata)[2] > 1) {
        ### Dimensionality reduction and account for LD ###
        if(method == "pca") {
          res <- prerocess.pca(data=preddata, scale=scale, 
                             cumvar.threshold=cumvar.threshold)
        } else if(method == "pls") {
          if(out.type == "D") {
            ##### PLSDA #####
            res <- preprocess.plsda(data=preddata, y=phenotype)
          } else if(out.type == "C"){
            ##### PLS #####
            res <- preprocess.pls(data=preddata, y=phenotype, 
                                cumvar.threshold=cumvar.threshold)
          }
        } else if(method == "lasso" | method == "ridge") {
          res <- preprocess.lasso.ridge(data=preddata, y=phenotype, 
                                      reg.family=reg.family, method=method)
        } else {
          res[["S"]] <- preddata
        }
        
        #### Regression after data preprocessing ####
        if(!(method %in% c("lasso", "ridge"))) {
          S <- res[["S"]]
          res <- simple.multvar.reg(y=phenotype, data=S, reg.family=reg.family)
          fit <- res$fit
          coef <- try(coef(summary(fit))[-1,1:2],TRUE)
          if(mode(fit)=="character"){length(coef) <-0 }
          if(length(coef)!=2){beta1<-coef[,1];se1 <- coef[,2]}
          if(length(coef)==2){beta1<-coef[1];se1 <- coef[2]}
          
          vv <- vcov(fit)[-1,-1]
          alpha <- (1/(se1^2)) #/sum(1/(se1^2))
        } else {
          fit <- res$fit
          coef <- coef(fit)[-1]
          
          if(sum(coef) == 0) {
              print("All coefficients in lasso/ridge regression 
                    are equal to 0. Trying ordinary regressing instead.")
              res <- simple.multvar.reg(y=phenotype, data=preddata, 
                                    reg.family=reg.family)
              S <- res$S
              fit <- res$fit
              coef <- try(coef(summary(fit))[-1,1:2],TRUE)
          
              if(mode(fit)=="character"){length(coef) <-0 }
              if(length(coef)!=2){beta1<-coef[,1];se1 <- coef[,2]}
              if(length(coef)==2){beta1<-coef[1];se1 <- coef[2]}
              vv <- vcov(fit)[-1,-1]
          } else {
            beta1 <- coef[which(coef != 0)]
            vv <- vcov_rigde(x=as.matrix(preddata), 
                           y=phenotype, 
                           rmod=fit)$vcov[which(coef != 0), which(coef != 0)]
            if(class(vv)[1] == "dgeMatrix") {
              se1 <- sqrt(diag(vv))
            } else {
              se1 <- sqrt(vv)
            }
            vv <- as.matrix(vv)
          }
          
          alpha <- as.matrix((1/(se1^2)), ncol=1) #/sum(1/(se1^2))
        }
        
        if(length(coef) != 0) {
          ###### QTest1 ######
          if(weight == FALSE) {
            var.pool0 <- t(alpha) %*% vv %*% alpha
            beta.pool0 <- t(alpha) %*% beta1
            z.score0 <- beta.pool0/sqrt(var.pool0)
            Q1 <- z.score0^2
            p.Q1 <- pchisq(Q1, df=1, lower.tail=FALSE)
          }
          
          if(weight == TRUE) {
            maf.S <- apply(S,2,function(v)mean(na.omit(v))/2)
            w.S0 <- qbeta(maf.S,1,25,lower.tail=FALSE)
            WS <- diag(w.S0)
            if(length(beta1)==1){
              WS <- w.S0
            }
            var.pool <- t(alpha) %*% WS %*% vv %*% WS %*% alpha
            beta.pool <- t(alpha) %*% WS %*% beta1
            z.score <- beta.pool/sqrt(var.pool)
            Q1 <- z.score^2
            p.Q1 <- pchisq(Q1, df=1, lower.tail=FALSE)
          }
          
          ## QTest2 ##
          Q2.eigen <- eigen(vv)
          U2 <- Q2.eigen$vectors
          l2 <- Q2.eigen$values
          na.l2 <- which(l2/mean(l2) < 0.01)
          if(length(na.l2) > 0){
            l2 <- l2[-na.l2]
            U2 <- U2[,-na.l2]
          }
          
          p2 <- pchisq((t(U2) %*% beta1)^2/l2, df=1, lower.tail=FALSE)
          a <- get.a(length(p2), STT)
          q2 <- 2*(qgamma(p2, a, 1, lower.tail=FALSE))
          Q2 <- sum(q2)
          p.Q2 <- pchisq(Q2, df=2*a*length(q2), lower.tail=FALSE)
          
          ##QTest3 ##
          if(weight==FALSE){
            cov.beta <- c(t(alpha) %*% vv)
            b.star <- beta1 - beta.pool0*cov.beta/var.pool0[1]
            vv.star <- vv - (cov.beta%*%t(cov.beta))/var.pool0[1]
          }
          if(weight==TRUE){
            w.vv<-WS%*%vv%*%WS
            c(t(alpha)%*%w.vv)->cov.beta
            b.star<-WS%*%beta1-beta.pool*cov.beta/var.pool[1]
            vv.star<-w.vv-(cov.beta%*%t(cov.beta))/var.pool[1]
          }
    
          if(length(beta1)!=1){
            Q3.eigen<-eigen(vv.star)
            U3<-Q3.eigen$vectors
            l3<-Q3.eigen$values
            na.l3<-which(l3/mean(l3)<0.001)
          
            if(length(na.l3)>0){l3<-l3[-na.l3];U3<-U3[,-na.l3]}
            L3<-try(diag(l3),TRUE);if(length(l3)==1){L3<-l3}
            q2.proj<-(t(U3)%*%b.star)^2/l3; 
            p2.1<-pchisq(q2.proj,df=1,lower.tail=FALSE)
            a <- get.a(length(p2.1),STT)
            Q2.proj<-sum(2*qgamma(p2.1,a,1,lower.tail=FALSE))
            p.Q2.proj<-pchisq(Q2.proj,df=2*a*length(l3),lower.tail=FALSE)
            Q2.1<-qchisq(p.Q2.proj,df=1,lower.tail=FALSE)
      
            pi0<-seq(0,1,by=0.1);p.Q3.can<-Q3<-rep(1,11)
            p.Q3.can[1]<-p.Q2.proj; Q3[1] <- Q2.1
            p.Q3.can[11]<-p.Q1; Q3[11]<-Q1
            for(h in 2:10){
              Q3[h]<-pi0[h]*Q1+(1-pi0[h])*Q2.1
              p.Q3.can[h]<-davies(Q3[h],c(pi0[h],(1-pi0[h])),c(1,1))$Qq
              if(p.Q3.can[h]<=0|p.Q3.can[h]>1){
                p.Q3.can[h]<-imhof(Q3[h],c(pi0[h],(1-pi0[h])),c(1,1))$Qq
              }
              if(p.Q3.can[h]<=0|p.Q3.can[h]>1){
                p.Q3.can[h]<-liu(Q3[h],c(pi0[h],(1-pi0[h])),c(1,1))[1]
              }
            }
      
            Q3final <- Q3[which.min(p.Q3.can)]
            p.Q3 <- (sum(
              null.dist.Q3[null.dist.Q3[,1] > 
                        -log10(min(p.Q3.can)),2])+1)/(sum(null.dist.Q3[,2])+1)
          }
    
          if(length(beta1)==1){p.Q3<-p.Q1; Q3final<-Q1}
          
          rslt <- list( data.frame(Q1, Q2, Q3=Q3final), 
                      data.frame(p.Q1,p.Q2,p.Q3))
          names(rslt)<-c("Qstatistic", "p.value")
        }
  
        if(length(coef)==0) { 
          rslt <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                      data.frame(p.Q1=1,p.Q2=1,p.Q3=1) )
        }
        options (warn=-1)
      } else {
        # Simple logistic regression:
        if(out.type == "D") {
          res <- glm(phenotype ~ ., data=data.frame(preddata), 
                   family = binomial(link=logit))
        } else {
          res <- glm(phenotype ~ ., data=data.frame(preddata), 
                   family = gaussian)
        }
        reg.coef <- coef(summary(res))
    
        if(dim(reg.coef)[1] == 2) {
          rslt<-list( data.frame(Q1=reg.coef[2,3], 
                               Q2=reg.coef[2,3], 
                               Q3=reg.coef[2,3]), 
                    data.frame(p.Q1=reg.coef[2,4],
                               p.Q2=reg.coef[2,4],
                               p.Q3=reg.coef[2,4]))
        } else if(dim(reg.coef)[1] == 1) {
          rslt <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                      data.frame(p.Q1=NA,p.Q2=NA,p.Q3=NA) )
        } else {
          rslt <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                      data.frame(p.Q1=NA,p.Q2=NA,p.Q3=NA) )
        }
      }
    },error=function(e) {
      print(e)
      rslt <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                  data.frame(p.Q1=NA,p.Q2=NA,p.Q3=NA) )
      names(rslt)<-c("Qstatistic", "p.value")
    }, finally=rslt)

    if(!is.na(rslt)) {
      names(rslt)<-c("Qstatistic", "p.value")
    } else {
      rslt <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                  data.frame(p.Q1=NA,p.Q2=NA,p.Q3=NA) )
      names(rslt)<-c("Qstatistic", "p.value")
    }
  
    return(rslt)
}




