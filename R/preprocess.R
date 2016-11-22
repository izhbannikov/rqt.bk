#### Dimensionality reduction methods ####

#' Calculates variance-covariance matrix for LASSO/ridge regression.
#' 
vcov_rigde <- function(x, y,  rmod) {
  
  ridge_se <- function(xs,y,yhat,my_mod){
    # Note, you can't estimate an intercept here
    n <- dim(xs)[1]
    k <- dim(xs)[2]
    sigma_sq <- sum((y - yhat)^2)/ (n-k)
    lam <- my_mod$lambda.min
    if(is.null(my_mod$lambda.min)==TRUE){lam <- 0}
    i_lams <- Matrix(diag(x=1,nrow=k,ncol=k),sparse=TRUE)
    xpx <- t(xs) %*% xs
    xpxinvplam <- solve(xpx + lam*i_lams)
    var_cov <- sigma_sq * (xpxinvplam %*% xpx %*% xpxinvplam)
    se_bs <- sqrt(diag(var_cov))
    
    print('NOTE: These standard errors are very biased.')
    return(list(vcov=var_cov, se=se_bs))
  }
  
  # Predictions
  r_yhat   <- predict(rmod,newx=x,s='lambda.min')
  ro_yhat  <- predict(rmod,newx=x)
  # Variance-covariance matrix and Standard Erros
  rmod_ses <- ridge_se(x,y,r_yhat,rmod)
  
  return(rmod_ses)
}

#### Calculating PCA ####
prerocess.pca <- function(data, scale, cumvar.threshold) {
  res.pca <- prcomp(data, scale=scale)
  # Eigenvalues
  eig <- (res.pca$sdev)^2
  # Variances in percentage
  variance <- eig*100/sum(eig)
  # Cumulative variances
  cumvar <- cumsum(variance)
  eig.decathlon2.active <- data.frame(eig = eig, variance = variance, cumvariance = cumvar)
  
  ######### Filtering by threshold ##############
  S <- res.pca$x[,which(eig.decathlon2.active$cumvar <= cumvar.threshold)] %*% t(res.pca$rotation[,which(eig.decathlon2.active$cumvar <= cumvar.threshold)])
  
  ########## And add the center (and re-scale) back to data ###########
  if(scale != FALSE){
    S <- scale(S, center = FALSE , scale=1/res.pca$scale)
  }
  if(center != FALSE){
    S <- scale(S, center = -1 * res.pca$center, scale=FALSE)
  }
  
  return(list(S = S))
}

## PLS-DA ##
preprocess.plsda <- function(data, phenotype) {
  numcomp <- ifelse(dim(data)[2] < 10, dim(data)[2], NA)
  res.plsda <- opls(x = data, y=as.factor(phenotype), predI=numcomp, plotL = FALSE, log10L=F, algoC = "nipals")
  #d.plsda <- cbind(y=phenotype, res.plsda$scoreMN)
  return(list(S=res.plsda$scoreMN))
}

## PLS ##
preprocess.pls <- function(data, phenotype, cumvar.threshold) {

  inpdata <- data.frame("y"=phenotype, data)
  res.pls.tmp <- plsr(y ~ ., data = inpdata, validation = "LOO")
  numcomp <- 1
  for(i in 1:res.pls.tmp$ncomp) {
    if (sum(res.pls.tmp$Xvar[1:i])/res.pls.tmp$Xtotvar > cumvar.threshold/100) {
      numcomp <- i
      break
    } else if(i == res.pls.tmp$ncomp) {
      numcomp <- res.pls.tmp$ncomp
    }
  }
  
  res.pls <- plsr(y ~ ., ncomp=numcomp, data = inpdata, validation = "LOO")
  #d.plsr <- cbind(y=inpdata$y, res.pls$scores)
  return(list(S=data.frame(matrix(res.pls$scores, ncol = dim(res.pls$scores)[2], nrow=dim(res.pls$scores)[1], byrow=TRUE))))
}

preprocess.lasso.ridge <- function(data, phenotype, reg.family, method) {
  #### LASSO/Ridge ####
  tryCatch({
    if(reg.family == "binomial") {
      fit <- cv.glmnet(x=as.matrix(data),alpha=ifelse(method=="lasso", 1, 0), y=as.factor(phenotype), family=reg.family)
    } else {
      fit <- cv.glmnet(x=as.matrix(data),alpha=ifelse(method=="lasso", 1, 0), y=phenotype, family=reg.family)
    }
  } , error=function(e) {
    stop(print(e))
  })
  
  return(list(fit=fit))
} 

simple.multvar.reg <- function(phenotype, data, reg.family) {
  if(reg.family == "binomial") {
    fit <- try(glm(phenotype ~ ., data=data.frame(data), family = binomial(link=logit)),TRUE)
  } else if(reg.family == "gaussian") {
    fit <- try(glm(phenotype ~ ., data=data.frame(data), family = reg.family),TRUE)
  } else {
    stop(paste("Unknown reg.family:", reg.family))
  }
  na.S <- try(which(is.na(coef(fit)[-1]) == TRUE),TRUE)
  
  if(length(na.S) > 0){
    S <- try(as.matrix(data[,-na.S]),TRUE)
    fit <- try(glm(phenotype ~ . ,data=data.frame(S)),TRUE)
  }
  return(list(S=S, fit=fit))
}
