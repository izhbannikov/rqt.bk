#### Dimensionality reduction methods ####

#' Calculates a standard error (SE) from glmnet.
#'@param xs A matrix of predictors.
#'@param y A vector of outputs (dependent variable).
#'@param yhat A vector of predicted values
#'from \code{glmnet.predict(...)}.
#'@param my_mod An object returned from glmnet.
#'@param verbose Indicates verbosing output,
#'FALSE by default.
#'@return A list of two: variance-covariance matrix, 
#'standard deviations of coefficients.
ridge_se <- function(xs,y,yhat,my_mod, verbose=FALSE){
  # Note, you can't estimate an intercept here
  n <- dim(xs)[1]
  k <- dim(xs)[2]
  sigma_sq <- sum((y - yhat)^2)/ (n-k)
  lam <- my_mod$lambda.min
  if(is.null(my_mod$lambda.min)) {
    lam <- 0
  }
  i_lams <- Matrix(diag(x=1,nrow=k,ncol=k),sparse=TRUE)
  xpx <- t(xs) %*% xs
  xpxinvplam <- solve(xpx + lam*i_lams)
  var_cov <- sigma_sq * (xpxinvplam %*% xpx %*% xpxinvplam)
  se_bs <- sqrt(diag(var_cov))
  
  if(verbose) 
    print('NOTE: These standard errors are very biased.')
  
  return(list(vcov=var_cov, se=se_bs))
}

#' Calculates variance-covariance matrix for LASSO/ridge regression.
#'@param x A matrix of predictors.
#'@param y A vector of outputs (dependent variable).
#'@param rmod An object returned from glmnet.
#'@param verbose Indicates verbosing output. Default: FALSE.
#'@return A list of two: variance-covariance matrix, 
#'standard deviations of coefficients.
vcov_rigde <- function(x, y,  rmod, verbose=FALSE) {
    # Predictions
    r_yhat   <- predict(rmod,newx=x,s='lambda.min')
    ro_yhat  <- predict(rmod,newx=x)
    # Variance-covariance matrix and Standard Erros
    rmod_ses <- ridge_se(x,y,r_yhat,rmod,verbose=verbose)
    return(rmod_ses)
}

#' Preprocess input data with Principal Component Analysis method (PCA)
#' @param data An input matrix with values of 
#' independent variables (predictors).
#' @param y A phenotype - column-vector, needed for LASSO/ridge and 
#' \code{NULL} by default.
#' @param method A dimensionality reduction method.
#' Default: \code{pca}.
#' @param reg.family A regression family. 
#' Default: \code{"binomial"}.
#' @param scaleData A logical variable, indicates wheither or 
#' not scaling should be performed. Default: \code{FALSE}.
#' @param cumvar.threshold A threshold value for explained variance.
#' Default: \code{75}
#' @param verbose Indicates verbosing output. Default: FALSE.
#' @param out.type An output (phenotype) type. 
#' Default: \code{"D"}.
#' @return A list of one: "S" - a data frame of predictor values.
preprocess <- function(data, y=NULL,
                       method="pca",
                       reg.family="binomial", 
                       scaleData=FALSE, 
                       cumvar.threshold=75,
                       out.type="D",
                       verbose=FALSE) {
    
    
    switch(method, 
        pca={
            ct <- cumvar.threshold
            res.pca <- prcomp(data, scale=scaleData)
            # Eigenvalues
            eig <- (res.pca$sdev)^2
            # Variances in percentage
            variance <- eig*100/sum(eig)
            # Cumulative variances
            cumvar <- cumsum(variance)
            eig.table <- data.frame(eig = eig, 
                                    variance = variance, 
                                    cumvariance = cumvar)
          
            ######### Filtering by threshold ##############
            if(length(eig.table$cumvar[eig.table$cumvar <= ct]) == 0) {
                if(verbose) {  
                    print("Warning: cumvar.threshold is too low and will be set to")
                    print(paste("first component of cumulative variance:", 
                          eig.table$cumvar[1]))
                }
                cumvar.threshold <- eig.table$cumvar[1]
            }
          
            #S <- res.pca$x[,which(eig.decathlon2.active$cumvar <= 
            #    cumvar.threshold)] %*% 
            #t(res.pca$rotation[,which(eig.decathlon2.active$cumvar <= 
            #    cumvar.threshold)])
          
            S <- res.pca$x[,which(eig.table$cumvar <= ct)]
          
            ########## And add the center (and re-scale) back to data ###########
            #if(scale){
            #   S <- scale(S, center = FALSE , scale=1/res.pca$scale)
            #}
            #if(center){
            #  S <- scale(S, center = -1 * res.pca$center, scale=FALSE)
            #}
            
            indexes <- which(eig.table$cumvar <= ct)
            return(list(S = S, indexes=indexes))
        },
        pls={
          
            ct <- cumvar.threshold/100
            
            if(scaleData){
                data.scaled <- scale(data, center = TRUE)
            } else {
                data.scaled <- data
            }
            
            npred <- round(dim(data.scaled)[2]*ct)
            numcomp <- ifelse(dim(data.scaled)[2] < 10, dim(data.scaled)[2], npred)
            
            if(out.type == "D") { # PLS-DA
                model <- try(opls(x = data.scaled, y=as.factor(y), predI=numcomp, 
                            plotL = FALSE, log10L=FALSE, algoC = "nipals"), 
                       silent = TRUE)
          
                if(inherits(model, "try-error")) {
                    for(i in 2:6) {
                        ct <- ct/i
                        npred <- round(dim(data.scaled)[2]*ct)
              
                        model <- opls(x = data.scaled, y=as.factor(y), predI=npred, 
                                plotL = FALSE, 
                                log10L=FALSE, algoC = "nipals", silent = TRUE)
              
                        if(!inherits(model, "try-error")) {
                            break
                        }
                    }
                }
            } else if(out.type == "C") { # PLS
                model <- try(opls(x = data.scaled, y=y, predI=numcomp, 
                                plotL = FALSE, log10L=FALSE, algoC = "nipals"), 
                           silent = TRUE)
              
                if(inherits(model, "try-error")) {
                    cumvar.threshold <- cumvar.threshold/2
                    npred <- round(dim(data.scaled)[2]*ct)
                    model <- opls(x = data.scaled, y=as.factor(y), 
                              predI=npred, plotL = FALSE, 
                              log10L=FALSE, algoC = "nipals", 
                              silent = TRUE)
                }
            }
          
            S <- model@scoreMN
            return(list(S = S))
        },
        lasso={
            #### LASSO/Ridge ####
            tryCatch({
                if(reg.family == "binomial") {
                  fit <- cv.glmnet(x=as.matrix(data),
                               alpha=ifelse(method=="lasso", 1, 0),
                               y=as.factor(y), family=reg.family)
                } else {
                  fit <- cv.glmnet(x=as.matrix(data),
                               alpha=ifelse(method=="lasso", 1, 0),
                               y=y, family=reg.family)
                }
            } , error=function(e) {
                print(e)
            })
          
            return(list(fit=fit))
        },
        ridge={
            tryCatch({
                if(reg.family == "binomial") {
                    fit <- cv.glmnet(x=as.matrix(data),
                             alpha=ifelse(method=="lasso", 1, 0),
                             y=as.factor(y), family=reg.family)
                } else {
                    fit <- cv.glmnet(x=as.matrix(data),
                             alpha=ifelse(method=="lasso", 1, 0),
                             y=y, family=reg.family)
                }
            } , error=function(e) {
                print(e)
            })
  
            return(list(fit=fit))
        },
        {
            stop("Unknown method provided.")
        }
  )
}


#' Applies linear of logistic regregression to the data.
#' @param null.model A fitted null model
#' @param Z A genotype matrix
#' @param verbose Indicates verbosing output. Default: FALSE.
#' @return A list of two: "S" - a dataframe with predictors and "fit"
#' - an object returned by "glm" function.
simple.multvar.reg <- function(null.model, Z, verbose=FALSE) {
    # Fit regression according to provided null model #
    fit <- glm(null.model ~ ., data=data.frame(Z))
    na.S <- which(is.na(coef(fit)[-1]) == TRUE)

    if( length(na.S) > 0 & (dim(data.frame(Z))[2] > 1) ){
        S <- as.matrix(Z[,-na.S])
        fit <- glm(null.model ~ . ,data=data.frame(S))
    } else {
        S <- Z
    }
    return(list(S=S, fit=fit, na.S=na.S))
}

#' Applies linear of logistic regregression to the data.
#' @param y A vector with values of dependent variable (outcome).
#' @param x A data.frame of covariates.
#' @param reg.family A regression family. 
#' Can be either "binomial" or "gaussian."
#' @param verbose Indicates verbosing output. Default: FALSE.
#' @return A list of two: "S" - a dataframe with predictors and "fit"
#' - an object returned by "glm" function.
build.null.model <- function(y, x, reg.family="binomial", verbose=FALSE) {
    if(reg.family == "binomial") {
        if(length(x) != 0) {
            fit <- glm(y ~ ., data=data.frame(x), 
                           na.action=na.exclude,
                         family = binomial) #family = binomial(link=logit)),TRUE)
        } else {
            fit <- glm(y ~ 1, na.action=na.exclude,
                         family = binomial) #family = binomial(link=logit)),TRUE)
        }
    } else if(reg.family == "gaussian") {
        if(length(x) != 0) {
            fit <- glm(y ~ x, data=data.frame(x), 
                   na.action=na.exclude,
                   family = reg.family)
        } else {
            fit <- glm(y ~ 1, 
                     na.action=na.exclude,
                     family = reg.family)
        }
    } else {
        stop(paste("Unknown reg.family:", reg.family))
    }
  
    return(resid(fit))
}



