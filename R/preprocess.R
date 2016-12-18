#### Dimensionality reduction methods ####

#' Calculates variance-covariance matrix for LASSO/ridge regression.
#'@param x A matrix of predictors.
#'@param y A vector of outputs (dependent variable).
#'@param rmod An object returned from glmnet.
#'@return A list of two: variance-covariance matrix, 
#'standard deviations of coefficients.
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

#' Preprocess input data with Principal Component Analysis method (PCA)
#' @param data An input matrix with values of 
#' independent variables (predictors).
#' @param scale A logical variable, indicates wheither or 
#' not scaling should be performed.
#' @param cumvar.threshold A threshold value for explained variance.
#' @return A list of one: "S" - a data frame of predictor values.
prerocess.pca <- function(data, scale, cumvar.threshold) {
    res.pca <- prcomp(data, scale=scale)
    # Eigenvalues
    eig <- (res.pca$sdev)^2
    # Variances in percentage
    variance <- eig*100/sum(eig)
    # Cumulative variances
    cumvar <- cumsum(variance)
    eig.decathlon2.active <- data.frame(eig = eig, 
        variance = variance, 
        cumvariance = cumvar)
  
    ######### Filtering by threshold ##############
    if(length(which(eig.decathlon2.active$cumvar <= cumvar.threshold)) == 0) {
        print("Warning: cumvar.threshold is too low and will be set to")
        print(paste("first component of cumulative variance:", 
                    eig.decathlon2.active$cumvar[1]))
        cumvar.threshold <- eig.decathlon2.active$cumvar[1]
    }
    
    S <- res.pca$x[,which(eig.decathlon2.active$cumvar <= 
        cumvar.threshold)] %*% 
    t(res.pca$rotation[,which(eig.decathlon2.active$cumvar <= 
        cumvar.threshold)])
  
    ########## And add the center (and re-scale) back to data ###########
    if(scale){
        S <- scale(S, center = FALSE , scale=1/res.pca$scale)
    }
    #if(center){
    #  S <- scale(S, center = -1 * res.pca$center, scale=FALSE)
    #}
  
    return(list(S = S))
}

#' Preprocess input data with Partial Linear Square 
#' Regregression Discriminant Analysis method (PLSDA)
#' @param data An input matrix with values of independent 
#' variables (predictors).
#' @param y A vector with values of dependent variable (outcome).
#' @return A list of one: "S" - a data frame of predictor values.
preprocess.plsda <- function(data, y) {
    numcomp <- ifelse(dim(data)[2] < 10, dim(data)[2], NA)
    model <- try(opls(x = data, y=as.factor(y), predI=numcomp, 
        plotL = FALSE, log10L=FALSE, algoC = "nipals"), 
        silent = TRUE)
    if(inherits(model, "try-error") &&
        substr(unclass(attr(model, "condition"))$message, 1, 85) == 
"No model was built because the first predictive component was already not significant") {
        model <- opls(x = data, y=as.factor(y), predI=1, plotL = FALSE, 
            log10L=FALSE, algoC = "nipals")
    }
    return(list(S=model@scoreMN))
}

#' Preprocess input data with Partial Linear Square Regregression 
#' method (PLS)
#' @param data An input matrix with values of 
#' independent variables (predictors).
#' @param y A vector with values of dependent variable (outcome).
#' @param cumvar.threshold A threshold value for explained variance.
#' @return A list of one: "S" - a data frame of predictors.
preprocess.pls <- function(data, y, cumvar.threshold) {
    inpdata <- data.frame("y"=y, data)
    res.pls.tmp <- plsr(y ~ ., data = inpdata, validation = "LOO")
    numcomp <- 1
    for(i in 1:res.pls.tmp$ncomp) {
        if(sum(res.pls.tmp$Xvar[1:i])/res.pls.tmp$Xtotvar > 
            cumvar.threshold/100) {
            numcomp <- i
            break
        } else if(i == res.pls.tmp$ncomp) {
            numcomp <- res.pls.tmp$ncomp
        }
    }
  
    res.pls <- plsr(y ~ ., ncomp=numcomp, 
        data = inpdata, validation = "LOO")
    #d.plsr <- cbind(y=inpdata$y, res.pls$scores)
    S <- data.frame(matrix(res.pls$scores, ncol = dim(res.pls$scores)[2],
        nrow=dim(res.pls$scores)[1], byrow=TRUE))
  
    return(list(S=S))
}

#' Preprocess input data with LASSO/Ridge regregression method
#' @param data An input matrix with values of 
#' independent variables (predictors).
#' @param y A vector with values of dependent variable (outcome).
#' @param reg.family A regression family. 
#' Can be either "binomial" or "gaussian."
#' @param method A method. Can be either "lasso" or "ridge."
#' @return fit An object returned by "cv.glmnet" function.
preprocess.lasso.ridge <- function(data, y, reg.family, method) {
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
        stop(print(e))
    })
    return(list(fit=fit))
} 

#' Applies linear of logistic regregression to the data.
#' @param data An input matrix with values of 
#' independent variables (predictors).
#' @param y A vector with values of dependent variable (outcome).
#' @param reg.family A regression family. 
#' Can be either "binomial" or "gaussian."
#' @return A list of two: "S" - a dataframe with predictors and "fit"
#' - an object returned by "glm" function.
simple.multvar.reg <- function(y, data, reg.family) {
    if(reg.family == "binomial") {
        fit <- try(glm(y ~ ., data=data.frame(data), 
            family = binomial(link=logit)),TRUE)
    } else if(reg.family == "gaussian") {
        fit <- try(glm(y ~ ., data=data.frame(data), 
            family = reg.family),TRUE)
    } else {
        stop(paste("Unknown reg.family:", reg.family))
    }
    na.S <- try(which(is.na(coef(fit)[-1]) == TRUE),TRUE)

    if(length(na.S) > 0){
        S <- try(as.matrix(data[,-na.S]),TRUE)
        fit <- try(glm(y ~ . ,data=data.frame(S)),TRUE)
    } else {
        S <- data
    }
    return(list(S=S, fit=fit))
}
