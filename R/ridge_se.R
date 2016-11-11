require(glmnet)
require(ridge)
set.seed(5)
n <- 1e4
k <- 10
x <- matrix(rnorm(n*k),ncol=k)
beta <- rpois(k,lambda=4)/10
print(beta)
yhat <- x%*%beta 
y <- 1*yhat+1*rnorm(n)
#========================================================
ridge_se <- function(xs,y,yhat,my_mod){
  # Note, you can't estimate an intercept here
  n <- dim(xs)[1]
  k <- dim(xs)[2]
  sigma_sq <- sum((y-yhat)^2)/ (n-k)
  lam <- my_mod$lambda.min
  if(is.null(my_mod$lambda.min)==TRUE){lam <- 0}
  i_lams <- Matrix(diag(x=1,nrow=k,ncol=k),sparse=TRUE)
  xpx <- t(xs)%*%xs
  xpxinvplam <- solve(xpx+lam*i_lams)
  var_cov <- sigma_sq * (xpxinvplam %*% xpx %*% xpxinvplam)
  se_bs <- sqrt(diag(var_cov))
  print(var_cov)
  print('NOTE: These standard errors are very biased.')
  return(se_bs)
}
#========================================================
lmod <- glmnet(x=x,y=y,intercept=FALSE,
               lambda=0,alpha=0,standardize=TRUE)
rmod <- cv.glmnet(x=x,y=y,
                  nlambda=100,alpha=0,nfolds=5,
                  ,standardize=TRUE,intercept=FALSE)
print(rmod$lambda.min)
#========================================================
lmod2 <-summary( lm(y~scale(x)-1))
ridge <- linearRidge(y~x-1,scaling='scale',lambda=rmod$lambda.min)
#========================================================
# Predictions
r_yhat   <- predict(rmod,newx=x,s='lambda.min')
l_yhat   <- predict(lmod,newx=x,s='lambda.min')
ro_yhat  <- predict(rmod,newx=x)
# Standard Erros
rmod_ses        <- ridge_se(x,y,r_yhat,rmod)
lmod_ses_GLMNET <- ridge_se(x,y,l_yhat,lmod)
romod_ses       <- data.frame(summary(ridge)$summaries$summary1[1])
names(romod_ses)<- c('Estimate','Scale_Estimate','Std_Error','T_value','P_Value')
lmod_ses        <- data.frame(lmod2[4])
names(lmod_ses) <- c('Estimate','Std_Error','T_value','P_Value')
lmod_ses <- data.frame(Names=row.names(lmod_ses),
                       OLS=lmod_ses$Std_Error)
#========================================================
std_errs <- lmod_ses
std_errs$Ridge <- rmod_ses
std_errs$OLS_GLMNET <- lmod_ses_GLMNET
std_errs$OLS_Diffs <- with(std_errs,(OLS-OLS_GLMNET))
std_errs$Orig_Ridge <- romod_ses$Std_Error
std_errs$Ridge_Diff <- with(std_errs,Ridge-Orig_Ridge)
# It looks like GLMNET is giving the accurate standard errors, at least within 4 decimal places
print(std_errs)