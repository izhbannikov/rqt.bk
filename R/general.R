#' @import plyr
#' @import CCP
#' @import Matrix
#' @import CompQuadForm
#' @import parallel
#' @import pls
#' @import glmnet
#' @import ropls
#' @import methods
#' @import metap
#' @import SummarizedExperiment
#' @import car
#' @importFrom stats binomial cor glm
#' @importFrom stats pchisq pgamma prcomp 
#' @importFrom stats qbeta qchisq qgamma
#' @importFrom stats resid var vcov
#' @importFrom stats na.exclude na.omit
NULL

#' Empirical null distribution for Q3 test.
#' 
#' @return None
rQTtest.prepare <- function() {
  assign("null.dist.Q3", null.dist.Q3, envir=baseenv())
}

