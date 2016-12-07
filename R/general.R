#' @import plyr
#' @import CCP
#' @import Matrix
#' @import CompQuadForm
#' @import parallel
#' @import pls
#' @import glmnet
#' @import ropls
NULL

#' Empirical null distribution for Q3 test.
#' 
#' @return None
rQTtest.prepare <- function() {
  assign("null.dist.Q3", null.dist.Q3, envir=baseenv())
}

