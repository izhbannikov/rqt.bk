# Clases

#'The rqt class
#'
#'This class stores parameters and results of the rtq algorithms
#'
setClass("rqt", slots=c(phenotype="vector", genotype="data.frame", covariates="data.frame", results="list"))

#' The rqt class constructor
#' 
#' This function generates rqt class objects
#' 
rqtClass <- function(phenotype=NULL, genotype=NULL, covariates=NULL, results=NULL) {
  if(is.null(phenotype)) {
    phenotype <- c()
  }
  if(is.null(genotype)) {
    genotype <- data.frame()
  }
  if(is.null(covariates)) {
    covariates <- data.frame()
  }
  if(is.null(results)) {
    results <- list()
  }
    
  new("rqt", phenotype=phenotype, genotype=genotype, covariates=covariates, results=results)
}

#' Basic methods for class rqt
#' 
#' This document lists a series of basic methods for the class rqt
#'
setMethod("print", "rqt", function(x, pval=.05) {
  # TODO
})

#' @rdname rqt-methods
#' @aliases show.rqt
#' @return show returns summary information about the object of class rqt
#' @export
setMethod("show", "rqt", function(object) {
  print(object)
})

#' @rdname rqt-methods
#' @aliases summary.rqt
#' @return summary returns the integrated results from the analysis
#' @export
setMethod("summary", "rqt", function(object) {
  # TODO
})


