# Clases

#'The rqt class
#'
#'This class stores parameters and results of the rtq algorithms
#'
#'@section Slots:
#'    \describe{
#'    \item{\code{phenotype}:}{Phenotype (a vector of length 
#'        \code{N}, where \code{N} - number of individuals).}
#'      \item{\code{genotype}:}{Genotype - an object of class 
#'      \code{SummarizedExperiment}. Should contain one assay 
#'      (matrix, \code{N} 
#'      by \code{M} where \code{N} - number of individuals, \code{M}
#'       - number of genetic variants).}
#'      \item{\code{covariates}:}{data frame \code{N} 
#'      by \code{K} where \code{N} - number of individuals, \code{K}
#'       - number of covariates)}
#'      \item{\code{results}:}{A list of two: 
#'      test statistics (\code{Q1}, \code{Q2}, \code{Q3}), 
#'      p-values (\code{p1.Q1}, \code{p2.Q2}, \code{p3.Q3})}
#'}
#'@rdname rqt-class
#'@details see \code{\link[=print.rqt]{rqt-methods} for related methods}
setClass("rqt", slots=c(phenotype="vector", genotype="SummarizedExperiment", covariates="data.frame", results="list"))

#' The rqt class constructor
#' 
#' This function generates rqt class objects
#' @param  phenotype Phenotype (a vector of length 
#'        \code{N}, where \code{N} - number of individuals).
#' @param genotype Genotype - an object of class 
#'      \code{SummarizedExperiment}. Should contain one assay 
#'      (matrix, \code{N} 
#'      by \code{M} where \code{N} - number of individuals, \code{M}
#'       - number of genetic variants).
#' @param covariates Covariates, a data frame \code{N} 
#'      by \code{K} where \code{N} - number of individuals, \code{K}
#'       - number of covariates
#' @param results A list of two: test statistics: 
#' (\code{Q1}, \code{Q2}, \code{Q3}), 
#' p-values: (\code{p1.Q1}, \code{p2.Q2}, \code{p3.Q3})
#' @return Object of class rqt
#' @examples
#'data <- data.matrix(read.table(system.file("extdata/test.bin1.dat",
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqtClass(phenotype=pheno, genotype=geno.obj)
#' print(obj)
#' @rdname rqt-class
#' @aliases rqtClass
#' @export
rqtClass <- function(phenotype=NULL, genotype=NULL, covariates=NULL, 
    results=NULL) {
    
    if(is.null(phenotype)) {
        phenotype <- c()
    }
  
    if(is.null(genotype)) {
        genotype <- matrix()
        colnames(genotype) <- "geno"
        genotype <- SummarizedExperiment(genotype)
    } 
  
    if(is.null(covariates)) {
        covariates <- data.frame()
    }
    
    if(is.null(results)) {
        results <- list()
    }

    new("rqt", phenotype=phenotype, genotype=genotype, covariates=covariates, results=results)
}

setGeneric("phenotype", function(x) standardGeneric("phenotype"))
setGeneric("genotype", function(x) standardGeneric("genotype"))
setGeneric("covariates", function(x) standardGeneric("covariates"))
setGeneric("results", function(x) standardGeneric("results"))
setGeneric("results<-", function(x, value) standardGeneric("results<-"))

#' @rdname rqt-methods
#' @aliases phenotype.rqt
#' @return pheotype returns the phenotype
#' @docType methods
#' @export
setMethod("phenotype", "rqt", function(x) {
    return(slot(x, "phenotype"))
})

#' @rdname rqt-methods
#' @aliases genotype.rqt
#' @return genotype returns the genotype
#' @docType methods
#' @export
setMethod("genotype", "rqt", function(x) {
    return(slot(x, "genotype"))
})

#' @rdname rqt-methods
#' @aliases covariates.rqt
#' @return covariates returns the covariates
#' @docType methods
#' @export
setMethod("covariates", "rqt", function(x) {
    return(slot(x, "covariates"))
})

#' @rdname rqt-methods
#' @aliases results.rqt
#' @return results returns the results
#' @docType methods
#' @export
setMethod("results", "rqt", function(x) {
    return(slot(x, "results"))
})

#' @rdname rqt-methods
#' @aliases results<-.rqt
#' @docType methods
#' @export
setMethod("results<-", "rqt", function(x, value) {
    slot(x, "results") <- value
    x
})


#' Basic methods for class rqt
#' 
#' This document lists a series of basic methods for the class rqt
#' 
#'
#' @rdname rqt-methods
#' @param x An object of \code{rqt} class
#' @return print returns summary information about the rqt object
#' @rdname rqt-methods
#' @aliases print.rqt
#' @docType methods
setMethod("print", "rqt", function(x) {
    cat("Phenotype:\n")
    print(head(phenotype(x)))
    cat("...\n\n")
    cat("Genotype:\n")
    print(head(assays(genotype(x))[[1]]))
    cat("...\n\n")
    cat("Covariates:\n")
    print(head(covariates(x)))
    cat("\n\n")
    cat("Results:\n\n")
    print(results(x))
})

#' @rdname rqt-methods
#' @aliases show.rqt
#' @return show returns summary information about 
#' the object of class rqt
#' @export
setMethod("show", "rqt", function(object) {
    print(object)
})

#' @rdname rqt-methods
#' @aliases summary.rqt
#' @return summary returns the integrated results from the analysis
#' @docType methods
#' @export
setMethod("summary", "rqt", function(object) {
    cat("Phenotype:\n")
    print(summary(phenotype(object)))
    cat("...\n\n")
    cat("Genotype:\n")
    print(summary(assays(genotype(object))[[1]]))
    cat("...\n\n")
    cat("Covariates:\n")
    cat(summary(covariates(object)))
    cat("\n\n")
    cat("Results:\n\n")
    print(summary(results(object)))
})