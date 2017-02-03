#' This function performs a gene-level test based on combined effect sizes.
#' 
#' @param obj Object of class rqt
#' @param ... Additional parameters to pass to the function
#' @return Updated rqt object with result slot
#' @export
#' @docType methods
#' @rdname rQTest-methods
setGeneric("rQTest", function(obj, ...) standardGeneric("rQTest"))

#' This function performs a gene-level test based on combined effect sizes.
#' @param perm Integer indicating the number of permutations 
#' to compute p-values. Default: 0.
#' @param STT Numeric indicating soft truncation threshold (STT) 
#' to convert to gamma parameter (must be <= 0.4). 
#' Needed for an optimal parameter a in Gamma-distribution. Default: 0.2. 
#' See, for example, Fridley, et al 2013: "Soft truncation thresholding 
#' for gene set analysis of RNA-seq data: Application to a vaccine study".
#' @param weight Logical value. Indicates using weights (see Lee et al 2016). 
#' Default: FALSE.
#' @param cumvar.threshold Numeric value indicating 
#' the explained variance threshold for PCA-like methods. Default: 75
#' @param method Method used to reduce multicollinerity and account for LD. 
#' Default: PLS-DA.
#' @param out.type Character, indicating a type of phenotype. 
#' Possible values: \code{D} (dichotomous or binary), 
#' \code{C} (continous or quantitative).
#' @param scaleData A logic parameter (TRUE/FALSE) indicating scaling of 
#' the genotype dataset.
#' @param asym.pval Indicates Monte Carlo approximation for p-values. Default: FALSE.
#' @param verbose Indicates verbosing output. Default: FALSE.
#' @examples
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat",
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqtClass(phenotype=pheno, genotype=geno.obj)
#' res <- rQTest(obj, method="pca", out.type = "D")
#' print(res)
#' @rdname rQTest-methods
#' @export
setMethod("rQTest", signature="rqt", 
    function(obj, perm=0, STT=0.2, weight=FALSE, 
            cumvar.threshold=75, out.type="D", 
            method="pca", scaleData=FALSE, asym.pval=FALSE,
            verbose=FALSE) {
            # Prepare test: load distribution table and prepare #
            # some other information #
        if(cumvar.threshold > 100) {
            warning("Warning: cumvar.threshold > 100 and will be set to 100.")
            cumvar.threshold <- 100
        }
      
        # Load data #
        phenotype <- phenotype(obj)
        genotype <- assays(genotype(obj))[[1]]
        covariates <- covariates(obj)
        
        # Dimensions
        phenoSize <- length(phenotype)
        genoSize <- dim(genotype)
        
        if(weight & scaleData) {
          if(verbose) {
              print("Warining! You can not use scaling in presence of weights!")
              print("Parameter weight will be set to FALSE.")
          }
          weight <- FALSE
        }
        
        # Start the tests #
        if(perm==0){
            rslt0 <- QTest.one(phenotype=phenotype, 
                    genotype=genotype, 
                    covariates=covariates, 
                    STT=STT, 
                    weight=weight,
                    cumvar.threshold=cumvar.threshold, 
                    out.type=out.type, method=method, 
                    scaleData = scaleData, verbose=verbose)
            
            if(asym.pval) {
                rsltMC <- do.call(rbind, 
                                  lapply(1:dim(genotype)[1], 
                                         function(k){
                    yP <- phenotype[sample(1:phenoSize, phenoSize, 
                                           replace=FALSE)]
                    t.res <- QTest.one(phenotype=yP,genotype=genotype, 
                                 covariates=covariates,STT=STT,
                                 weight=weight,
                                 cumvar.threshold=cumvar.threshold, 
                                 out.type=out.type, method=method, 
                                 scaleData = scaleData)
                    if(is.na(t.res)) {
                        t.res <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                               data.frame(p.Q1=1,p.Q2=1,p.Q3=1) )
                        names(t.res) <- c("Qstatistic", "p.value")
                    }
                    tt.res <- t.res$Qstatistic
                }))
            
                if(!is.na(rslt0)) {
                    nn <- genoSize[1]+1
                    rslt <- list("Qstatistic"= data.frame(
                            Q1=rslt0$Qstatistic$Q1,
                            Q2=rslt0$Qstatistic$Q2, 
                            Q3=rslt0$Qstatistic$Q3),
                           "p.value" = data.frame(
                             p.Q1 = (length(rsltMC[,1][rsltMC[,1] >= 
                                            rslt0$Qstatistic$Q1])+1)/nn,
                             p.Q2 = (length(rsltMC[,2][rsltMC[,2] >= 
                                            rslt0$Qstatistic$Q2])+1)/nn,
                             p.Q3 = (length(rsltMC[,3][rsltMC[,3] >= 
                                            rslt0$Qstatistic$Q3])+1)/nn),
                           beta = rslt0$beta)
                } else {
                    rslt <- list( Qstatistic=data.frame(Q1=NA, Q2=NA, Q3=NA),
                                p.value=data.frame(p.Q1=1,p.Q2=1,p.Q3=1) )
                }
            } else {
                rslt <- rslt0
            }
        } else {
            rsltPP <- do.call(rbind, lapply(1:perm, function(k){
                yP <- phenotype[sample(1:phenoSize, phenoSize,replace=FALSE)]
                t.res <- QTest.one(phenotype=yP,genotype=genotype, 
                    covariates=covariates,STT=STT,
                    weight=weight,
                    cumvar.threshold=cumvar.threshold, 
                    out.type=out.type, method=method, 
                    scaleData = scaleData, verbose=verbose)
                if(is.na(t.res)) {
                    t.res <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                        data.frame(p.Q1=1,p.Q2=1,p.Q3=1) )
                    names(t.res) <- c("Qstatistic", "p.value")
                }
                pv <- t.res$p.value
            }))
              
            rslt0 <- as.numeric(QTest.one(phenotype=phenotype, 
                genotype=genotype, 
                covariates=covariates, 
                STT=STT, weight=weight, 
                cumvar.threshold=cumvar.threshold, 
                out.type=out.type, method=method, 
                scaleData = scaleData))
              
            if(!is.na(rslt0)) {
                rslt <- list(Qstatistic= data.frame(Q1=NA, Q2=NA, Q3=NA),
                    p.value = data.frame(
                    p.Q1 = (length(rsltPP[,1][rsltPP[,1] < 
                                              rslt0$p.value[1]])+1)/(perm+1),
                    p.Q2 = (length(rsltPP[,2][rsltPP[,2] < 
                                              rslt0$p.value[2]])+1)/(perm+1),
                    p.Q3 = (length(rsltPP[,3][rsltPP[,3] < 
                                              rslt0$p.value[3]])+1)/(perm+1)),
                    beta = rslt0$beta)
            } else {
                rslt <- list( Qstatistic=data.frame(Q1=NA, Q2=NA, Q3=NA),
                            p.value=data.frame(p.Q1=1,p.Q2=1,p.Q3=1) )
            }
        }
        
        results(obj) <- rslt
        return(obj)
})

#' Common methods for class rqt
#' 
#' @name rQTest-general
#' @rdname rQTest-general
#'
#' @aliases show.rqt
#' @aliases summary.rqt
#' @aliases print.rqt
#' 
#' @title General functions of \code{rqt} 
#' such as accessors and printing.
#' 
#' @description Common methods for class rqt. 
#' This document lists a series of basic methods for the class rqt
#' 
NULL


#' This function performs an access to phenotype
#' 
#' @rdname rQTest-phenotype
#' @export
setGeneric("phenotype", function(obj) standardGeneric("phenotype"))

#' This function performs an access to phenotype
#' 
#' @description A genotype accessor
#' @param obj An object of \code{rqt} class.
#' @return phenotype returns the genotype
#' @examples 
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat", 
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqtClass(phenotype=pheno, genotype=geno.obj)
#' phenotype(obj)
#' @rdname rQTest-phenotype
#' @export
setMethod("phenotype", "rqt", function(obj) {
  return(slot(obj, "phenotype"))
})


#' This function performs an access to genotype.
#' 
#' @rdname rQTest-genotype
#' @export
setGeneric("genotype", function(obj) standardGeneric("genotype"))

#' This function performs an access to genotype.
#' 
#' @description A genotype accessor
#' @param obj An object of \code{rqt} class.
#' @return genotype returns the genotype
#' @examples 
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat", 
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqtClass(phenotype=pheno, genotype=geno.obj)
#' genotype(obj)
#' @rdname rQTest-genotype
#' @export
setMethod("genotype", "rqt", function(obj) {
  return(slot(obj, "genotype"))
})


#' This function performs an access to covariates
#' 
#' @rdname rQTest-covariates
#' @export
setGeneric("covariates", function(obj) standardGeneric("covariates"))

#' This function performs an access to covariates
#' 
#' @description An accessor to covariates
#' @param obj An object of \code{rqt} class.
#' @return covariates returns the covariates
#' @examples 
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat", 
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqtClass(phenotype=pheno, genotype=geno.obj)
#' covariates(obj)
#' @rdname rQTest-covariates
#' @export
setMethod("covariates", "rqt", function(obj) {
  return(slot(obj, "covariates"))
})

#' This function performs an access to covariates
#' 
#' @rdname rQTest-results
#' @export
setGeneric("results", function(obj) standardGeneric("results"))


#' This function performs an access to covariates
#' 
#' @description An accessor to results
#' @param obj An object of \code{rqt} class.
#' @return results returns the results
#' @examples 
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat", 
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqtClass(phenotype=pheno, genotype=geno.obj)
#' res <- rQTest(obj, method="pca", out.type = "D")
#' results(res)
#' @rdname rQTest-results
#' @export
setMethod("results", "rqt", function(obj) {
  return(slot(obj, "results"))
})

setGeneric("results<-", function(obj, value) standardGeneric("results<-"))

setMethod("results<-", "rqt", function(obj, value) {
  slot(obj, "results") <- value
  obj
})

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

setMethod("show", "rqt", function(object) {
  print(object)
})

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
