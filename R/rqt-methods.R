setClass("rqt", slots=c(phenotype="vector", 
                        genotype="SummarizedExperiment", 
                        covariates="data.frame", 
                        results="list"))

#' This function performs a gene-level test based on combined effect sizes.
#' 
#' @param obj Object of class \code{rqt}
#' @param ... Additional parameters to pass to the function
#' @return Updated rqt object with result slot
#' @export
#' @docType methods
#' @rdname rqt-geneTest
setGeneric("geneTest", function(obj, ...) standardGeneric("geneTest"))

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
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' res <- geneTest(obj, method="pca", out.type = "D")
#' print(res)
#' @rdname rqt-geneTest
#' @export
setMethod("geneTest", signature = "rqt", 
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
            rslt0 <- geneTestOne(phenotype=phenotype, 
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
                    t.res <- geneTestOne(phenotype=yP,genotype=genotype, 
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
                t.res <- geneTestOne(phenotype=yP,genotype=genotype, 
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
              
            rslt0 <- as.numeric(geneTestOne(phenotype=phenotype, 
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

#' This function performs a gene-level meta-analysis based on 
#' combined effect sizes.
#' 
#' @param objects List of objects of class rqt
#' @param ... Additional parameters to pass to the function
#' @return A list of two: (i) final.pvalue - 
#' a final p-value across all studies; 
#' (ii) pvalueList - p-values for each study; 
#' @export
#' @docType methods
#' @rdname rqt-geneTestMeta
setGeneric("geneTestMeta", function(objects, ...) standardGeneric("geneTestMeta"))

#' This function performs a gene-level test based on combined effect sizes.
#' 
#' @param perm Integer indicating the number of permutations 
#' to compute p-values. Default: 0.
#' @param STT Numeric indicating soft truncation threshold (STT) 
#' to convert to gamma parameter (must be <= 0.4). 
#' Needed for an optimal parameter a in Gamma-distribution. Default: 0.2. 
#' See, for example, Fridley, et al 2013: "Soft truncation thresholding 
#' for gene set analysis of RNA-seq data: Application to a vaccine study".
#' @param weight Logical value. 
#' Indicates using weights (see Lee et al 2016). 
#' Default: FALSE.
#' @param cumvar.threshold Numeric value indicating the explained 
#' variance threshold for PCA-like methods. Default: 90.
#' @param method Method used to reduce multicollinerity and account 
#' for LD. Default: PCA.
#' @param out.type Character, indicating a type of phenotype. 
#' Possible values: D (dichotomous or binary), 
#' C (continous or qualitative).
#' @param scaleData A logic parameter (TRUE/FALSE) indicating 
#' scaling of the genotype dataset.
#' @param asym.pval Indicates Monte Carlo approximation for p-values. 
#' Default: FALSE.
#' @param comb.test Statistical test for combining p-values.
#' @param verbose Indicates verbosing output. Default: FALSE.
#' @rdname rqt-geneTestMeta
#' @examples
#'data1 <- data.matrix(read.table(system.file("extdata/phengen2.dat",
#'                                            package="rqt"), skip=1))
#'pheno <- data1[,1]
#'geno <- data1[, 2:dim(data1)[2]]
#'colnames(geno) <- paste(seq(1, dim(geno)[2]))
#'geno.obj <- SummarizedExperiment(geno)
#'obj1 <- rqt(phenotype=pheno, genotype=geno.obj)
#'
#'data2 <- data.matrix(read.table(system.file("extdata/phengen3.dat",
#'                                            package="rqt"), skip=1))
#'pheno <- data2[,1]
#'geno <- data2[, 2:dim(data2)[2]]
#'colnames(geno) <- paste(seq(1, dim(geno)[2]))
#'geno.obj <- SummarizedExperiment(geno)
#'obj2 <- rqt(phenotype=pheno, genotype=geno.obj)
#'
#'data3 <- data.matrix(read.table(system.file("extdata/phengen.dat",
#'                                            package="rqt"), skip=1))
#'pheno <- data3[,1]
#'geno <- data3[, 2:dim(data3)[2]]
#'colnames(geno) <- paste(seq(1, dim(geno)[2]))
#'geno.obj <- SummarizedExperiment(geno)
#'obj3 <- rqt(phenotype=pheno, genotype=geno.obj)
#'
#'res.meta <- geneTestMeta(list(obj1, obj2, obj3))
#'print(res.meta)
setMethod("geneTestMeta", signature="list", 
          function(objects, perm=0, STT=0.2, weight=FALSE, 
                   cumvar.threshold=75, out.type="D", 
                   method="pca", scaleData=FALSE, asym.pval=FALSE,
                   comb.test="wilkinson",
                   verbose=FALSE) {
            
            if(cumvar.threshold > 100) {
              warning("Warning: cumvar.threshold > 100 
                      and will be set to 100.")
              cumvar.threshold <- 100
            }
            if(class(objects) != "list") {
              stop("objects must be a list of rqt class objects!")
            }
            ### Meta-analysis ###
            numStudies <- length(objects)
            pv <- rep(NA_real_, numStudies) 
            for(i in 1:numStudies) {
              res <- geneTest(objects[[i]],
                              STT=STT, 
                              weight=weight, 
                              cumvar.threshold=cumvar.threshold, 
                              out.type=out.type, 
                              method=method, 
                              perm = perm, 
                              scaleData=scaleData,
                              asym.pval=asym.pval,
                              verbose=verbose)
              
              if(length(results(res)) != 0) {
                pv[i] <- results(res)$p.value$p.Q3
              }
              
            }
            
            ### Combining p-values via some comb.test ###
            comb.res <- list()
            switch(comb.test, 
                   wilkinson={
                     # Wilkinson
                     comb.res <- wilkinsonp(pv)
                   },
                   fisher={
                     # Fisher
                     chi.comb <- sum(-2*log(pv[!is.na(pv)]))
                     df <- 2*length(pv)
                     comb.res[["p"]] <- 1-pchisq(q=chi.comb, df=df)
                   },
                   minimump={
                     # minimump
                     comb.res <- minimump(pv)
                   },
                   sump={
                     # sump
                     comb.res <- sump(pv)
                   },
                   sumlog={
                     # sumlog
                     comb.res <- sumlog(pv)
                   },
                   meanp={
                     comb.res <- meanp(pv)
                   },
                   logitp={
                     comb.res <- logitp(pv)
                   },
                   votep={
                     comb.res <- votep(pv)
                   },
                   {
                     # Wilkinson
                     comb.res <- wilkinsonp(pv)
                   }
            )
            
            #### End of combining p-values ####
            ### End of meta-analysis ###
            
            return(list(final.pvalue=comb.res[["p"]], 
                        pvalueList=pv))
          })


#' Common methods for class rqt
#' 
#' @name rqt-general
#' @rdname rqt-general
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
#' @rdname rqt-phenotype
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
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' phenotype(obj)
#' @rdname rqt-phenotype
#' @export
setMethod("phenotype", signature(obj = "rqt"), function(obj) {
  return(slot(obj, "phenotype"))
})


#' This function performs an access to genotype.
#' 
#' @rdname rqt-genotype
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
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' genotype(obj)
#' @rdname rqt-genotype
#' @export
setMethod("genotype", signature(obj = "rqt"), function(obj) {
  return(slot(obj, "genotype"))
})


#' This function performs an access to covariates
#' 
#' @rdname rqt-covariates
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
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' covariates(obj)
#' @rdname rqt-covariates
#' @export
setMethod("covariates", signature(obj = "rqt"), function(obj) {
  return(slot(obj, "covariates"))
})

#' This function performs an access to covariates
#' 
#' @rdname rqt-results
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
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' res <- geneTest(obj, method="pca", out.type = "D")
#' results(res)
#' @rdname rqt-results
#' @export
setMethod("results", signature(obj = "rqt"), function(obj) {
  return(slot(obj, "results"))
})

setGeneric("results<-", function(obj, value) standardGeneric("results<-"))
setReplaceMethod("results", signature(obj = "rqt"), 
                 function(obj, value) {
                     slot(obj, "results") <- value
                     obj
                 }
)

setMethod("show", signature(object = "rqt"), function(object) {
  cat("Phenotype:\n")
  print(head(phenotype(object)))
  cat("...\n\n")
  cat("Genotype:\n")
  print(head(assays(genotype(object))[[1]]))
  cat("...\n\n")
  cat("Covariates:\n")
  print(head(covariates(object)))
  cat("\n\n")
  cat("Results:\n\n")
  print(results(object))
})

setMethod("summary", signature(object = "rqt"), function(object) {
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
