#' This function performs a gene-level test based on combined effect sizes.
#' 
#' @param x Object of class rqt
#' @param ... Additional parameters to pass to the function
#' @return Updated rqt object with result slot
#' @export
#' @docType methods
#' @rdname rQTest-methods
setGeneric("rQTest", function(x, ...) standardGeneric("rQTest"))

#' This function performs a gene-level test based on combined effect sizes.
#' 
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
#' @param scale A logic parameter (TRUE/FALSE) indicating scaling of 
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
    function(x, perm=0, STT=0.2, weight=FALSE, 
            cumvar.threshold=75, out.type="D", 
            method="pca", scale=FALSE, asym.pval=FALSE,
            verbose=FALSE) {
            # Prepare test: load distribution table and prepare #
            # some other information #
        if(cumvar.threshold > 100) {
            warning("Warning: cumvar.threshold > 100 and will be set to 100.")
            cumvar.threshold <- 100
        }
      
        # Load data #
        phenotype <- phenotype(x)
        genotype <- assays(genotype(x))[[1]]
        covariates <- covariates(x)
        
        if((weight == TRUE) & (scale == TRUE)) {
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
                    scale = scale, verbose=verbose)
            
            if(asym.pval) {
                rsltMC <- do.call(rbind, lapply(1:dim(genotype)[1], function(k){
                    yP <- phenotype[sample(1:length(phenotype),
                                     length(phenotype),
                                     replace=FALSE)]
                    t.res <- QTest.one(phenotype=yP,genotype=genotype, 
                                 covariates=covariates,STT=STT,
                                 weight=weight,
                                 cumvar.threshold=cumvar.threshold, 
                                 out.type=out.type, method=method, 
                                 scale = scale)
                    if(is.na(t.res)) {
                        t.res <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                               data.frame(p.Q1=1,p.Q2=1,p.Q3=1) )
                        names(t.res) <- c("Qstatistic", "p.value")
                    }
                    tt.res <- t.res$Qstatistic
                }))
            
                if(!is.na(rslt0)) {
                    rslt <- list("Qstatistic"= data.frame(Q1=rslt0$Qstatistic$Q1, 
                                                    Q2=rslt0$Qstatistic$Q2, 
                                                    Q3=rslt0$Qstatistic$Q3),
                           "p.value" = data.frame(
                             p.Q1 = (length(which(rsltMC[,1] >= 
                                                    rslt0$Qstatistic$Q1))+1)/(dim(genotype)[1]+1),
                             p.Q2 = (length(which(rsltMC[,2] >= 
                                                    rslt0$Qstatistic$Q2))+1)/(dim(genotype)[1]+1),
                             p.Q3 = (length(which(rsltMC[,3] >= 
                                                    rslt0$Qstatistic$Q3))+1)/(dim(genotype)[1]+1)),
                           beta = rslt0$beta)
                } else {
                    rslt <- NA
                }
            } else {
                rslt <- rslt0
            }
        } else {
            rsltPP <- do.call(rbind, lapply(1:perm, function(k){
                yP <- phenotype[sample(1:length(phenotype),
                    length(phenotype),
                    replace=FALSE)]
                t.res <- QTest.one(phenotype=yP,genotype=genotype, 
                    covariates=covariates,STT=STT,
                    weight=weight,
                    cumvar.threshold=cumvar.threshold, 
                    out.type=out.type, method=method, 
                    scale = scale, verbose=verbose)
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
                scale = scale))
              
            if(!is.na(rslt0)) {
                rslt <- list(Qstatistic= data.frame(Q1=NA, Q2=NA, Q3=NA),
                    p.value = data.frame(
                    p.Q1 = (length(which(rsltPP[,1] < 
                    rslt0$p.value[1]))+1)/(perm+1),
                    p.Q2 = (length(which(rsltPP[,2] < 
                    rslt0$p.value[2]))+1)/(perm+1),
                    p.Q3 = (length(which(rsltPP[,3] < 
                    rslt0$p.value[3]))+1)/(perm+1)),
                    beta = rslt0$beta)
            } else {
                rslt <- NA
            }
        }
        
        results(x) <- rslt
        return(x)
})
