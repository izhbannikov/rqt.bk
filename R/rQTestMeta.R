#' This function performs a gene-level meta-analysis based on 
#' combined effect sizes.
#' 
#' @param x List of objects of class rqt
#' @param ... Additional parameters to pass to the function
#' @return A list of two: (i) final.pvalue - 
#' a final p-value across all studies; 
#' (ii) pvalueList - p-values for each study; 
#' @export
#' @docType methods
#' @rdname rQTestMeta-methods
setGeneric("rQTestMeta", function(x, ...) standardGeneric("rQTestMeta"))

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
#' @rdname rQTestMeta-methods
#' @examples
#'data1 <- data.matrix(read.table(system.file("extdata/phengen2.dat",
#'                                            package="rqt"), skip=1))
#'pheno <- data1[,1]
#'geno <- data1[, 2:dim(data1)[2]]
#'colnames(geno) <- paste(seq(1, dim(geno)[2]))
#'geno.obj <- SummarizedExperiment(geno)
#'obj1 <- rqtClass(phenotype=pheno, genotype=geno.obj)
#'
#'data2 <- data.matrix(read.table(system.file("extdata/phengen3.dat",
#'                                            package="rqt"), skip=1))
#'pheno <- data2[,1]
#'geno <- data2[, 2:dim(data2)[2]]
#'colnames(geno) <- paste(seq(1, dim(geno)[2]))
#'geno.obj <- SummarizedExperiment(geno)
#'obj2 <- rqtClass(phenotype=pheno, genotype=geno.obj)
#'
#'data3 <- data.matrix(read.table(system.file("extdata/phengen.dat",
#'                                            package="rqt"), skip=1))
#'pheno <- data3[,1]
#'geno <- data3[, 2:dim(data3)[2]]
#'colnames(geno) <- paste(seq(1, dim(geno)[2]))
#'geno.obj <- SummarizedExperiment(geno)
#'obj3 <- rqtClass(phenotype=pheno, genotype=geno.obj)
#'
#'res.meta <- rQTestMeta(list(obj1, obj2, obj3))
#'print(res.meta)
setMethod("rQTestMeta", signature="list", 
    function(x, perm=0, STT=0.2, weight=FALSE, 
        cumvar.threshold=75, out.type="D", 
        method="pca", scaleData=FALSE, asym.pval=FALSE,
        comb.test="wilkinson",
        verbose=FALSE) {
            
        if(cumvar.threshold > 100) {
            warning("Warning: cumvar.threshold > 100 
                       and will be set to 100.")
            cumvar.threshold <- 100
        }
        if(class(x) != "list") {
            stop("x must be a list with rqt class objects!")
        }
        ### Meta-analysis ###
        numStudies <- length(x)
        pv <- rep(NA_real_, numStudies) 
        for(i in 1:numStudies) {
            res <- rQTest(x[[i]],
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
        switch(comb.test, 
               wilkinson={
                   # Wilkinson
                   final.pvalue <- wilkinsonp(pv)$p
               },
               fisher={
                   # Fisher
                   chi.comb <- sum(-2*log(pv[!is.na(pv)]))
                   df <- 2*length(pv)
                   final.pvalue <- 1-pchisq(q=chi.comb, df=df)
               },
               minimump={
                   # minimump
                   final.pvalue <- minimump(pv)$p
               },
               sump={
                   # sump
                   final.pvalue <- sump(pv)$p
               },
               sumlog={
                   # sumlog
                   final.pvalue <- sumlog(pv)$p
               },
               meanp={
                   final.pvalue <- meanp(pv)$p
               },
               logitp={
                   final.pvalue <- logitp(pv)$p
               },
               votep={
                   final.pvalue <- votep(pv)$p
               },
               {
                   # Wilkinson
                   final.pvalue <- wilkinsonp(pv)$p
               }
        )
        
        
        #### End of combining p-values ####
        ### End of meta-analysis ###
           
        return(list(final.pvalue=final.pvalue, 
            pvalueList=pv))
})
