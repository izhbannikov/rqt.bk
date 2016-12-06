#' This function performs a gene-level meta-analysis based on 
#' combined effect sizes.
#' 
#' @param x List of objects of class rqt
#' @param ... Additional parameters to pass to the function
#' @return A list of four: (i) final.pvalue - 
#' a final p-value across all studies; 
#' (ii) pvalueList - p-values for each study; 
#' (iii) df - degrees of freedom; 
#' (iv) chi.comb - meta-analysis Chi-square statistics
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
#' @param weight Logical value. Indicates using weights (see Lee et al 2016). 
#' Default: FALSE.
#' @param cumvar.threshold Numeric value indicating the explained 
#' variance threshold for PCA-like methods. Default: 90.
#' @param method Method used to reduce multicollinerity and account 
#' for LD. Default: PCA.
#' @param out.type Character, indicating a type of phenotype. 
#' Possible values: D (dichotomous or binary), 
#' C (continous or qualitative).
#' @param scale A logic parameter (TRUE/FALSE) indicating 
#' scaling of the genotype dataset.
#' @rdname rQTestMeta-methods
#' @examples
#'data1 <- read.table(system.file("data/phengen2.dat",package="rqt"), skip=1)
#'obj1 <- rqtClass(phenotype=data1[,1], genotype=data1[, 2:dim(data1)[2]])
#'data2 <- read.table(system.file("data/phengen3.dat",package="rqt"), skip=1)
#'obj2 <- rqtClass(phenotype=data2[,1], genotype=data2[, 2:dim(data2)[2]])
#'data3 <- read.table(system.file("data/phengen.dat",package="rqt"), skip=1)
#'obj3 <- rqtClass(phenotype=data3[,1], genotype=data3[, 2:dim(data3)[2]])
#'# Gene-level meta-analysis:
#'res <- rQTestMeta(list(obj1, obj2, obj3))
#'print(res)
setMethod("rQTestMeta", signature="list", 
          function(x, perm=0, STT=0.2, weight=FALSE, 
                   cumvar.threshold=90, out.type="D", 
                   method="pca", scale=FALSE) {
            
            if(cumvar.threshold > 100) {
              cat("Warning: cumvar.threshold > 100 
                  and will be set to 100.")
              cumvar.threshold <- 100
            }
            if(class(x) != "list") {
              stop("x must be a list with rqt class objects!")
            }
            ### Meta-analysis ###
            pv <- c()
            for(i in 1:length(x)) {
              res <- rQTest(x[[i]],
                            STT=STT, 
                            weight=weight, 
                            cumvar.threshold=cumvar.threshold, 
                            out.type=out.type, 
                            method=method, 
                            perm = perm, 
                            scale=scale)
              
              if(length(res@results) != 0) {
                pv <- c(pv, res@results$p.value$p.Q3)
              }
           }
            
           ### Combining p-values via Fisher's test ###
           chi.comb <- sum(-2*log(pv[which(!is.na(pv))]))
           df <- 2*length(pv)
           final.pvalue <- 1-pchisq(q = chi.comb, df=df)
           #### End of combining p-values ####
           ### End of meta-analysis ###
           
           return(list(final.pvalue=final.pvalue, 
                       pvalueList=pv, df=df, 
                       chi.comb=chi.comb))
})
