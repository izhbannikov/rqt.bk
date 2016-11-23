#' This function performs a gene-level test based on combined effect sizes.
#' 
#' @param data A matrix n by 1 where n - number of individuals).
#' @param covariates A matrix of covariates. Default: NULL.
#' @param perm Integer indicating the number of permutations to compute p-values. Default: 0.
#' @param STT Numeric indicating soft truncation threshold (STT) to convert to gamma parameter (must be <= 0.4). 
#' Needed for an optimal parameter a in Gamma-distribution. Default: 0.2. 
#' See, for example, Fridley, et al 2013: "Soft truncation thresholding for gene set analysis of RNA-seq data: Application to a vaccine study".
#' @param weight Logical value. Indicates using weights (see Lee et al 2016). Default: FALSE.
#' @param cumvar.threshold Numeric value indicating the explained variance threshold for PCA-like methods. Default: 90.
#' @param method Method used to reduce multicollinerity and account for LD. Default: PCA.
#' @param out.type Character, indicating a type of phenotype. Possible values: D (dichotomous or binary), 
#' C (continous or qualitative).
#' @param scale A logic parameter (TRUE/FALSE) indicating scaling of the genotype dataset.
#' @return A list of two: test statistics (Q1, Q2, Q3), p-values (p1.Q1, p2.Q2, p3.Q3).
#' @description 
#' @export
#' @examples
#' 
#'library(rqt)
#'data1 <- read.table(system.file("data/phengen2.dat",package="rqt"), skip=1)
#'data2 <- read.table(system.file("data/phengen3.dat",package="rqt"), skip=1)
#'data3 <- read.table(system.file("data/phengen.dat",package="rqt"), skip=1)
#'# Combining datasets:
#'indata <- list(data1, data2, data3)
#'# Gene-level meta-analysis:
#'res <- rQTestMeta(indata)
#'res
rQTestMeta <- function(data, 
                       covariates=NULL,
                       perm=0, 
                       STT=0.2, 
                       weight=FALSE, 
                       cumvar.threshold=90, 
                       out.type="D", 
                       method="pca",
                       scale=FALSE) {
  ### Meta-analysis ###
  pv <- c()
  for(i in 1:length(data)) {
    dimensions <- dim(data[[i]])
    pheno <- data[[i]][,1]
    geno <- data[[i]][,2:dimensions[2]]
    if(!is.null(covariates)) {
      res <- rQTest(phenotype = pheno, genotype = geno, 
                  covariates = covariates[[i]],
                  STT=STT, weight=weight, cumvar.threshold=cumvar.threshold, out.type=out.type, method=method, 
                  perm = perm, scale=scale)
    } else {
      res <- rQTest(phenotype = pheno, genotype = geno, 
                    covariates = covariates,
                    STT=STT, weight=weight, cumvar.threshold=cumvar.threshold, out.type=out.type, method=method, 
                    perm = perm, scale=scale)
    }
    
    if(!is.na(res)) {
      pv <- c(pv, res$p.value$p.Q3)
    }
  }
  
  ### Combining p-values via Fisher's test ###
  chi.comb <- sum(-2*log(pv[which(!is.na(pv))]))
  df <- 2*length(pv)
  final.pvalue <- 1-pchisq(q = chi.comb, df=df)
  #### End of combining p-values ####
  ### End of meta-analysis ###
  
  return(list(final.pvalue=final.pvalue, pvalueList=pv, df=df, chi.comb=chi.comb))
}
