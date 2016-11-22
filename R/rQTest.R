
#' Empirical null distribution for Q3 test.
#' 
rQTtest.prepare <- function() {
  null.dist.Q3 <- read.table(system.file("data/n.log10.minp.1e09.txt",package="rqt"), header=T)
  assign("null.dist.Q3", null.dist.Q3, envir=baseenv())
}

#rQTtest.prepare <- function() {
#  #null.dist.Q3 <- read.table("W:/data/work/iz12/rqt/package/rqt-master/data/n.log10.minp.1e09.txt", header=T)
#  null.dist.Q3 <- read.table("W:/data/work/iz12/rqt/package/rqt-master/data/n.log10.minp.1e09.txt", header=T)
#  assign("null.dist.Q3", null.dist.Q3, envir=baseenv())
#}


#' This function performs a gene-level test based on combined effect sizes.
#' 
#' @param phenotype Phenotype (a matrix n by 1 where n - number of individuals).
#' @param genotype Genotype (matrix n by m where n - number of individuals, m - number of genetic variants).
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
#' library(rqt)
#' data <- read.table(system.file("data/phengen2.dat",package="rqt"), skip=1)
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' res <- rQTest(phenotype=pheno, genotype=geno)
rQTest <- function(phenotype, genotype, covariates=NULL, perm=0, STT=0.2,weight=FALSE, cumvar.threshold=90, out.type="D", method="pca", scale=FALSE) {
  # Prepare test: load distribution table and prepare some other information#
  rQTtest.prepare()
  
  num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
  
  if(perm==0){
    rslt <- try(QTest.one(phenotype=phenotype, genotype=genotype, covariates=covariates, STT=STT, cumvar.threshold=cumvar.threshold, out.type=out.type, method=method, scale = scale),TRUE)
  } else {
    rsltPP <- do.call(rbind, lapply(1:perm, function(k){
                                               yP <- phenotype[sample(1:length(phenotype),length(phenotype),replace=F)]
                                               t.res <- QTest.one(phenotype=yP,genotype=genotype, covariates=covariates,STT=STT,weight=weight, cumvar.threshold=cumvar.threshold, out.type=out.type, method=method, scale = scale)
                                               if(is.na(t.res)) {
                                                 t.res <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), data.frame(p.Q1=1,p.Q2=1,p.Q3=1) )
                                                 names(rslt)<-c("Qstatistic", "p.value")	
                                               } 
                                                          
                                               pv <- t.res$p.value
                                                          
                                            }
                                    )
                      )
    
    rslt0 <- try(as.numeric(QTest.one(phenotype=phenotype, genotype=genotype, covariates=covariates, STT=STT, weight=weight, cumvar.threshold=cumvar.threshold, out.type=out.type, method=method, scale = scale)$p.value),TRUE)
    
    if(!is.na(rslt0)) {
      
      rslt <- list(
          "Qstatistic"= data.frame(Q1=NA, Q2=NA, Q3=NA),
          "p.value" = data.frame(p.Q1 = (length(which(rsltPP[,1] < rslt0[1])) + 1) / (perm + 1),
                         p.Q2 = (length(which(rsltPP[,2] < rslt0[2])) + 1) / (perm + 1),
                         p.Q3 = (length(which(rsltPP[,3] < rslt0[3])) + 1) / (perm + 1))
      )
    } else {
      rslt <- NA
    }
  }
  
  rslt
}
