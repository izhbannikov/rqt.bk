
#' Empirical null distribution for Q3 test.
#' 
rQTtest.prepare <- function() {
  null.dist.Q3 <- read.table(system.file("data/n.log10.minp.1e09.txt",package="rqt"), header=T)
  assign("null.dist.Q3", null.dist.Q3, envir=baseenv())
}


#' This function performs a gene-level test based on combined effect sizes.
#' 
#' @param pheno Phenotype (vector).
#' @param geno Genotype (matrix n by m where n - number of individuals, m - number of genetic variants)/
#' @param perm Integer indicating the number of permutations to compute p-values. Default: 0.
#' @param STT  TODO
#' @param weight Logical. TODO
#' @param cumvar.threshold Numeric value indicating the explained variance threshold for PCA-like methods. Default: 90.
#' @param out.type Character, indicating a type of phenotype. Possible values: D (dichotomous or binary), 
#' C (continous or qualitative).
#' @param method Method used to reduce multicollinerity and account for LD. Default: PCA.
#' @return A list with p-values and corresponding statistics.
#' @description TODO
#' @export
#' @examples
#' library(rqt)
#' data <- read.table(system.file("data/phengen2.dat",package="rqt"), skip=1)
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' res <- rQTest(pheno=pheno, geno=geno)
rQTest <- function(pheno, geno, perm=0, STT=0.2,weight=FALSE, cumvar.threshold=90, out.type="D", method="pca") {
  # Prepare test: load distribution table and prepare some other information#
  rQTtest.prepare()
  
  num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
  
  if(perm==0){
    rslt <- try(QTest.one(y=pheno,covadat=NULL,newgeno=geno, STT=STT, cumvar.threshold=cumvar.threshold, out.type=out.type, method=method),TRUE)
  } else {
    rsltPP <- do.call(rbind, lapply(1:perm, function(k){
                                                          yP <- pheno[sample(1:length(pheno),length(pheno),replace=F)]
                                                          t.res <- QTest.one(y=yP,covadat=NULL, newgeno=geno,STT=STT,weight=weight, cumvar.threshold=cumvar.threshold, out.type=out.type, method=method)
                                                          if(is.na(t.res)) {
                                                            t.res <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), data.frame(p.Q1=1,p.Q2=1,p.Q3=1) )
                                                            names(rslt)<-c("Qstatistic", "p.value")	
                                                          } 
                                                          
                                                          pv <- t.res$p.value
                                                          print(pv)
                                              })
                      )
    
    rslt0 <- try(as.numeric(QTest.one(y=pheno, covadat=NULL, newgeno=geno, STT=STT, weight=weight, cumvar.threshold=cumvar.threshold, out.type=out.type, method=method)$p.value),TRUE)
    if(!is.na(rslt0)) {
      
      rslt <- data.frame(p.Q1 = (length(which(rsltPP[,1] < rslt0[1])) + 1) / (perm + 1),
                         p.Q2 = (length(which(rsltPP[,2] < rslt0[2])) + 1) / (perm + 1),
                         p.Q3 = (length(which(rsltPP[,3] < rslt0[3])) + 1) / (perm + 1))
    } else {
      rslt <- NA
    }
  }
  
  rslt
}
