
rQTestMeta <- function(data, 
                       perm=0, 
                       STT=0.2, 
                       weight=FALSE, 
                       cumvar.threshold=90, 
                       out.type="D", 
                       method="pca") {
  ### Meta-analysis ###
  pv <- c()
  for(d in data) {
    dimensions <- dim(d)
    pheno <- d[,1]
    geno <- d[,2:dimensions[2]]
    res <- rQTest(pheno = pheno, geno=geno, STT=STT, weight=weight, cumvar.threshold=cumvar.threshold, out.type=out.type, method=method)
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
