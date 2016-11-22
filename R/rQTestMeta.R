
rQTestMeta <- function(data, 
                       covariates,
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
    res <- rQTest(pheno = pheno, geno=geno, covariates[[i]], STT=STT, weight=weight, cumvar.threshold=cumvar.threshold, out.type=out.type, method=method, perm = perm, scale=scale)
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
