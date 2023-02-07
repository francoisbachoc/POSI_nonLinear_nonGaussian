rm( list=ls() )   

vName = c("logistic_lasso_pairwiseCorrelated_wellSpecified_zeroBeta_smallLambda_n30",
          "logistic_lasso_pairwiseCorrelated_wellSpecified_zeroBeta_smallLambda_n100",
          "logistic_lasso_pairwiseCorrelated_wellSpecified_sparseBeta_smallLambda_n30",
          "logistic_lasso_pairwiseCorrelated_wellSpecified_sparseBeta_smallLambda_n100",
          "logistic_lasso_pairwiseCorrelated_wellSpecified_scaledBeta_smallLambda_n30",
          "logistic_lasso_pairwiseCorrelated_wellSpecified_scaledBeta_smallLambda_n100",
          "logistic_lasso_pairwiseCorrelated_wellSpecified_scaledBeta_largeLambda_n30",
          "logistic_lasso_pairwiseCorrelated_wellSpecified_scaledBeta_largeLambda_n100",
          "logistic_lasso_pairwiseCorrelated_misspecified_denseBeta_smallLambda_n30",
          "logistic_lasso_pairwiseCorrelated_misspecified_denseBeta_smallLambda_n100",
          "logistic_lasso_pairwiseCorrelated_misspecified_denseBeta_largeLambda_n30",
          "logistic_lasso_pairwiseCorrelated_misspecified_denseBeta_largeLambda_n100")

for (name in vName) {
  load(paste0(name,".Rdata"))
  cat("############################## \n")
  cat(name,"\n")
  cat("n=",n,"p=",p,"lp=",lp,"alpha=",alpha,"rho=",rho,"lambda=",lambda,"\n")
  cat("lbeta0=",lbeta0,"\n")
  cat("############################## \n")
  failurePOSI = is.na(vLengthPOSI) | is.na(vCoveragePOSI)
  failureNaive = is.na(vLengthNaive) | is.na(vCoverageNaive)
  failureTG = is.na(vLengthTG) | is.na(vCoverageTG)
  cat("Proportion of computational failure POSI",mean(failurePOSI),"\n")
  cat("Proportion of computational failure Naive",mean(failureNaive),"\n")
  cat("Proportion of computational failure TG",mean(failureTG),"\n")
  cat("coverage proportion POSI",mean(vCoveragePOSI[!failurePOSI]),"\n")
  cat("coverage proportion Naive",mean(vCoverageNaive[!failureNaive]),"\n")
  cat("coverage proportion TG",mean(vCoverageTG[!failureTG]),"\n")
  cat("median length POSI",median(vLengthPOSI[!failurePOSI]),"\n")
  cat("median length Naive",median(vLengthNaive[!failureNaive]),"\n")
  cat("median length TG",median(vLengthTG[!failureTG]),"\n")
  cat("quantile 0.9 length POSI",quantile(probs=0.9,vLengthPOSI[!failurePOSI]),"\n")
  cat("quantile 0.9 length Naive",quantile(probs=0.9,vLengthNaive[!failureNaive]),"\n")
  cat("quantile 0.9 length TG",quantile(probs=0.9,vLengthTG[!failureTG]),"\n")
}

