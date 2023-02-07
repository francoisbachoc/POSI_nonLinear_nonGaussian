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

Mres = matrix(nrow=12,ncol=9) 

for (iname in 1:length(vName)) {
  name = vName[iname]
  load(paste0(name,".Rdata"))
  failurePOSI = is.na(vLengthPOSI) | is.na(vCoveragePOSI)
  failureNaive = is.na(vLengthNaive) | is.na(vCoverageNaive)
  failureTG = is.na(vLengthTG) | is.na(vCoverageTG)
  Mres[iname,1] = mean(vCoveragePOSI[!failurePOSI])
  Mres[iname,2] = mean(vCoverageTG[!failureTG])
  Mres[iname,3] = mean(vCoverageNaive[!failureNaive])
  Mres[iname,4] = median(vLengthPOSI[!failurePOSI])
  Mres[iname,5] = median(vLengthTG[!failureTG])
  Mres[iname,6] = median(vLengthNaive[!failureNaive])
  Mres[iname,7] = quantile(probs=0.9,vLengthPOSI[!failurePOSI])
  Mres[iname,8] = quantile(probs=0.9,vLengthTG[!failureTG])
  Mres[iname,9] = quantile(probs=0.9,vLengthNaive[!failureNaive])
}

round(Mres,2)

