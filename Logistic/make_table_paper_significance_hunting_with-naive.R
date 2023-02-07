rm( list=ls() )   


vName = c("logistic_significance_hunting_with-naive_zerobeta_nbest20_n30",
          "logistic_significance_hunting_with-naive_zerobeta_nbest20_n100",
          "logistic_significance_hunting_with-naive_non-zerobeta_nbest20_n30",
          "logistic_significance_hunting_with-naive_non-zerobeta_nbest20_n100",
          "logistic_significance_hunting_with-naive_zerobeta_nbest5_n30",
          "logistic_significance_hunting_with-naive_zerobeta_nbest5_n100",
          "logistic_significance_hunting_with-naive_non-zerobeta_nbest5_n30",
          "logistic_significance_hunting_with-naive_non-zerobeta_nbest5_n100")

Mres = matrix(nrow=8,ncol=6)

for (iname in 1:length(vName)) {
  name = vName[iname]
  load(paste0(name,".Rdata"))
  Mres[iname,1] = mean(vCoveragePOSI)
  Mres[iname,2] = mean(vCoverageNaive)
  Mres[iname,3] = median(vLengthPOSI)
  Mres[iname,4] = median(vLengthNaive) 
  Mres[iname,5] = quantile(probs=0.9,vLengthPOSI)
  Mres[iname,6] = quantile(probs=0.9,vLengthNaive)
}

round(Mres,2)
