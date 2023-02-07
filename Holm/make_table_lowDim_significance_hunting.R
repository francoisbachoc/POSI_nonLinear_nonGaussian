rm( list=ls() )   


vName = c("holm_lowDim_significance_hunting_normal_nbest20",
          "holm_lowDim_significance_hunting_normal_nbest20_nonZeroBeta0",
          "holm_lowDim_significance_hunting_normal_nbest5",
          "holm_lowDim_significance_hunting_normal_nbest5_nonZeroBeta0")

Mres = matrix(nrow=length(vName),ncol=6)

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
