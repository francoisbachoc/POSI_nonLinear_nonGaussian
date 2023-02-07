rm( list=ls() )   


vName = c("logistic_significance_hunting_with-naive_zerobeta_nbest20_n30",
          "logistic_significance_hunting_with-naive_zerobeta_nbest5_n30",
          "logistic_significance_hunting_with-naive_non-zerobeta_nbest20_n30",
          "logistic_significance_hunting_with-naive_non-zerobeta_nbest5_n30",
          "logistic_significance_hunting_with-naive_zerobeta_nbest20_n100",
          "logistic_significance_hunting_with-naive_zerobeta_nbest5_n100",
          "logistic_significance_hunting_with-naive_non-zerobeta_nbest20_n100",
          "logistic_significance_hunting_with-naive_non-zerobeta_nbest5_n100")

for (name in vName) {
  load(paste0(name,".Rdata"))
  cat("############################## \n")
  cat(name,"\n")
  cat("############################## \n")
  cat("n = ",n,"\n")
  cat("p = ",p,"\n")
  cat("alpha = ",p,"\n")
  cat("nbest = ",nbest,"\n")
  cat("lambda = ",lambda,"\n")
  cat("rho = ",rho,"\n")
  cat("coverage proportion POSI",mean(vCoveragePOSI),"\n")
  cat("median length POSI",median(vLengthPOSI),"\n")
  cat("quantile 0.9 length POSI",quantile(probs=0.9,vLengthPOSI),"\n")
  cat("coverage proportion Naive",mean(vCoverageNaive),"\n")
  cat("median length Naive",median(vLengthNaive),"\n")
  cat("quantile 0.9 length Naive",quantile(probs=0.9,vLengthNaive),"\n")
  cat("############################## \n")
}


