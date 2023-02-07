rm( list=ls() )   


vName = c("independent_normal",
          "independent_Laplace",
          "independent_uniform",
          "independent_skewednormal",
          "correlated_normal",
          "correlated_Laplace",
          "correlated_uniform",
          "correlated_skewednormal")

#by row: for each 8 settings: posi, TG and Naive
#by column: for the three step: coverage, median length, q0.9 length
#Then simultaneous coverage
Mres = matrix(nrow=3*length(vName),ncol = 10) 
v_explain_line = rep('expl',3*length(vName))

for (iname in 1:length(vName)) {
  name = vName[iname]
  load(paste0("holm_highDim_lar_",name,".Rdata"))
  for (i in 1:3) {
    Mres[3*(iname-1)+1,3*(i-1)+1] = mean(mCoveragePOSI[,i])
    Mres[3*(iname-1)+1,3*(i-1)+2] = median(mLengthPOSI[,i])
    Mres[3*(iname-1)+1,3*(i-1)+3] = quantile(probs=0.9,mLengthPOSI[,i])
  }
  v_explain_line[3*(iname-1)+1] = paste0(name,'_POSI')
  for (i in 1:3) {
    Mres[3*(iname-1)+2,3*(i-1)+1] = mean(mCoverageTG[,i])
    Mres[3*(iname-1)+2,3*(i-1)+2] = median(mLengthTG[,i])
    Mres[3*(iname-1)+2,3*(i-1)+3] = quantile(probs=0.9,mLengthTG[,i])
  }
  v_explain_line[3*(iname-1)+2] = paste0(name,'_TG')
  for (i in 1:3) {
    Mres[3*(iname-1)+3,3*(i-1)+1] = mean(mCoverageNaive[,i])
    Mres[3*(iname-1)+3,3*(i-1)+2] = median(mLengthNaive[,i])
    Mres[3*(iname-1)+3,3*(i-1)+3] = quantile(probs=0.9,mLengthNaive[,i])
  }
  v_explain_line[3*(iname-1)+3] = paste0(name,'_Naive')
  Mres[3*(iname-1)+1,10] = mean(rowSums(mCoveragePOSI)==3)
  Mres[3*(iname-1)+2,10] = mean(rowSums(mCoverageTG)==3)
  Mres[3*(iname-1)+3,10] = mean(rowSums(mCoverageNaive)==3)
}

M = data.frame(round(Mres,2),row.names = v_explain_line)
Mlatex = matrix(nrow=24,ncol=21,data="&")
Mlatex[,1] = v_explain_line
for (i in 1:10) {
  Mlatex[,3+2*(i-1)] = M[,i]
}
noquote(Mlatex)
round(Mres,2)
