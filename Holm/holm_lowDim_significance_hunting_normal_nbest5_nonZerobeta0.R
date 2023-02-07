rm( list=ls() )   
source("functions_holm.r")

##############################
#Specific settings
##############################
set.seed(1)
noise_gen = function(n,sigma) {
  sigma*rnorm(n)
}
generate_X = function(n,p) {
  X_table1_uniform(n,p)
}
name="holm_lowDim_significance_hunting_normal_nbest5_nonZeroBeta0"

##############################
#general settings
##############################
n = 100
p = 5
sigma = 1
beta_0 = c(2,-1,0,0,1)
nmc=1000
alpha=0.9
I=1000
nbest = 5
lambda = 2

##############################
#Simulations
##############################
vLengthPOSI = rep(-1,nmc)
vCoveragePOSI = rep(-1,nmc)
vLengthNaive = rep(-1,nmc)
vCoverageNaive = rep(-1,nmc)
mM = matrix(nrow=nmc,ncol=p,data=0)
vitarget = rep(-1,nmc)
vHatSigma = rep(-1,nmc)

for (i in 1:nmc) {
  cat("i = ",i,"\n")
  #data generation
  X = generate_X(n,p)
  theta = X %*% beta_0
  Y = theta + noise_gen(n,sigma)
  #model selection
  sig_hun = significance_hunting(X,Y,lambda,nbest)
  M = sig_hun$M
  mM[i,1:sum(M>0.5)] = M
  itarget = sig_hun$i
  vitarget[i] = itarget
  #specifying the target
  target = target_holm(X,M,theta,itarget)
  #POSI inference
  ICPOSI = IC_holm(X,M,Y,itarget,alpha,I)
  vLengthPOSI[i]=ICPOSI$upper - ICPOSI$lower
  vCoveragePOSI[i] = (ICPOSI$lower <= target) && (target <= ICPOSI$upper)
  vHatSigma[i] = hat_sigma_holm(X,M,Y)
  #Naive inference
  ICNaive = IC_naive_holm(X,M,Y,itarget,alpha)
  vLengthNaive[i]=ICNaive$upper - ICNaive$lower
  vCoverageNaive[i] = (ICNaive$lower <= target) && (target <= ICNaive$upper)
}

##############################
#Save of the results
##############################
save_name = paste0(name,".Rdata")
save(file=save_name,n,p,nbest,beta_0,lambda,X,vCoveragePOSI,vLengthPOSI,
     vCoverageNaive,vLengthNaive,mM,vHatSigma,vitarget)

