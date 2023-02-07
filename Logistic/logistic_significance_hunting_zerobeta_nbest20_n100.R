rm( list=ls() )   
source("functions_logistic.r")

##############################
#Specific settings
##############################
set.seed(1)
name="logistic_significance_hunting_zerobeta_nbest20_n100"

##############################
#general settings
##############################
n = 100
p = 5
beta_0 = rep(0,p)
nmc=1000
alpha=0.9
I=10000
lambda = 2
nbest = 20
maxL1=20
rho = 0.8

##############################
#Simulations
##############################
vLengthPOSI = rep(-1,nmc)
vCoveragePOSI = rep(-1,nmc)
mM = matrix(nrow=nmc,ncol=p,data=0)
vitarget = rep(-1,nmc)


for (i in 1:nmc) {
  cat("i = ",i,"\n")
  #data generation
  X = X_pairwise_correlated(n,p,rho)
  theta = h_logistic(X%*%beta_0)
  Y = rbinom(n=n,size=1,prob=theta)
  #model selection
  sig_hun = significance_hunting_logistic(X,Y,lambda,nbest,maxL1)
  M = sig_hun$M
  mM[i,1:sum(M>0.5)] = M
  itarget = sig_hun$i
  vitarget[i] = itarget
  #specifying the target
  target = target_logistic(X,M,theta,itarget)
  #POSI inference
  ICPOSI = IC_logistic(X,M,Y,itarget,alpha,I,maxL1)
  vLengthPOSI[i]=ICPOSI$upper - ICPOSI$lower
  vCoveragePOSI[i] = (ICPOSI$lower <= target) && (target <= ICPOSI$upper)
}

##############################
#Save of the results
##############################
save_name = paste0(name,".Rdata")
save(file=save_name,n,p,alpha,nbest,beta_0,rho,lambda,vCoveragePOSI,vLengthPOSI,mM,vitarget)

