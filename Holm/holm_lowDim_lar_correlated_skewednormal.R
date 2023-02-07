rm( list=ls() )   
source("functions_holm.r")

##############################
#Specific settings
##############################
set.seed(1)
noise_gen = function(n,sigma) {
  rsnorm(n=n, mean = 0, sd = sigma, xi = 5)
}
name="holm_lowDim_lar_correlated_skewednormal"
generate_X = function(n,p) {
  X_correlated(n,p,rho)
}

##############################
#general settings
##############################
n = 50
p = 10
sigma = 1
beta_0 = c(-4,4,0*(1:(p-2)))
nmc=500
alpha=0.9
nStep=3
I=1000
rho=0.1

##############################
#Simulations
##############################
mLengthPOSI = matrix(nrow=nmc,ncol=nStep,data=-1)
mCoveragePOSI = matrix(nrow=nmc,ncol=nStep,data=-1)
mCoverageAllPOSI = matrix(nrow=nmc,ncol=nStep,data=-1)
mLengthTG = matrix(nrow=nmc,ncol=nStep,data=-1)
mCoverageTG = matrix(nrow=nmc,ncol=nStep,data=-1)
mCoverageAllTG = matrix(nrow=nmc,ncol=nStep,data=-1)
mLengthNaive = matrix(nrow=nmc,ncol=nStep,data=-1)
mCoverageNaive = matrix(nrow=nmc,ncol=nStep,data=-1)
mCoverageAllNaive = matrix(nrow=nmc,ncol=nStep,data=-1)
mMstep = matrix(nrow=nmc,ncol=nStep,data=-1)
mHatSigma = matrix(nrow=nmc,ncol=nStep,data=-1)

for (i in 1:nmc) {
  cat("i = ",i,"\n")
  X = generate_X(n,p)
  theta = X %*% beta_0
  Y = theta + noise_gen(n,sigma)
  lar = lar(X,Y,intercept=FALSE)
  Msteps = lar$action
  mMstep[i,] = Msteps[1:nStep]
  #The 3 steps
  for (k in 1:nStep) {
    #specifying the target
    M = Msteps[1:k]
    itarget = M[k]
    target = target_holm(X,M,theta,itarget)
    #POSI inference
    ICPOSI = IC_holm(X,M,Y,itarget,alpha,I)
    mLengthPOSI[i,k]=ICPOSI$upper - ICPOSI$lower
    mCoveragePOSI[i,k] = (ICPOSI$lower <= target) && (target <= ICPOSI$upper)
    mHatSigma[i,k] = hat_sigma_holm(X,M,Y)
    #TG inference
    larinf = larInf(lar,alpha=1-alpha,k=k,type='active',gridrange = c(-500,500))
    ICTGlower = larinf$ci[k,1]
    ICTGupper = larinf$ci[k,2]
    mLengthTG[i,k]= ICTGupper - ICTGlower
    mCoverageTG[i,k] = (ICTGlower <= target) && (target <= ICTGupper)
    #Naive inference
    ICNaive = IC_naive_holm(X,M,Y,itarget,alpha)
    mLengthNaive[i,k]=ICNaive$upper - ICNaive$lower
    mCoverageNaive[i,k] = (ICNaive$lower <= target) && (target <= ICNaive$upper)
  }
  #The simultaneous coverage
  M = Msteps[1:nStep]
  larinf = larInf(lar,alpha=1-alpha,k=nStep,type='all',gridrange = c(-500,500))
  for (k in 1:nStep) {
    #specifying the target
    itarget = M[k]
    target = target_holm(X,M,theta,itarget)
    #POSI inference
    ICPOSI = IC_holm(X,M,Y,itarget,alpha,I)
    mCoverageAllPOSI[i,k] = (ICPOSI$lower <= target) && (target <= ICPOSI$upper)
    #TG inference
    ICTGlower = larinf$ci[k,1]
    ICTGupper = larinf$ci[k,2]
    mCoverageAllTG[i,k] = (ICTGlower <= target) && (target <= ICTGupper)
    #Naive inference
    ICNaive = IC_naive_holm(X,M,Y,itarget,alpha)
    mCoverageAllNaive[i,k] = (ICNaive$lower <= target) && (target <= ICNaive$upper)
  }
}

##############################
#Save of the results
##############################
save_name = paste0(name,".Rdata")
save(file=save_name,X,mCoveragePOSI,mLengthPOSI,mCoverageAllPOSI,
     mCoverageTG,mLengthTG,mCoverageAllTG,
     mCoverageNaive,mLengthNaive,mCoverageAllNaive,
     mMstep,mHatSigma)

