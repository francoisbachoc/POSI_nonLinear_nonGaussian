rm( list=ls() )   
source("functions_logistic.r")


##############################
#Specific settings
##############################
set.seed(1)
name="logistic_lasso_pairwiseCorrelated_wellSpecified_zeroBeta_smallLambda_n100"

##############################
#general settings
##############################
n = 100
p = 10
lp=10
lbeta0 = rep(0,10)
nmc=1000
alpha=0.9
I=1000
rho=0.2
lambda = n*0.012
indCoeff=1 #we cover this index in the model
maxL1=20 #We truncate absolute values of hatbeta larger than this
function_X = X_pairwise_correlated
##############################
#Simulations
##############################
vLengthPOSI = matrix(nrow=nmc,ncol=1,data=NaN)
vCoveragePOSI = matrix(nrow=nmc,ncol=1,data=NaN)
vLengthNaive = matrix(nrow=nmc,ncol=1,data=NaN)
vCoverageNaive = matrix(nrow=nmc,ncol=1,data=NaN)
mModel = matrix(nrow=nmc,ncol=p,data=NaN)
vLengthTG = matrix(nrow=nmc,ncol=1,data=NaN)
vCoverageTG = matrix(nrow=nmc,ncol=1,data=NaN)

for (imc in 1:nmc) {
  cat("imc = ",imc,"\n")
  #data generation
  lX = function_X(n,lp,rho)
  X = lX[,1:p]
  theta = h_logistic(lX%*%lbeta0)
  Y = rbinom(n=n,size=1,prob=theta)
  #model selection
  gfit = glmnet(x=X,y=Y,standardize=FALSE,family="binomial",intercept=FALSE)
  beta_hat = coef(gfit, s=lambda/n, exact=TRUE)
  if ( length(beta_hat) == p) { #if no intercept
    M = which(as.matrix(beta_hat)!=0)
  }
  if ( length(beta_hat) == p+1) { #if intercept
    M = which(as.matrix(beta_hat[2:(p+1)])!=0)
  }
  mModel[imc,] = c(M,seq(from=0,to=0,length=p-length(M)))
  itarget = M[indCoeff]
  target = target_logistic(X,M,theta,itarget)
  #POSI inference
  ICPOSItry = tryCatch.W.E(    #call to the POSI procedure with the tryCatch structure
    IC_logistic(X,M,Y,itarget,alpha,I,maxL1)
  )
  ICPOSI = ICPOSItry$value
  if ( is.numeric(ICPOSI$upper)) {
    vLengthPOSI[imc]=ICPOSI$upper - ICPOSI$lower
    vCoveragePOSI[imc] = (ICPOSI$lower <= target) && (target <= ICPOSI$upper)
  }
  #Naive inference
  ICNaivetry = tryCatch.W.E(    #call to the Naive procedure with the tryCatch structure
    IC_logistic_naive(X,M,Y,itarget,alpha,maxL1)
  )
  ICNaive = ICNaivetry$value
  if ( is.numeric(ICNaive$upper)) {
    vLengthNaive[imc]=ICNaive$upper - ICNaive$lower
    vCoverageNaive[imc] = (ICNaive$lower <= target) && (target <= ICNaive$upper)
  }
  #TG inference
  infTGtry = tryCatch.W.E(    #call to the TG procedure with the tryCatch structure
    fixedLassoInf(x=X,y=Y,beta=beta_hat,lambda=lambda,family="binomial",alpha=1-alpha)
  )
  infTG = infTGtry$value
  if ( is.numeric(infTG$ci)) {  #If TG procedure did not return errors
    vLengthTG[imc] = infTG$ci[indCoeff,2] - infTG$ci[indCoeff,1]
    vCoverageTG[imc] = (infTG$ci[indCoeff,1] <= target) && (target <= infTG$ci[indCoeff,2])
  }
}

##############################
#Save of the results
##############################
save_name = paste0(name,".Rdata")
save(file=save_name,X,vCoveragePOSI,vLengthPOSI,vCoverageNaive,vLengthNaive,mModel,vCoverageTG,vLengthTG,n,p,alpha,lp,lbeta0,rho,lambda)

