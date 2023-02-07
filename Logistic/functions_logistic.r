library(lars)
library(selectiveInference)
library(LaplacesDemon)
library(fGarch)
library(rmutil)
library(pracma)


subsets=function(n,r,v=1:n){
  # subsects() gets all subsets of size r from a set of size n
  if(r<=0) NULL else
    if(r>=n) v[1:n] else
      rbind(cbind(v[1],subsets(n-1,r-1,v[-1])),subsets(n-1,r,v[-1]))
}

alltbMj = function(X) {
  #X...........: the nXp design matrix (full rank, p <=n) 
  #return......: the matrix of size (p2^(p-1)) X p  of all the vectors bar{t}_M (normalized). 
  #
  n = dim(X)[1]
  p = dim(X)[2]
  mtMj = matrix(data=0,nrow=p*2^(p-1),ncol=n)
  index = 1 # index of the next tM to be calculated
  for (m.sz in 1:p) {   #loop on size of models
    if (m.sz<p){ #obtaining the matrix of all the models of a given size
      models.of.sz <- subsets(p,m.sz)
    }else{
      models.of.sz <- matrix(1:p,1,p)
    }
    for (i in 1:dim(models.of.sz)[1]) {   #loop on models of a given size
      model = models.of.sz[i,]
      for (j in 1:m.sz) {   #loop on regressors in a given model
        XM = as.matrix(X[,model])
        mtMj[index,] = (solve( t(XM)%*%XM ) %*% t(XM))[j,]  #computation of t_M for given j,M
        index = index+1
      }
    }
  }
  vNorm = sqrt(rowSums(mtMj^2)) #normalization of the non-zero vectors
  mtMj[vNorm!=0,] = mtMj[vNorm!=0,] / (matrix(data=vNorm[vNorm!=0],nrow=length(vNorm[vNorm!=0]),ncol=1)%*%matrix(data=1,nrow=1,ncol=n)) 
  return(mtMj)
}

maxOfInnerProducs = function(X,I) {
  #X...........: the nXp design matrix (full rank, p <=n)   
  #I...........: number of Monte Carlo sample
  #return......: vector of I realizations of max_{M,j} ( tMj' V ). V standard n-dim gaussian 
  #
  n = dim(X)[1]
  mbtMj = alltbMj(X) #matrix of all the tbMj over all submodels and coefficients
  U = matrix(nrow = n,ncol=I,data=rnorm(n*I))  #matrix of I realizations of N(0,I_n)
  return( apply(abs(mbtMj%*%U),2,max) )  #the vectors of the I max
}

Kone = function(X,alpha,I=1000) {
  #X...........: the nXp design matrix (full rank)  
  #alpha.......: confidence level
  #I...........: number of Monte Carlo sample for evaluating proba of excedeence of max of inner products
  #return......: Constant K1, estimated by Monte Carlo
  #
  n = dim(X)[1]
  vC = maxOfInnerProducs(X,I) #I realizations of maximum of inner product
  quantile(vC,alpha)
}

target_logistic = function(X,M,theta,j) {
  #Compute the pseudo target in logistic regression
  #X: n*p design matrix
  #M: the submodel (subset of 1:p)
  #theta: n*1 vector of the true P(y_i=1)
  #j in M: the coefficient index (in the full model) under consideration
  XM = X[,M]
  glmfit = glm.fit(x=XM, y=theta, family = binomial(), intercept = FALSE)
  hatbeta = glmfit$coefficients
  jdM = which(M==j)
  hatbeta[jdM]
}

hat_beta_logistic = function(X,M,Y,maxL1=10^6) {
  #Compute the restricted ML estimator in logistic regression
  #X: n*p design matrix
  #M: the submodel (subset of 1:p)
  #Y: n*1 vector of binary observations
  #maxL1: if the absolute L1 norm of an estimated coefficient. We projec the estimation on the
  #       L1 ball
  XM = X[,M]
  glmfit = glm.fit(x=XM, y=Y, family = binomial(), intercept = FALSE)
  hatbeta = glmfit$coefficients
  if (sum(abs(hatbeta)) > maxL1) {
    hatbeta = (maxL1/sum(abs(hatbeta))) *  hatbeta 
  } 
  hatbeta
}

st_dev_logistic = function(X,M,Y,j,maxL1=10^6) {
  #compute the upper bound of the asymptotic standard deviation of the maximum likelihood estimator of a pseudo-coefficient
  #X: n*p design matrix
  #M: submodel, subset of 1:p
  #Y: n*1 vector of binary observations
  #j: index of the coefficient (in the full model)
  #maxL1: see hat_beta_logistic
  jdM = which(M==j)
  XM = X[,M]
  betaM = hat_beta_logistic(X,M,Y,maxL1)
  gammaM = as.matrix(XM)%*%betaM
  D = diag( as.vector(-Y*phi_1_second_logistic(gammaM) - (1-Y)*phi_2_second_logistic(gammaM)))
  HM = t(XM)%*%D%*%XM
  iHM = solve(HM)
  uM = ( h_prime_logistic(gammaM) / ( h_logistic(gammaM)*(1-h_logistic(gammaM)) ) ) * ( Y - h_logistic(gammaM) )
  SM = iHM%*%t(XM)%*%diag(as.vector(uM^2))%*%XM%*%iHM
  sqrt(SM[jdM,jdM])
}

st_dev_logistic_naive = function(X,M,Y,j,maxL1=10^6) {
  #compute the estimate of the asymptotic standard deviation of the maximum likelihood estimator of a pseudo-coefficient
  #X: n*p design matrix
  #M: submodel, subset of 1:p
  #Y: n*1 vector of binary observations
  #j: index of the coefficient (in the full model)
  #maxL1: see hat_beta_logistic
  jdM = which(M==j)
  XM = X[,M]
  betaM = hat_beta_logistic(X,M,Y,maxL1)
  gammaM = as.matrix(XM)%*%betaM
  D = diag( as.vector(-Y*phi_1_second_logistic(gammaM) - (1-Y)*phi_2_second_logistic(gammaM)))
  HM = t(XM)%*%D%*%XM
  iHM = solve(HM)
  ubM2 = h_prime_logistic(gammaM)^2 / ( h_logistic(gammaM)*(1-h_logistic(gammaM)))
  SM = iHM%*%t(XM)%*%diag(as.vector(ubM2))%*%XM%*%iHM
  sqrt(SM[jdM,jdM])
}

asymptotic_cov_matrix_logistic_naive = function(X,Y,maxL1=10^6) {
  #compute the estimate of the asymptotic covariance matrix
  #of all the estimation errors
  #X: n*p design matrix
  #Y: n*1 vector of binary observations
  #maxL1: see hat_beta_logistic
  #
  #precompute empty asymptotic covariance matrix
  ksize = p*2^(p-1)
  asym_cov = matrix(nrow=k_,ncol=ksize,data=NaN)
  #Generate all submodels
  mM = matrix(nrow=2^p-1,ncol=p,data=0)
  cpt=0
  for (k in 1:p) {    #loop on submodel size
    mM[(cpt+1):(cpt+nchoosek(p,k)),1:k] = subsets(p,k)
    cpt = cpt+nchoosek(p,k)
  }
  #compute vector of total sizes
  #(number of pseudo parameters up to each model)
  vTotalSize = rep(-1,2^p-1)
  vTotalSize[1] = 0
  for (i in 2:dim(mM)[1]) {
    vTotalSize[i] = vTotalSize[i-1] + sum(mM[i-1,] >= 0.5) 
  }
  #compute cov matrix of Y
  hatbeta = hat_beta_logistic(X,1:p,Y,maxL1)
  hatgamma = as.matrix(X)%*%hatbeta
  cov_y = diag(as.vector(h_logistic(hatgamma)*(1-h_logistic(hatgamma))))
  #fill in asymptotic covariance matrix
  for (im1 in 1:dim(mM)[1]) {
    for (im2 in 1:dim(mM)[1]) {
      M1 = mM[im1,]
      M2 = mM[im2,]
      XM1 = X[,M1]
      XM2 = X[,M2]
      hatbetaM1 = hat_beta_logistic(X,M1,Y,maxL1)
      hatgammaM1 = as.matrix(XM1)%*%hatbetaM1
      hatbetaM2 = hat_beta_logistic(X,M2,Y,maxL1)
      hatgammaM2 = as.matrix(XM2)%*%hatbetaM2
      D1 = diag( as.vector(-Y*phi_1_second_logistic(hatgammaM1) - (1-Y)*phi_2_second_logistic(hatgammaM1)))
      D2 = diag( as.vector(-Y*phi_1_second_logistic(hatgammaM2) - (1-Y)*phi_2_second_logistic(hatgammaM2)))
      HM1 = t(XM1)%*%D1%*%XM1
      HM2 = t(XM2)%*%D2%*%XM2
      iHM1 = solve(HM1)
      iHM2 = solve(HM2)
      indexCoeffM1 = (vTotalSize[im1]+1):(vTotalSize[im1]+sum(M1>=0.5))
      indexCoeffM2 = (vTotalSize[im2]+1):(vTotalSize[im2]+sum(M2>=0.5))
      asym_cov[indexCoeffM1,indexCoeffM2] = iHM1%*%t(XM1)%*%cov_y%*%XM2%*%iHM2
    }
  }
  asym_cov
}

KSigma = function(Sigma,alpha,I,nugget=10^(-6)) {
  #Compute the alpha quantile of the max absolute value of Z
  #where Z sim N(0 , Sigma)
  d = dim(Sigma)[1]
  mz = matrix(nrow=I,ncol=d,data=rnorm(I*d))%*%chol(Sigma+nugget*min(diag(Sigma))*diag(d))
  vmax = apply(abs(mz),1,'max')
  quantile(vmax,alpha)
}

get_corr_matrix = function(Sigma) {
  #get a correlation matrix from a covariance matrix
  corrM = Sigma
  for (i in 1:dim(Sigma)[1]) {
    for (j in 1:dim(Sigma)[2]) {
      corrM[i,j] = Sigma[i,j] / (sqrt(Sigma[i,i]*Sigma[j,j]  ))
    }
  }
  corrM
}

Balpha = function(q,N,alpha,I=1000) {
  #compute the upper bound for quantile of inner products between standard Gaussian vector and unit vectors
  #cf K4 in Bachoc Leep Poetscher 2014 and Balpha in Bachoc Preinerstorfer Steinberger 2016
  #q...........: dimension of the Gaussian vector
  #N...........: number of unit vectors
  #alpha.......: confidence level
  #I...........: numerical precision
  #return......: upper bound, numerically approximated
  #
  vC = sqrt(qbeta(p=seq(from=0,to=1/N,length=I),shape1=1/2,shape2=(q-1)/2, lower.tail=FALSE)) #vector of quantiles of Beta distribution
  fconfidence = function(K){mean(pchisq(q=(K/vC)^2,df=q,lower.tail=FALSE))-(1-alpha)} #MC evaluation of confidence level for a constant K
  uniroot(fconfidence,interval=c(1,100))$root 
}

IC_logistic = function(X,M,Y,j,alpha,I,maxL1=10^6) {
  #POSI Confidence interval for logistic regression with all submodels and coefficients allowed
  #X...........: the nXp design matrix (full rank)  
  #M: the submodel (subset of 1:p)
  #Y: nx1 response vector
  #j in M: the coefficient index (in the full model) under consideration
  #alpha: nominal coverage probability
  #I: number of Monte Carlo samples to compute 
  #maxL1: see hat_beta_logistic
  #return: lower and upper bound for the confidence interval
  p = dim(X)[2]
  N = p*2^(p-1)
  jdM = which(M==j)
  center = as.vector(hat_beta_logistic(X,M,Y,maxL1))[jdM]
  Ba = Balpha(q=p,N=N,alpha=alpha,I=1000)
  radius = st_dev_logistic(X,M,Y,j,maxL1)*Ba
  list(lower=center-radius,upper=center+radius)
}

IC_plug_in_logistic = function(X,M,Y,j,alpha,I,maxL1=10^6,nugget=10^(-6)) {
  #POSI Confidence interval for logistic regression with all submodels and coefficients allowed
  #plug in of estimates in the asymptotic covariance matrix of the
  #covariance matrix of all the estimators
  #X...........: the nXp design matrix (full rank)  
  #M: the submodel (subset of 1:p)
  #Y: nx1 response vector
  #j in M: the coefficient index (in the full model) under consideration
  #alpha: nominal coverage probability
  #I: number of Monte Carlo samples to compute 
  #maxL1: see hat_beta_logistic
  #nugget: see KSigma
  #return: lower and upper bound for the confidence interval
  p = dim(X)[2]
  N = p*2^(p-1)
  jdM = which(M==j)
  center = as.vector(hat_beta_logistic(X,M,Y,maxL1))[jdM]
  asym_cov_mat = asymptotic_cov_matrix_logistic_naive(X,Y,maxL1)
  Kposi = KSigma(Sigma = get_corr_matrix(asym_cov_mat),alpha,I)
  radius = st_dev_logistic(X,M,Y,j,maxL1)*Kposi
  list(lower=center-radius,upper=center+radius)
}

IC_logistic_naive = function(X,M,Y,j,alpha,maxL1=10^6) {
  #naive Confidence interval for logistic regression
  #(takes quantile of plug-in-estimated asymptotic
  #distribution assuming selected model is fixed)
  #X...........: the nXp design matrix (full rank)  
  #M: the submodel (subset of 1:p)
  #Y: nx1 response vector
  #j in M: the coefficient index (in the full model) under consideration
  #alpha: nominal coverage probability
  #I: number of Monte Carlo samples to compute 
  #maxL1: see hat_beta_logistic
  #return: lower and upper bound for the confidence interval
  p = dim(X)[2]
  N = p*2^(p-1)
  jdM = which(M==j)
  center = as.vector(hat_beta_logistic(X,M,Y,maxL1))[jdM]
  radius = st_dev_logistic_naive(X,M,Y,j,maxL1)*qnorm(p=1-(1-alpha)/2)
  list(lower=center-radius,upper=center+radius)
}

h_logistic = function(gamma) {
  exp(gamma)/(1+exp(gamma))
}

h_prime_logistic = function(gamma) {
  exp(gamma)/((1+exp(gamma))^2)
}

phi_1_logistic = function(gamma) {
  log(h_logistic(gamma))
}

phi_2_logistic = function(gamma) {
  log(1-h_logistic(gamma))
}

phi_1_second_logistic = function(gamma) {
  -exp(gamma)/((1+exp(gamma))^2)
}

phi_2_second_logistic = function(gamma) {
  -exp(gamma)/((1+exp(gamma))^2)
}

X_pairwise_correlated = function(n,p,rho) {
  #  Generate a n*p matrix where each line is a Gaussian vector with variances 1 and
  # pairwise correlation rho
  Sigma = matrix(nrow=p,ncol=p,data=rho) + (1-rho) * diag(p)
  X = matrix(nrow=n,ncol=p,data=rnorm(n*p))
  X = X%*%chol(Sigma)
}

X_equi_correlated = function(n,p,rho = sqrt(1/p)) {
  #Generate a n*p matrix where each line is a Gaussian vector with variances 1,
  #independence of the p-1 first components, and correlation rho
  #between component p and the others
  Sigma = diag(p)
  Sigma[1:(p-1),p] = rho
  Sigma[p,1:(p-1)] = rho
  X = matrix(nrow=n,ncol=p,data=rnorm(n*p))
  X = X%*%chol(Sigma)
}

tryCatch.W.E <- function(expr) {
  W <- NULL
  w.handler <- function(w)  { # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}

significance_hunting_logistic = function(X,Y,lambda,nbest,maxL1=10^6) {
  #Does the significance hunting procedure suggested by David
  #Rank all possible models with respect to penalized log lokelihood
  #with penalty lambda*card(model)
  #Then for the nbest models, pick the one with largest absolute asymptotic t-statistic
  #return the corresponding model and coefficient
  #
  #INPUTS
  #X:        n*p design matrix
  #Y:        n*1 response vector
  #lambda:   penalty parameter
  #nbest:    nbest of best models among which t-statistics is searched for
  #maxL1:    see hat_beta_logistic
  #
  #OUTPUTS
  #M: selected model
  #i: corresponding coefficient (w.r.t. full model) achieving largest absolute t-statistic
  #
  #Generating all submodels
  mM = matrix(nrow=2^p-1,ncol=p,data=0)
  cpt=0
  for (k in 1:p) {    #loop on submodel size
    mM[(cpt+1):(cpt+nchoosek(p,k)),1:k] = subsets(p,k)
    cpt = cpt+nchoosek(p,k)
  }
  vPenLik = rep(0,dim(mM)[1])
  #
  #Computing penalized likelihood for all submodels
  for (i in 1:length(vPenLik)) {  #loop on all submodels
    M = mM[i,]
    XM = X[,M]
    hatbetaM = hat_beta_logistic(X,M,Y,maxL1)
    hatgammaM = as.matrix(XM)%*%hatbetaM
    vPenLik[i] = - sum(Y*phi_1_logistic(hatgammaM)) - sum((1-Y)*phi_2_logistic(hatgammaM)) + lambda * sum(M>=0.5)
  }
  #
  #Comuting T statistics for submodels with largest penalized likelihood
  mMbest = mM[sort.int(vPenLik,index.return = TRUE)$ix[1:nbest],]
  mTbest = matrix(nrow=nbest,ncol=p,data=-1) 
  for (iM in 1:nbest) {   #loop on the best submodels
    M = mMbest[iM,]
    hat_betaM = hat_beta_logistic(X,M,Y,maxL1)
    for (ic in M[M>0.5]) {  #loop on variable in submodel
      st_dev_ = st_dev_logistic_naive(X,M,Y,ic,maxL1)
      mTbest[iM,ic] = abs(hat_betaM[which(M==ic)])/st_dev_
    } 
  }
  vMaxAbsT = apply(mTbest,1,'max')
  imax = which.max(vMaxAbsT)
  M = mMbest[imax,]
  M = M[which(M>0.5)]
  list(M = M,i=which.max(mTbest[imax,]))
}
