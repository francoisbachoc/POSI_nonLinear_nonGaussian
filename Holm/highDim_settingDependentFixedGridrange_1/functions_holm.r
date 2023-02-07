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


hat_beta_holm = function(X,M,Y) { 
  #X...........: the nXp design matrix (full rank)  
  #M: the submodel (subset of 1:p)
  #Y: nx1 response vector
  #return: least square estimator in the selected model
  XM = X[,M]
  solve(t(XM)%*%XM)%*%t(XM)%*%Y
}

target_holm = function(X,M,theta,j) { 
  #X...........: the nXp design matrix (full rank)  
  #M: the submodel (subset of 1:p)
  #theta: nx1 unknown mean
  #j in M: the coefficient index (in the full model) under consideration
  #return: least square estimator in the selected model
  XM = X[,M]
  jdM = which(M==j)
  (solve(t(XM)%*%XM)%*%t(XM)%*%theta)[jdM]
}

hat_sigma_holm = function(X,M,Y) { 
  #X...........: the nXp design matrix (full rank)  
  #M: the submodel (subset of 1:p)
  #Y: nx1 response vector
  #return: standard deviation estimator in the selected model
  XM = X[,M]
  sqrt( (1/(n-length(M))) * sum(  (Y-XM%*%hat_beta_holm(X,M,Y))^2  ) )
}

st_dev_holm = function(X,M,j) {
  #X...........: the nXp design matrix (full rank)  
  #M: the submodel (subset of 1:p)
  #j in M: the coefficient index (in the full model) under consideration
  #return: lower and upper bound for the confidence interval
  jdM = which(M==j)
  XM = X[,M]
  sqrt(solve(t(XM)%*%XM)[jdM,jdM])
}

IC_holm = function(X,M,Y,j,alpha,I) {
  #POSI Confidence interval in the homoscedastic model with all submodels and coefficients allowed
  #X...........: the nXp design matrix (full rank)  
  #M: the submodel (subset of 1:p)
  #Y: nx1 response vector
  #j in M: the coefficient index (in the full model) under consideration
  #alpha: nominal coverage probability
  #I: number of MC samples to compute K1
  #return: lower and upper bound for the confidence interval
  jdM = which(M==j)
  center = as.vector(hat_beta_holm(X,M,Y))[jdM]
  radius = hat_sigma_holm(X,M,Y)*st_dev_holm(X,M,j)*Kone(X,alpha,I)
  list(lower=center-radius,upper=center+radius)
}

IC_s_holm = function(X,M,Y,j,alpha,I,s) {
  #POSI Confidence interval in the homoscedastic model
  #   with submodels up to size s allowed and all coefficients
  #   uses the upper bound Balpha on the POSI constant
  #X...........: the nXp design matrix (all submatrices with less than s columns have full rank)  
  #M: the submodel (subset of 1:p of size less than s)
  #Y: nx1 response vector
  #j in M: the coefficient index (in the full model) under consideration
  #alpha: nominal coverage probability
  #I: number of MC samples to compute Balpha
  #return: lower and upper bound for the confidence interval
  N = 0
  for (k in 1:s) {
    N = N + k*nchoosek(p,k)
  }
  Kupper = Balpha(q=n,N=N,alpha=alpha,I=I) 
  jdM = which(M==j)
  center = as.vector(hat_beta_holm(X,M,Y))[jdM]
  radius = hat_sigma_holm(X,M,Y)*st_dev_holm(X,M,j)*Kupper
  list(lower=center-radius,upper=center+radius)
}

IC_naive_holm = function(X,M,Y,j,alpha) {
  #naive Confidence interval in the homoscedastic model 
  #X...........: the nXp design matrix (full rank)  
  #M: the submodel (subset of 1:p)
  #Y: nx1 response vector
  #j in M: the coefficient index (in the full model) under consideration
  #alpha: nominal coverage probability
  #return: lower and upper bound for the confidence interval
  jdM = which(M==j)
  center = as.vector(hat_beta_holm(X,M,Y))[jdM]
  radius = hat_sigma_holm(X,M,Y)*st_dev_holm(X,M,j)*qnorm(p=(1-alpha)/2,lower.tail=FALSE)
  list(lower=center-radius,upper=center+radius)
}

X_table1_uniform = function(n,p) {
  #generate the matrix X as described in Table 1 of Uniform asymptotic inference and the bootstrap
  #     after model selection
  X = matrix(nrow=n,ncol=p)
  for (i in 1:p) {
    ind = sample(1:3,1)
    if (ind == 1) {
      X[,i] = rnorm(n)
    }
    if (ind == 2) {
      X[,i] = rbern(n=n,prob=0.5)
    } 
    if (ind == 3) {
      X[,i] = rsnorm(n=n, mean = 0, sd = 1, xi = 5)
    }
  }
  X = X / (matrix (nrow=n,ncol=1,data=1) %*% matrix(nrow=1,ncol=p,data=sqrt(colSums(X^2))) )
  X
}

X_correlated = function(n,p,rho) {
  #Generate a random matrix where each line is independent and Gaussian with AR(1) structure
  #rho: correlation (exp(-rho|i-j|))
  #each column is finally renormalized to have unit norm
  Sigma = exp(-rho*abs(outer(seq(1,p),seq(1,p),'-')))
  X = matrix(nrow=n,ncol=p,data=rnorm(n*p))
  X = X%*%chol(Sigma)
  X = X / (matrix (nrow=n,ncol=1,data=1) %*% matrix(nrow=1,ncol=p,data=sqrt(colSums(X^2))) )
  X
}

significance_hunting = function(X,Y,lambda,nbest) {
  #Does the significance hunting procedure suggested by David
  #Rank all possible models with respect to penalized log lokelihood
  #with penalty lambda*card(model)
  #Then for the nbest models, pick the one with largest absolute t-statistic
  #return the corresponding model and coefficient
  #
  #INPUTS
  #X:        n*p design matrix
  #Y:        n*1 response vector
  #lambda:   penalty parameter
  #nbest:    nbest of best models among which t-statistics is searched for
  #
  #OUTPUTS
  #M: selecteed model
  #i: corresponding coefficient (w.r.t. full model) achieving largest absolute t-statistic
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
    hatsigma = hat_sigma_holm(X,M,Y)
    hatbeta = hat_beta_holm(X,M,Y)
    vPenLik[i] = log(hatsigma)+(1/2)*(sum((Y-X[,M]%*%hatbeta)^2)/hatsigma^2) + lambda * sum(M>=0.5)
  }
  #
  #Comuting T statistics for submodels with largest penalized likelihood
  mMbest = mM[sort.int(vPenLik,index.return = TRUE)$ix[1:nbest],]
  mTbest = matrix(nrow=nbest,ncol=p,data=-1) 
  for (iM in 1:nbest) {   #loop on the best submodels
    M = mMbest[iM,]
    hatbeta = hat_beta_holm(X,M,Y)
    hatsigma = hat_sigma_holm(X,M,Y)
    for (ic in M[M>0.5]) {  #loop on variable in submodel
      st_dev_ = st_dev_holm(X,M,ic)
      mTbest[iM,ic] = abs(hatbeta[which(M==ic)])/(st_dev_*hatsigma)
    } 
  }
  vMaxAbsT = apply(mTbest,1,'max')
  imax = which.max(vMaxAbsT)
  M = mMbest[imax,]
  M = M[which(M>0.5)]
  list(M = M,i=which.max(mTbest[imax,]))
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