source("functions_logistic.r")

#test that the pseudo target is the true coefficient for well-specified models
n=25
p=5
beta=c(3,0,-2,0,0)
M = c(1,3,5)
j=3
X = matrix(nrow=n,ncol=p,data=runif(n*p))
theta = h_logistic(X%*%beta)
target = target_logistic(X,M,theta,j)

#test of consistency of restricted maximum likelihood in well-specified submodel
n=10000
p=5
beta=c(3,0,-2,0,0)
M = c(1,3,5)
j=3
maxL1=10
X = matrix(nrow=n,ncol=p,data=runif(n*p))
theta = h_logistic(X%*%beta)
Y = rbinom(n=n,size=1,prob=theta)
target = target_logistic(X,M,theta,j)
hattarget = hat_beta_logistic(X,M,Y,maxL1)[which(M==j)]

#test of consistency of restricted maximum likelihood in misspecified submodel
n=10000
p=5
beta=c(3,1,-2,-3,1)
M = c(1,3,4)
j=3
maxL1=10
X = matrix(nrow=n,ncol=p,data=runif(n*p))
theta = h_logistic(X%*%beta)
Y = rbinom(n=n,size=1,prob=theta)
target = target_logistic(X,M,theta,j)
hattarget = hat_beta_logistic(X,M,Y,maxL1)[which(M==j)]


#Monte Carlo test of the asymptotic standard deviation for well-specified model
nmc=100
n=1000
p=5
beta=c(3,0,-2,0,0)
M = c(1,3,5)
j=5
maxL1=20
X = matrix(nrow=n,ncol=p,data=runif(n*p))
theta = h_logistic(X%*%beta)
target = target_logistic(X,M,theta,j)
vError = 0*(1:nmc)
vsd = 0*(1:nmc)
for (imc in 1:nmc) {
  Y = rbinom(n=n,size=1,prob=theta)
  hattarget = hat_beta_logistic(X,M,Y,maxL1)[which(M==j)]
  vError[imc] = target - hattarget
  vsd[imc] = st_dev_logistic(X,M,Y,j,maxL1)
}
mean(vsd)
sd(vsd)
sd(vError)
hist(vError,breaks=10)

#Test that Balpha > K1 
Kone = function(X,alpha,I=1000) {
  n = dim(X)[1]
  vC = maxOfInnerProducs(X,I) #I realizations of maximum of inner product
  quantile(vC,alpha)
}
n=15
p=5
X = matrix(nrow=n,ncol=p,data=runif(n*p))
I=1000
alpha=0.9
K1 = Kone(X,alpha,I)
Ba = Balpha(p,p*2^{p-1},alpha,I=I) 

#Test that Balpha = 1.64 when alpha=0.9 and p=1 
p=1
I=1000
alpha=0.9
Ba = Balpha(p,p*2^{p-1},alpha,I=I) 

#Test that Balpha == quantile of sqrt(X2) when number of unit vectors large
p=4
I=1000
alpha=0.9
N=100000
Ba = Balpha(p,N,alpha,I=I) 
qX2 = sqrt(qchisq(p=alpha,df=p))

#Test of the asymptotic for Ba
p=1520
N=1.5^p
I=10000
alpha=0.9
Ba = Balpha(p,N,alpha,I=I) 
Baass = sqrt(p*(1-N^(-2/(p-1))))


#Monte Carlo test of the naive asymptotic standard deviation for well-specified model
nmc=100
n=1000
p=5
beta=c(3,0,-2,0,0)
M = c(1,3,5)
j=5
maxL1=20
X = matrix(nrow=n,ncol=p,data=runif(n*p))
theta = h_logistic(X%*%beta)
target = target_logistic(X,M,theta,j)
vError = 0*(1:nmc)
vsd = 0*(1:nmc)
for (imc in 1:nmc) {
  Y = rbinom(n=n,size=1,prob=theta)
  hattarget = hat_beta_logistic(X,M,Y,maxL1)[which(M==j)]
  vError[imc] = target - hattarget
  vsd[imc] = st_dev_logistic_naive(X,M,Y,j,maxL1)
}
mean(vsd)
sd(vsd)
sd(vError)
hist(vError,breaks=10)

#Monte Carlo test that naive confidence intervals have
#asymptotically nominal coverage for well-specified model
nmc=1000
alpha=0.9
n=1000
p=5
beta=c(-1,5,0,0,2)
M = c(1,2,4,5)
j=2
maxL1=20
X = matrix(nrow=n,ncol=p,data=runif(n*p))
theta = h_logistic(X%*%beta)
target = target_logistic(X,M,theta,j)
vCoverage = 0*(1:nmc)
for (imc in 1:nmc) {
  Y = rbinom(n=n,size=1,prob=theta)
  IC = IC_logistic_naive(X,M,Y,j,alpha,maxL1)
  vCoverage[imc] = (IC$lower <= target) && (target <= IC$upper)
}
mean(vCoverage)

#Test that POSI intervals (modified by plug-in estimate of asymptotic covariance)
#with significance hunting with nbest = card all models
#and with beta_0 = 0
#has coverage equal to the nominal level
p=3
n=100
I=10000
alpha=0.6
X = matrix(data=runif(n*p),nrow=n,ncol=p)
beta_0 = c(0,0,0)
theta = h_logistic(X%*%beta_0)
nmc=1000
vLength = (-1)*(1:nmc)
vCoverage = (-1)*(1:nmc)
lambda=4
nbest=2^p-1
periodMsg = 100
maxL1=20
for (i in 1:nmc) {
  if (floor(i/periodMsg) == i/periodMsg) {
    cat("i = ",i,"\n")
  }
  Y = rbinom(n=n,size=1,prob=theta)
  sig_hun = significance_hunting_logistic(X,Y,lambda,nbest,maxL1)
  M = sig_hun$M
  itarget = sig_hun$i
  target = target_logistic(X,M,theta,itarget)
  IC = IC_plug_in_logistic(X,M,Y,itarget,alpha,I)
  vLength[i]=IC$upper - IC$lower
  vCoverage[i] = (IC$lower <= target) && (target <= IC$upper)
}
mean(vCoverage)


#Test that POSI intervals with significance hunting with nbest = card all models
#and with beta_0 different from 0
#has coverage larger than the nominal level
p=3
n=100
I=10000
alpha=0.9
X = matrix(data=runif(n*p),nrow=n,ncol=p)
beta_0 = c(1,0,0)
theta = h_logistic(X%*%beta_0)
nmc=1000
vLength = (-1)*(1:nmc)
vCoverage = (-1)*(1:nmc)
lambda=4
nbest=2^p-1
periodMsg = 100
maxL1=20
for (i in 1:nmc) {
  if (floor(i/periodMsg) == i/periodMsg) {
    cat("i = ",i,"\n")
  }
  Y = rbinom(n=n,size=1,prob=theta)
  sig_hun = significance_hunting_logistic(X,Y,lambda,nbest,maxL1)
  M = sig_hun$M
  itarget = sig_hun$i
  target = target_logistic(X,M,theta,itarget)
  IC = IC_plug_in_logistic(X,M,Y,itarget,alpha,I)
  vLength[i]=IC$upper - IC$lower
  vCoverage[i] = (IC$lower <= target) && (target <= IC$upper)
}
mean(vCoverage)