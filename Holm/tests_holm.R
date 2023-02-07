source("functions_holm.r")

#Test of alltbmj. We need to find 2^(k-1) times each k dimensional base vector
X=diag(3)
Mtbmj = alltbMj(X)

#Test of K1 as quantile of max of $p$ gaussian variables for orthogonal designs
p=10
n=20
I=10000
alpha=0.9
X = rbind(diag(p),matrix(data=0,nrow=n-p,ncol=p))
K1 = Kone(X,alpha,I=I)
sample_max = apply(matrix(nrow=I,ncol=p,data=abs(rnorm(n=I*p))),1,max)
K1theo = quantile(x=sample_max,probs=alpha)

#test of coverage for random selected model and Gaussian data
n=50
p=10
X = matrix(nrow=n,ncol=p,data=runif(n*p))
alpha=0.9
nmc=500
theta=runif(n)
sigma=1
I=1000
vLength = (-1)*(1:nmc)
vCoverage = (-1)*(1:nmc)
periodMsg = 100
for (i in 1:nmc) {
  if (floor(i/periodMsg) == i/periodMsg) {
    cat("i = ",i,"\n")
  }
  Y = theta + 4*rnorm(n)
  M = sample(1:p,size=1,replace=FALSE)
  itarget=M[1]
  target = target_holm(X,M,theta,itarget)
  IC = IC_holm(X,M,Y,itarget,alpha,I)
  vLength[i]=IC$upper - IC$lower
  vCoverage[i] = (IC$lower <= target) && (target <= IC$upper)
}

#Test that POSI intervals with significance hunting with nbest = card all models
#and with beta_0 = 0
#has coverage equal to the nominal level
p=3
n=100
I=10000
alpha=0.9
X = matrix(data=runif(n*p),nrow=n,ncol=p)
beta_0 = c(0,0,0)
theta = X%*%beta_0
nmc=1000
vLength = (-1)*(1:nmc)
vCoverage = (-1)*(1:nmc)
sigma=1
lambda=4
nbest=2^p-1
periodMsg = 100
for (i in 1:nmc) {
  if (floor(i/periodMsg) == i/periodMsg) {
    cat("i = ",i,"\n")
  }
  Y = theta + sigma*rnorm(n)
  sig_hun = significance_hunting(X,Y,lambda,nbest)
  M = sig_hun$M
  itarget = sig_hun$i
  target = target_holm(X,M,theta,itarget)
  IC = IC_holm(X,M,Y,itarget,alpha,I)
  vLength[i]=IC$upper - IC$lower
  vCoverage[i] = (IC$lower <= target) && (target <= IC$upper)
}

#Test that POSI intervals with significance hunting with nbest = card all models
#and with beta_0 different from 0
#has coverage larger than the nominal level
p=3
n=100
I=10000
alpha=0.9
X = matrix(data=runif(n*p),nrow=n,ncol=p)
beta_0 = c(2,0,-1)
theta = X%*%beta_0
nmc=1000
vLength = (-1)*(1:nmc)
vCoverage = (-1)*(1:nmc)
sigma=1
lambda=4
nbest=2^p-1
periodMsg = 100
for (i in 1:nmc) {
  if (floor(i/periodMsg) == i/periodMsg) {
    cat("i = ",i,"\n")
  }
  Y = theta + sigma*rnorm(n)
  sig_hun = significance_hunting(X,Y,lambda,nbest)
  M = sig_hun$M
  itarget = sig_hun$i
  target = target_holm(X,M,theta,itarget)
  IC = IC_holm(X,M,Y,itarget,alpha,I)
  vLength[i]=IC$upper - IC$lower
  vCoverage[i] = (IC$lower <= target) && (target <= IC$upper)
}

#Test that naive CI  = posi CI when p=1
n=50
p=1
X = matrix(data=runif(n*p),nrow=n,ncol=p)
M = 1
j=1
alpha=0.95
I=10000
Y = runif(n=n)
IC1 = IC_holm(X,M,Y,j,alpha,I)
IC2 = IC_naive_holm(X,M,Y,j,alpha)


#Test that naive CI length <= posi CI length 
n=50
p=3
X = matrix(data=runif(n*p),nrow=n,ncol=p)
M = 1
j=1
alpha=0.95
I=10000
Y = runif(n=n)
IC1 = IC_holm(X,M,Y,j,alpha,I)
IC2 = IC_naive_holm(X,M,Y,j,alpha)
length1 =IC1$upper - IC1$lower
length2 =IC2$upper - IC2$lower




#Test that naive intervals with fixed model (containing true mean) and component
#has coverage equal to the nominal level
p=3
n=1000
alpha=0.9
X = matrix(data=runif(n*p),nrow=n,ncol=p)
beta_0 = c(1,0,3)
M = c(1,3)
itarget=3
theta = X%*%beta_0
nmc=1000
vLength = (-1)*(1:nmc)
vCoverage = (-1)*(1:nmc)
sigma=1
periodMsg = 100
for (i in 1:nmc) {
  if (floor(i/periodMsg) == i/periodMsg) {
    cat("i = ",i,"\n")
  }
  Y = theta + sigma*rnorm(n)
  target = target_holm(X,M,theta,itarget)
  IC = IC_naive_holm(X,M,Y,itarget,alpha)
  vCoverage[i] = (IC$lower <= target) && (target <= IC$upper)
}

#Test if naive intervals with significance hunting
#has coverage lower than the nominal level 
p=3
n=100
alpha=0.9
X = matrix(data=runif(n*p),nrow=n,ncol=p)
beta_0 = c(0,0,0)
theta = X%*%beta_0
nmc=1000
vLength = (-1)*(1:nmc)
vCoverage = (-1)*(1:nmc)
sigma=1
lambda=4
nbest=2^p-1
periodMsg = 100
for (i in 1:nmc) {
  if (floor(i/periodMsg) == i/periodMsg) {
    cat("i = ",i,"\n")
  }
  Y = theta + sigma*rnorm(n)
  sig_hun = significance_hunting(X,Y,lambda,nbest)
  M = sig_hun$M
  itarget = sig_hun$i
  target = target_holm(X,M,theta,itarget)
  IC = IC_naive_holm(X,M,Y,itarget,alpha)
  vCoverage[i] = (IC$lower <= target) && (target <= IC$upper)
}




# Test that for n<p, the length of the CI with PoSI upper bound based on 1-sparse submodels are approximately equal to
# those equal to 2*st_dev_holm*hat_sigma*sqrt(n*(1-p^(-2/(n-1))))
p=10000
n=500
X = matrix(data=runif(n*p),nrow=n,ncol=p)
Y = runif(n)
M=6
itarget=6
alpha=0.9
s=1
I=1000
IC = IC_s_holm(X,M,Y,itarget,alpha,I,s)
length1 = IC$upper - IC$lower
length2 = 2*hat_sigma_holm(X,M,Y)*st_dev_holm(X,M,itarget)*sqrt(n*(1-p^(-2/(n-1))))


# Test that for n<p, the length of the CI with PoSI upper bound based on 2-sparse submodels are approximately equal to
# those equal to 2*st_dev_holm*hat_sigma*sqrt(n*(1-|M_s|^(-2/(n-1))))
# with |M_s| = p+2*p*(p-1)
p=10000
n=500
X = matrix(data=runif(n*p),nrow=n,ncol=p)
Y = runif(n)
M=6
itarget=6
alpha=0.9
s=2
I=1000
IC = IC_s_holm(X,M,Y,itarget,alpha,I,s)
length1 = IC$upper - IC$lower
N = p +2*p*(p-1)
length2 = 2*hat_sigma_holm(X,M,Y)*st_dev_holm(X,M,itarget)*sqrt(n*(1-N^(-2/(n-1))))





