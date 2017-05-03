
# Some R packages that need to be loaded:
# 'MASS'
# 'matrix'
# 'combinat'
# 'mvtnorm'
# 'numDeriv' version 2016.8-1 (I believe) from: http://cran.at.r-project.org/web/packages/numDeriv/index.html

# It requires the functions in 'boundary_functions.R'

# 'Posdef' generates a positive definite matrix, to be used 
# as variance-covariance matrix of the score statistic (denoted M_{0}
# in the thesis). 
# 'n' is the number of columns/rows, 'ev' is the vector of eigenvalues
# Source: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.html

M0<-Posdef(n=5,ev = runif(5, 0, 10)) 

# For a given variance-covariance matrix of the score statistic ('M0'),
# a given number of single parameters on the boundary ('lq'),
# and a given number of simulated observations of a multivariate standard
# normal random variable ('nz"), chisq.mix(M0,lq,nz) returns a vector with the
# mixing probabilities of a mixture of Chi-square distributions: the n-th
# component of the vector refers to a chi-square distribution with n-1 degrees
# of freedom. This mixture of distributions is the distribution of the log-
# likelihood ratio statistic when testing a hypothesis of the kind defined
# in section.... of the thesis.

mix<-chisq.mix(M0=M0,lq=4,nz=100000)


# 'cdf.mix' Returns the cdf of the mixture of chi-square distributions with
# mixing probabilities given by input vector 'mix':

cdf.mix(x=5.029245,mix=mix)


# 'rchisq.mix' returns 'n' simulated observations from a mixture of 
# chi-squared rv's when the mixing probabilities are given by 'mix':

vv<-rchisq.mix(n=1000000,mix=mix)

# 'rdev' returns 'n' simulated observations of the log-likelihood ratio
# statistic when the information matrix is given by M0 and
# there are 'lq' single parameters on the boundary:

dev<-rdev(n=100000,M0=M0,lq=4)

# to check whether 'vv' and 'dev' are likely to have been generated from
# the same distribution:

qqplot(vv,dev,ylim=c(0,20),xlim=c(0,20))
abline(a=0,b=1,col=2)

# to check whether a specific quantile of 'vv'
# is close to the corresponding quantile of 'dev':

quantile(dev,0.95)
quantile(vv,0.95)
cdf.mix(quantile(dev,0.95),mix=mix)

# Kolmogorov-Smirnov should not be used as there are many ties:
ks.test(dev,vv)

# creating better qq-plots of 'vv' vs 'dev':

def.par <- par(no.readonly = TRUE)

qqplot(quantile(vv,seq(0,1,0.01)),quantile(dev,seq(0,1,0.01)),ylim=c(0,12),xlim=c(0,12),ylab="expression (4.3)", xlab=expression("mixture of"~chi^2))
title(main="Q-Q Plot of simulated deviances",cex.main=1.2,line=2)
abline(a=0,b=1,col=2)

# Estimating M0 for the comparison ISO vs IM1:

load("C:/Users/ruibarrigana/Google Drive/PhD/RWORK/IIM paper code/Wang_Hey_sets.RData")
load("C:/Users/rui/Google Drive/PhD/RWORK/IIM paper code/Wang_Hey_sets.RData")

# Alternated loci for ingroup differences (x1a,x3bc):
x1<-x1a
x3<-x3bc

r1<-r.a
r3<-r.bc

fit.IM.1<-nlminb(rep(1.4,6),negll.IM.1,
                 lower=c(0,0,0.000001,0.000001,0.000001,0.000001))

M0<-jacobian(grad.negll, fit.IM.1$par, method="Richardson",
         side=c(1,NA,NA,NA,NA,NA))


M0<-nearPD(M0)$mat
M0<-matrix(M0,ncol=ncol(M0)) # changing the'type'class' of the matrix

M0<-M0[1:2,1:2]

save(M0.random_4_by_4, file = "M0.random_4_by_4.RData")

M0.jacobian<-M0
M0<-M0.jacobian
M0<-M0.random_4_by_4


#### SIMULATED DATA: 
### simulate 500 data sets from the iso.3 model (isolation model with symmetric
### population sizes);
### fit iso.3 model and symmetric IM model with bidirectional migration to each data set;
### calculate the 500 likelihood ratio test statistics;
### pick one data set and estimate the matrix M0 (hessian): reliable estimates are those
### based on data sets which have interior point MLE's (both migration rates 
### strictly greater than zero);
### estimate the mixture of chi squares based on the estimate of M0;
### simulate 500 observations from this mixture and build qq-plots to check accuracy.


source('C:/Users/rui/Google Drive/PhD/RWORK/IIM paper code/GIM_r_functions.R')

#simulating the data:

data1<-vector("list",500)
for(i in 1:500){
  data1[[i]]<-sim.GIM.4(a=1.25,theta=2,b=1,c1=1,c2=1,tau1=0,tau0=1.25,M1=0,M2=0,
                        M1c=0,M2c=0,N=c(10000,10000,20000))
  print(i)
}

r1<-rep(1,10000)
r2<-rep(1,10000)
r3<-rep(1,20000)

fit.data1.iso.3<-vector("list",500)
fit.data1.IM.5<-vector("list",500)

for (i in 1:500){
  x1<-data1[[i]]$x1
  x2<-data1[[i]]$x2
  x3<-data1[[i]]$x3
  
  fit.data1.IM.5[[i]]<-gim("IM.5",rep(0.5,5))
  print(i)
  print(fit.data1.IM.5[[i]]$message)    
  
}

for (i in 1:500){
  x1<-data1[[i]]$x1
  x2<-data1[[i]]$x2
  x3<-data1[[i]]$x3
  
  fit.data1.iso.3[[i]]<-gim("iso.3")
  print(i)
  print(fit.data1.iso.3[[i]]$message)    
  
}
save(list =c("data1","fit.data1.iso.3","fit.data1.IM.5"),
     file = "C:/Users/rui/Google Drive/PhD/RWORK/Thesis code/data1.RData")



#looking for cases of non-convergence:

ind.IM.5<-NULL
for (i in 1:500){
  if(fit.data1.IM.5[[i]]$message!="relative convergence (4)"){ind.IM.5<-c(ind.IM.5,i)}
  
}

length(ind.IM.5)

## after obtaining convergence for the above cases, calculate the 500 lrt's:

lrt.iso3_im5<-vector("numeric", 500)

for (i in 1:500){
  lrt.iso3_im5[i]<-2*(fit.data1.iso.3[[i]]$objective
                      -fit.data1.IM.5[[i]]$objective)
  
}

## estimating M0 and the chisquare mix based on the data set nr 7
## (the first data set with both estimates of M1 and M2 strictly positive)

x1<-data1[[7]]$x1
x2<-data1[[7]]$x2
x3<-data1[[7]]$x3

# Putting M1 and M2 as the first two parameters (required by the functions
# that calculate the chisquare mixture)

yy<-as.numeric(fit.data1.IM.5[[7]]$par)
mle<-c(yy[4:5],yy[1:3])
negll.IM.5<-function(params){
  
  a<-params[3]/params[4]
  theta<-params[4]
  b<-1
  tau<-params[5]/theta
  M<-vector(length=2,mode="numeric")
  M[1]<-params[1]
  M[2]<-params[2]
  
  
  
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4) 	#boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    loglike1<-log(colSums(A[,1]*((-nu*(theta*R1)^X1)/
                                   ((theta*R1-nu)^(X1+1))*
                                   (1-ppois(X1,tau*(-nu+theta*R1)))+
                                   exp(-tau*(theta*R1-nu))*
                                   (a*theta*R1)^X1/
                                   (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A[,2]*((-nu*(theta*R2)^X2)/
                                   ((theta*R2-nu)^(X2+1))*
                                   (1-ppois(X2,tau*(-nu+theta*R2)))+
                                   exp(-tau*(theta*R2-nu))*
                                   (a*theta*R2)^X2/
                                   (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A[,3]*((-nu*(theta*R3)^X3)/
                                   ((theta*R3-nu)^(X3+1))*
                                   (1-ppois(X3,tau*(-nu+theta*R3)))+
                                   exp(-tau*(theta*R3-nu))*
                                   (a*theta*R3)^X3/
                                   (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  if (M[1]==0 & M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  if (M[1]>0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
    nu11<--1/b
    A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
    nu12<--M[1]/2
    A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
    nu13<--(1+M[1])
    
    A1<-c(A11,A12,A13)
    nu1<-c(nu11,nu12,nu13)
    
    A2<-1
    nu2<--1/b
    
    A31<-b*M[1]/(b*M[1]-2)
    nu31<--1/b
    A32<-2/(2-b*M[1])
    nu32<--M[1]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
    
  }
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-log(((theta*R1[1,])^X1[1,])/				
                    ((1+theta*R1[1,])^(X1[1,]+1))*							
                    (1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
                    exp(-tau*(1+theta*R1[1,]))*
                    (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))
    
    loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
                    ((1+b*theta*R2[1,])^(X2[1,]+1))*							
                    (1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
                    exp(-tau*(1/b+theta*R2[1,]))*
                    (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))
    
    loglike3<-log(exp(-tau*theta*R3[1,])*			
                    (a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))
    
  }
  
  -sum(c(loglike1,loglike2,loglike3),na.rm=TRUE)
}

hess.IM.5<-hessian(negll.IM.5,mle)
#checking positive definiteness:
inv.hess.IM.5<-solve(hess.IM.5)
sqrt(diag(inv.hess.IM.5))

mix<-chisq.mix(M0=hess.IM.5/40000,lq=2,nz=100000)

vv<-rchisq.mix(n=1000,mix=mix)

qqplot(vv,lrt.iso3_im5,ylim=c(0,13),xlim=c(0,13),
       main="Q-Q Plot of simulated data",xlab=expression("estimated mixture of"~chi^2),
       ylab="likelihood ratio statistics")
abline(a=0,b=1,col=2)

qqplot(quantile(vv,seq(0,1,0.01)),quantile(lrt.iso3_im5,seq(0,1,0.01)),ylim=c(0,9),xlim=c(0,9),
       main="Q-Q Plot of simulated data",xlab=expression("percentiles of mixture of"~chi^2),
       ylab="LRT statistic distribution percentiles")
abline(a=0,b=1,col=2)

qqplot(qchisq(ppoints(100), df = 2),quantile(lrt.iso3_im5,seq(0,1,0.01)),ylim=c(0,13),xlim=c(0,13),
       main="Q-Q Plot",xlab=expression(chi[2]^2~~ "theoretical percentiles"),
       ylab="LRT statistic distribution percentiles")
abline(a=0,b=1,col=2)


## repeat the procedure to add another 500 data sets and increase accuracy

data11<-vector("list",500)
for(i in 1:500){
  data11[[i]]<-sim.GIM.4(a=1.25,theta=2,b=1,c1=1,c2=1,tau1=0,tau0=1.25,M1=0,M2=0,
                        M1c=0,M2c=0,N=c(10000,10000,20000))
  print(i)
}

r1<-rep(1,10000)
r2<-rep(1,10000)
r3<-rep(1,20000)

fit.data11.iso.3<-vector("list",500)
fit.data11.IM.5<-vector("list",500)

for (i in 1:500){
  x1<-data11[[i]]$x1
  x2<-data11[[i]]$x2
  x3<-data11[[i]]$x3
  
  fit.data11.IM.5[[i]]<-gim("IM.5",rep(0.5,5))
  print(i)
  print(fit.data11.IM.5[[i]]$message)    
  
}

for (i in 1:500){
  x1<-data11[[i]]$x1
  x2<-data11[[i]]$x2
  x3<-data11[[i]]$x3
  
  fit.data11.iso.3[[i]]<-gim("iso.3")
  print(i)
  print(fit.data11.iso.3[[i]]$message)    
  
}
save(list =c("data11","fit.data11.iso.3","fit.data11.IM.5"),
     file = "C:/Users/rui/Google Drive/PhD/RWORK/Thesis code/data11.RData")

#looking for cases of non-convergence:

ind.IM.5<-NULL
for (i in 1:500){
  if(fit.data11.IM.5[[i]]$message!="relative convergence (4)"){ind.IM.5<-c(ind.IM.5,i)}
  
}

length(ind.IM.5)

## after obtaining convergence for the above cases, calculate the 500 lrt's:


for (i in 1:500){
  lrt.iso3_im5<-c(lrt.iso3_im5,2*(fit.data11.iso.3[[i]]$objective
                      -fit.data11.IM.5[[i]]$objective))
  
}



## simulate instead from an asymmetric iso model and try to estimate
## the dbn of the likelihood ratio statistic using Godambe information

#simulating the data:

data2<-vector("list",500)
for(i in 1:500){
  data2[[i]]<-sim.GIM.4(a=1.25,theta=2,b=1.5,c1=1,c2=1,tau1=0,tau0=1,M1=0,M2=0,
                        M1c=0,M2c=0,N=c(10000,10000,20000))
  print(i)
}

r1<-rep(1,10000)
r2<-rep(1,10000)
r3<-rep(1,20000)

fit.data2.iso.1<-vector("list",500)
fit.data2.IM.1<-vector("list",500)

for (i in 161:500){
  x1<-data2[[i]]$x1
  x2<-data2[[i]]$x2
  x3<-data2[[i]]$x3
  
  fit.data2.IM.1[[i]]<-gim("IM.1",rep(1.2,6))
  print(i)
  print(fit.data2.IM.1[[i]]$message)    
  
}

for (i in 1:500){
  x1<-data2[[i]]$x1
  x2<-data2[[i]]$x2
  x3<-data2[[i]]$x3
  
  fit.data2.iso.1[[i]]<-gim("iso.1")
  print(i)
  print(fit.data2.iso.1[[i]]$message)    
  
}

#looking for cases of non-convergence:

ind.IM.5<-NULL
for (i in 1:500){
  if(fit.data2.IM.5[[i]]$message=="relative convergence (4)"){
    ind.IM.5<-c(ind.IM.5,i)
    print(i)
    }
  
}

length(ind.IM.5)

## after obtaining convergence for the above cases, calculate the 500 lrt's:

lrt.data2.iso3_im5<-vector("numeric", 500)

for (i in 1:500){
  lrt.data2.iso3_im5[i]<-2*(fit.data2.iso.3[[i]]$objective
                      -fit.data2.IM.5[[i]]$objective)
  
}

########

## estimating M0 for the IM1 model:

#putting the migration parameter as the first parameter
negll.IM.1<-function(params){
  
  M<-vector(length=2,mode="numeric")
  M[1]<-params[1]
  M[2]<-params[2]
  a<-params[3]/params[4]
  theta<-params[4]
  b<-params[5]
  tau<-params[6]/theta
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4)   #boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    loglike1<-log(colSums(A[,1]*((-nu*(theta*R1)^X1)/
                                   ((theta*R1-nu)^(X1+1))*
                                   (1-ppois(X1,tau*(-nu+theta*R1)))+
                                   exp(-tau*(theta*R1-nu))*
                                   (a*theta*R1)^X1/
                                   (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A[,2]*((-nu*(theta*R2)^X2)/
                                   ((theta*R2-nu)^(X2+1))*
                                   (1-ppois(X2,tau*(-nu+theta*R2)))+
                                   exp(-tau*(theta*R2-nu))*
                                   (a*theta*R2)^X2/
                                   (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A[,3]*((-nu*(theta*R3)^X3)/
                                   ((theta*R3-nu)^(X3+1))*
                                   (1-ppois(X3,tau*(-nu+theta*R3)))+
                                   exp(-tau*(theta*R3-nu))*
                                   (a*theta*R3)^X3/
                                   (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  if (M[1]==0 & M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  if (M[1]>0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
    nu11<--1/b
    A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
    nu12<--M[1]/2
    A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
    nu13<--(1+M[1])
    
    A1<-c(A11,A12,A13)
    nu1<-c(nu11,nu12,nu13)
    
    A2<-1
    nu2<--1/b
    
    A31<-b*M[1]/(b*M[1]-2)
    nu31<--1/b
    A32<-2/(2-b*M[1])
    nu32<--M[1]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
    
  }
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-log(((theta*R1[1,])^X1[1,])/  			
                    ((1+theta*R1[1,])^(X1[1,]+1))*							
                    (1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
                    exp(-tau*(1+theta*R1[1,]))*
                    (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))
    
    loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
                    ((1+b*theta*R2[1,])^(X2[1,]+1))*							
                    (1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
                    exp(-tau*(1/b+theta*R2[1,]))*
                    (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))
    
    loglike3<-log(exp(-tau*theta*R3[1,])*			
                    (a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))
    
  }
  
  -sum(c(loglike1,loglike2,loglike3),na.rm=TRUE)
}
lb.IM.1<-c(0,0,0.000001,0.000001,0.000001,0.000001)

x1<-data2[[1]]$x1
x2<-data2[[1]]$x2
x3<-data2[[1]]$x3
  
gim("IM.6",rep(1.2,4))
fit.data2.IM.5[[1]]$par
  
mle<-as.numeric(gim("IM.6",rep(1.2,4))$par)

hess.IM.6<-hessian(negll.IM.6,mle)
#checking positive definiteness:
inv.hess.IM.6<-solve(hess.IM.6)
sqrt(diag(inv.hess.IM.6))

#estimating the Godambe information:
negll.IM.6.x1<-function(params){
  
  M<-vector(length=2,mode="numeric")
  M[2]<-params[1]
  a<-params[2]/params[3]
  theta<-params[3]
  b<-1
  tau<-params[4]/theta


  
  
  if (M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)

    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
      }
  
  
  if (M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    
    
    loglike1<-log(((theta*R1[1,])^X1[1,])/				
                    ((1+theta*R1[1,])^(X1[1,]+1))*							
                    (1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
                    exp(-tau*(1+theta*R1[1,]))*
                    (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))
      }
  
  -sum(loglike1,na.rm=TRUE)
}
negll.IM.6.x2<-function(params){
  
  M<-vector(length=2,mode="numeric")
  M[2]<-params[1]
  a<-params[2]/params[3]
  theta<-params[3]
  b<-1
  tau<-params[4]/theta

  
  
  if (M[2]>0){
    
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    
    
    
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
     }
  
  
  if (M[2]==0){
    
    
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    
    
    
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    
    
    loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
                    ((1+b*theta*R2[1,])^(X2[1,]+1))*							
                    (1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
                    exp(-tau*(1/b+theta*R2[1,]))*
                    (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))
    
    }
  
  -sum(loglike2,na.rm=TRUE)
}
negll.IM.6.x3<-function(params){
  
  M<-vector(length=2,mode="numeric")
  M[2]<-params[1]
  a<-params[2]/params[3]
  theta<-params[3]
  b<-1
  tau<-params[4]/theta
  
  
  
  if (M[2]>0){
    
   R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
   X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  
  if (M[2]==0){
    
   R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
   X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike3<-log(exp(-tau*theta*R3[1,])*			
                    (a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))
    
  }
  
  -sum(loglike3,na.rm=TRUE)
}

r1.all<-r1
r2.all<-r2
r3.all<-r3

x1.all<-x1
x2.all<-x2
x3.all<-x3

grad1<-NULL
for (i in 1:length(x1.all)){   
  r1<-r1.all[i]
  x1<-x1.all[i]
  grad1<-cbind(grad1,grad(negll.IM.6.x1,mle))
  print(i)
}

outer1<-matrix(rep(0,length(mle)^2),nrow=length(mle))

for (i in 1:length(x1.all)){
  outer1<-outer1+outer(grad1[,i],grad1[,i])    
}

grad2<-NULL
for (i in 1:length(x2.all)){
  r2<-r2.all[i]
  x2<-x2.all[i]
  grad2<-cbind(grad2,grad(negll.IM.6.x2,mle))
}

outer2<-matrix(rep(0,length(mle)^2),nrow=length(mle))

for (i in 1:length(x2.all)){
  outer2<-outer2+outer(grad2[,i],grad2[,i])    
}

grad3<-NULL
for (i in 1:length(x3.all)){
  r3<-r3.all[i]
  x3<-x3.all[i]
  grad3<-cbind(grad3,grad(negll.IM.6.x3,mle))
  print(i)
}

outer3<-matrix(rep(0,length(mle)^2),nrow=length(mle))

for (i in 1:length(x3.all)){
  outer3<-outer3+outer(grad3[,i],grad3[,i])    
}

sig.IM.6<-outer1+outer2+outer3

dev<-rdev(n=500,M0=hess.IM.6/40000,sig=sig.IM.6/40000,lq=1)

hist(dev)
sort(dev)
hist(lrt.data2.iso3_im5)

####### qqplot of section 4.1



#percentiles of the distribution of the likelihood ratio statistic under 
#misspecification, using the approximation in Jesus & Chandler, equation 3.6, 
#compared with the percentiles of the chisquare dbn with 2 degrees of freedom

# constants a, b, c estimated using code in "paper_revision_code.R"

a<-1.39
b<-1.81
c<-0.08

par(mar=c(5.1, 4.5, 4.1, 2.1))
qqplot(a*qchisq(seq(0.01,1,0.01),b)+c,qchisq(seq(0.01,1,0.01),2),
       ylim=c(0,12),xlim=c(0,12),
       xlab =expression(paste("distribution of 1.39 ", italic("X"))~+0.08~paste(",",italic(" X ~"))~chi[1.81]^2),
       ylab=expression(chi[2]^2~~"distribution"))
title(main="Q-Q Plot",cex.main=1.2,line=2)
abline(0,1,col=2)

