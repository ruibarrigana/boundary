
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

mix<-chisq.mix(M0=M0,lq=3,nz=100000)


# 'cdf.mix' Returns the cdf of the mixture of chi-square distributions with
# mixing probabilities given by input vector 'mix':

cdf.mix(x=5.029245,mix=mix)


# 'rchisq.mix' returns 'n' simulated observations from a mixture of 
# chi-squared rv's when the mixing probabilities are given by 'mix':

vv<-rchisq.mix(n=1000000,mix=mix)

# 'rdev' returns 'n' simulated observations of the log-likelihood ratio
# statistic when the information matrix is given by M0 and
# there are 'lq' single parameters on the boundary:

dev<-rdev(n=100000,M0=M0,lq=3)

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


