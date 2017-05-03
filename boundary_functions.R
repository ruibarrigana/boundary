# A function to generate a positive definite matrix, to be used 
# as variance-covariance matrix of the score statistic (denoted M_{0}
# in the thesis).
# Source: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.html

Posdef<-function(n,ev=runif(n,0,10)) {
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

# The function 'region' below will make use of this function:
polar.cone<-function(V,lq){
  q<-1:lq
  combis<-combn(q,lq-1)
  combis2<-rbind(combis,combis[1,])
  U<-matrix(rep(NA,nrow(V)*lq),ncol=lq)
  for (i in 1:lq){
    U[,i]<-Null(cbind(V[,combis2[,i]],V[,-q]))  
    dum<-solve(V)%*%U[,i]
    if(sum(dum[q]>0)==lq){ #check whether U[,i] is inside the cone (it shouldn't be)
      U[,i]<--U[,i]
    }
    if(sum(U[,i]*V[,q[-combis[,i]]])>0){ #check that the angle between U[,i] and each of the other edges of the cone is higher than 90 degrees  
      U[,i]<--U[,i]
    }
  }
  list(pc=cbind(U,V[,-q]),combis=combis)
}

# if we multiply each the first q columns of the polar cone matrix by an appropriate constant
# we get the (generalised) inverse of V[,1:lq] (the matrix that generates the cone), 
# - see Dattorro, Convex Optimization, page 151, expression 421

pcone<-polar.cone(V,2)$pc[,c(2,1)]
t(-pcone)%*%V[,1:2] # non scaled; -pcone is the dual cone;
#scaling the cols:
col.scalars<-diag(diag(t(pcone)%*%V[,1:2])^-1)
Xinv.tr<-pcone%*%col.scalars
t(Xinv.tr)%*%V[,1:2]



# The function 'region' which determines which region R_k
# the vector z falls into, to be used by the function 'chisq.mix' below:

region<-function(z,V,U,lq,combis){
  
  q<-1:lq
  dum<-solve(U)%*%z
    if(sum(dum[q]>0)==lq){
      return(0)
    }
  
  dum<-solve(V)%*%z
  if(sum(dum[q]>0)==lq){
    return(lq)
  }
  
  for(k in 1:(lq-1)){
    combis.k.pol<-combn(q,k)
    
    for(i in 1:ncol(combis.k.pol)){
      lst<-vector("list",k)
      for(j in 1:k){
        lst[[j]]<-combis[,combis.k.pol[j,i]] 
      }
      
      UV<-cbind(U[,combis.k.pol[,i]],V[,Reduce(intersect,lst)],V[,-q])
      dum<-solve(UV)%*%z
      if(sum(dum[q]>0)==lq){
        return(lq-k)
      }
    }
    
  }
  
}


# A function that returns the mixing probabilities for the chi-square
# distribution:

chisq.mix<-function(M0,lq,nz){
  
  sqL<-diag(sqrt(eigen(M0)$values))
  Pt<-t(eigen(M0)$vectors)
  V<-sqL%*%Pt 
  # the above matrix V corresponds to the square root of M0 in the thesis
  
  U<-polar.cone(V,lq)$pc
  combis<-polar.cone(V,lq)$combis
  
  z<-mvrnorm(nz,rep(0,ncol(M0)),diag(rep(1,ncol(M0))))
  which.region<-apply(z,1,region,lq=lq, V=V, U=U,combis=combis)
  
  out<-rep(NA,lq+1)
  
  for(i in 0:(lq)){
    out[i+1]<-sum(which.region==i)/nz
  }
  out
}


# 'cdf.mix' Returns the cdf of the mixture of chi-square distributions with
# mixing probabilities given by input vector 'mix':

cdf.mix<-function(x,mix){
  cdfs<-rep(NA,length(mix))
  for(i in 1:length(mix)){
    cdfs[i]<-pchisq(x,i-1)
  }
  sum(mix*cdfs)
} 

# 'sq.dist' returns the squared distance of x to its projection on the parameter
# space of x, when
# x is multivariate normal with var-cov matrix equal to the inverse of M0: 

sq.dist<-function(y,x,M0){
  d<-x-y
  t(d)%*%M0%*%d
}


# 'rchisq.mix' returns 'n' simulated observations from a mixture of 
# chi-squared rv's when the mixing probabilities are given by 'mix':

rchisq.mix<-function(n,mix){
  vv<-rep(NA,n) 
  for(i in 1:n){
    ru<-runif(1)
    for(j in 1:length(mix)){
    if(ru<sum(mix[1:j])) break
    }
    vv[i]<-rchisq(1,j-1)
  } 
  vv
}

# 'rdev' returns 'n' simulated observations of the log-likelihood ratio
# statistic when the information matrix is given by M0 and
# there are 'lq' single parameters on the boundary:

# rdev<-function(n,M0.null,sig_hat.null=M0.null,M0.alt,sig_hat.null=M0.alt,lq){
#   invM0<-solve(M0)
#   p<-ncol(M0)
#   dev<-rep(NA,n)
#   for(i in 1:n){
#     x<-mvrnorm(1,rep(0,p),invM0%*%sig%*%t(invM0))
#     alt.min<-nlminb(rep(2,p),sq.dist,x=x,M0=M0,lower=c(rep(0,lq),rep(-Inf,p-lq)),
#                     upper=Inf)$objective
#     null.min<-nlminb(rep(2,p),sq.dist,x=x,M0=M0,lower=c(rep(0,lq),rep(-Inf,p-lq)),
#                      upper=c(rep(0,lq),rep(Inf,p-lq)))$objective
#     dev[i]<-null.min-alt.min
#   }
#   dev
# }

rdev<-function(n,M0,lq){
  invM0<-solve(M0)
  p<-ncol(M0)
  dev<-rep(NA,n)
  for(i in 1:n){
    x<-mvrnorm(1,rep(0,p),invM0)
    alt.min<-nlminb(rep(2,p),sq.dist,x=x,M0=M0,lower=c(rep(0,lq),rep(-Inf,p-lq)),
                    upper=Inf)$objective
    null.min<-nlminb(rep(2,p),sq.dist,x=x,M0=M0,lower=c(rep(0,lq),rep(-Inf,p-lq)),
                     upper=c(rep(0,lq),rep(Inf,p-lq)))$objective
    dev[i]<-null.min-alt.min
  }
  dev
}

# The next functions are needed to estimate M0 for the ISO vs IM1 comparison:

negll.IM.1<-function(params){
  
  M<-vector(length=2,mode="numeric")
  M[1]<-params[1]
  M[2]<-params[2]
  a<-params[3]/params[4]
  theta<-params[4]
  b<-params[5]/theta
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
  
  -sum(c(loglike1,loglike2,loglike3),na.rm=TRUE)/length(c(loglike1,loglike2,loglike3))
}


grad.negll<-function(x){
  grad(negll.IM.1,x, side=c(1,NA,NA,NA,NA,NA))
}




polar.cone<-function(V,lq){
  q<-1:lq
  combis<-combn(q,lq-1)
  combis2<-rbind(combis,combis[1,])
  U<-matrix(rep(NA,nrow(V)*lq),ncol=lq)
  for (i in 1:lq){
    U[,i]<-Null(cbind(V[,combis2[,i]],V[,-q]))  
    dum<-solve(V)%*%U[,i]
    if(sum(dum[q]>0)==lq){
      U[,i]<--U[,i]
    }
    if(sum(U[,i]*V[,q[-combis[,i]]])>0){
      U[,i]<--U[,i]
    }
  }
  list(pc=cbind(U,V[,-q]),combis=combis)
}