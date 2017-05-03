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

# If we multiply each the first q columns of the polar cone matrix by an appropriate constant
# we get the (generalised) inverse of V[,1:lq] (the matrix that generates the cone
# -see Dattorro, Convex Optimization, page 151, expression 421),
# This is what following piece of (commented) code does:

# pcone<-polar.cone(V,2)$pc[,c(2,1)]
# t(-pcone)%*%V[,1:2] # non scaled; -pcone is the dual cone;
# #scaling the cols:
# col.scalars<-diag(diag(t(pcone)%*%V[,1:2])^-1)
# Xinv.tr<-pcone%*%col.scalars
# t(Xinv.tr)%*%V[,1:2]




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




