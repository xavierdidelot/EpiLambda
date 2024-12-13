#' Simulation from lambda-coalescent model
#'
#' @param n Number of lineages
#' @param lambda Lambda function of (n,k)
#'
#' @return Simulated phylogeny
#' @export
#'
lambda_simtree=function(n,lambda) {
  dates=rep(0,2*n)
  parents=rep(NA,2*n)
  nn=n+1
  tocoal=1:n
  curtime=0
  ll=0
  while (length(tocoal)>1) {
    l=length(tocoal)
    v=lambda(l,k=2:l)*choose(l,2:l)
    dt=rexp(1,rate=sum(v))
    ll=ll+dexp(dt,sum(v),log=T)
    if (l==2) siz=2 else {siz=sample(2:l,1,prob=v);ll=ll+log(v[siz-1]/sum(v))}
    if (siz<l) new=nn+1 else new=n+1#this is the root
    w=sample(tocoal,siz)
    ll=ll-lchoose(l,siz)
    curtime=curtime+dt
    dates[new]=curtime
    parents[w]=new
    tocoal=c(setdiff(tocoal,w),new)
    nn=nn+1
  }
  edge=cbind(parents[1:(nn-1)],1:(nn-1))
  edge=edge[-(n+1),]
  t=list()
  t$Nnode=nn-1-n
  t$edge=edge
  t$tip.label=1:n
  t$edge.length=dates[t$edge[,1]]-dates[t$edge[,2]]
  t$ll=ll
  class(t)<-'phylo'
  return(t)
}

#' Log-likelihood function for lambda-coalescent model
#'
#' @param t Tree
#' @param lambda Lambda function of (n,k)
#'
#' @return Log-likelihood
#' @export
#'
lambda_loglik=function(t,lambda) {
  ci=coalescent.intervals(multi2di(t))
  lins=ci$lineages
  dts=ci$interval.length
  rem=which(dts==0)
  if (length(rem)>0) {
    lins=lins[-rem]
    dts=dts[-rem]
  }
  ll=0
  for (i in 1:length(lins)) {
    l=lins[i]
    v=lambda(l,k=2:l)*choose(l,2:l)
    dt=dts[i]
    ll=ll+dexp(dt,sum(v),log=T)
    if (i==length(lins)) siz=l else siz=l-lins[i+1]+1
    if (l>2) {
      ll=ll+log(v[siz-1]/sum(v))
      ll=ll-lchoose(l,siz)
    }
  }
  return(ll)
}

#' Lambda function for Beta-coalescent
#'
#' @param n Number of lineages
#' @param k Number of lineages to coalesce
#' @param alpha Alpha parameter
#'
#' @return Rate of coalescence
#' @export
#'
beta_lambda=function(n,k,alpha) {
  if (alpha==2) return(1*(k==2))#Kingman coalescent
  if (alpha==0) return(1*(k==n))#Star phylogeny
  return(beta(k-alpha,n-k+alpha)/beta(2-alpha,alpha))
}

#' Size distribution for Beta-coalescent
#'
#' @param n Number of lineages
#' @param alpha Alpha parameter
#'
#' @return Distribution of size of the coalescence event
#' @export
#'
beta_psize=function(n,alpha) {
  v=beta_lambda(n,k=2:n,alpha)*choose(n,2:n)
  v=v/sum(v)
  return(v)
}

#' Simulation from beta-coalescent model
#'
#' @param n Number of lineages
#' @param alpha Alpha parameter
#'
#' @return Simulated phylogeny
#' @export
#'
beta_simtree=function(n,alpha) {
  lambda_simtree(n,function(n,k) beta_lambda(n,k,alpha=alpha))
}

#' Log-likelihood function for beta-coalescent model
#'
#' @param t Tree
#' @param alpha Alpha parameter
#'
#' @return Log-likelihood
#' @export
#'
beta_loglik=function(t,alpha) {
  lambda_loglik(t,function(n,k) beta_lambda(n,k,alpha=alpha))
}

#' Simulation from new lambda-coalescent model
#'
#' @param n Number of lineages
#' @param nt Population size parameter
#' @param r Dispersion parameter
#'
#' @return Simulated phylogeny
#' @export
#'
new_simtree=function(n,nt,r) {
  lambda=function(n,k) return(-log1p(-negbin_exclusive(k=k,n=n,nt=nt,r=r)))
  lambda_simtree(n,lambda)
}

#' MLE for beta-coalescent model
#'
#' @param t Tree
#'
#' @return Log-likelihood
#' @export
#'
beta_mle=function(t) {
  f=function(alpha) beta_loglik(t,alpha)
  o=optim(1,f,method='Brent',lower=0,upper=2,control=list(fnscale=-1))
  return(o$par)
}

#' Size distribution for new lambda-coalescent model
#'
#' @param n Number of lineages
#' @param nt Population size parameter
#' @param r Dispersion parameter
#'
#' @return Distribution of size of the coalescence event
#' @export
#
new_psize=function(n,nt,r) {
  lambda=function(n,k) return(-log1p(-negbin_exclusive(k=k,n=n,nt=nt,r=r)))
  v=lambda(n,k=2:n)*choose(n,2:n)
  v=v/sum(v)
  return(v)
}

#' Log-likelihood function for new lambda-coalescent model
#'
#' @param t Tree
#' @param nt Population size parameter
#' @param r Dispersion parameter
#'
#' @return Log-likelihood
#' @export
#'
new_loglik=function(t,nt,r) {
  lambda=function(n,k) return(-log1p(-negbin_exclusive(k=k,n=n,nt=nt,r=r)))
  lambda_loglik(t,lambda)
}

#' MLE for new lambda-coalescent model
#'
#' @param t Tree
#' @param nt Optional, if provided only the mle of r will be computed
#' @param r Optional, if provided only the mle of nt will be computed
#' @param ntbounds Optional, bounds for estimate of nt
#' @param rbounds Optional, bounds for estimate of r
#'
#' @return Maximum likelihood estimator of parameters of the new lambda-coalescent model
#' @export
#'
new_mle=function(t,nt,r,ntbounds,rbounds) {
  if (missing(ntbounds)) ntbounds=c(1,1000)
  if (missing(rbounds)) rbounds=c(0,10)
  if (missing(nt) && missing(r)) {
    #Select good starting values for both parameters
    nt=NA
    r=NA
    best=-Inf
    ss=seq(0.05,0.95,0.1)
    for (nts in ss*(ntbounds[2]-ntbounds[1])+ntbounds[1]) for (rs in ss*(rbounds[2]-rbounds[1])+rbounds[1]) {
      v=new_loglik(t,nts,rs)
      if (v>best) {best=v;nt=nts;r=rs}
    }
  }
  if (!missing(r) && !missing(nt)) {
    #Optimise both parameters
    f1=function(p) {
      if (p[1]<=ntbounds[1] || p[1]>=ntbounds[2] || p[2]<=rbounds[1] || p[2]>=rbounds[2]) return(-1e10)
      r=new_loglik(t,p[1],p[2])
      if (is.na(r)||is.infinite(r)) return(-1e10)
      return(r)
    }
    o=optim(c(nt,r),f1,method='L-BFGS-B',lower=c(ntbounds[1],rbounds[1]),upper=c(ntbounds[2],rbounds[2]),control=list(fnscale=-1))
    return(o$par)
  }
  if (missing(r)) {
    #Optimise only nt
    f2=function(r) {re=new_loglik(t,nt,r);if (is.na(re)||is.infinite(re)) return(-1e10) else return(re)}
    o=optim(1,f2,method='Brent',lower=rbounds[1],upper=rbounds[2],control=list(fnscale=-1))
    return(o$par)
  }
  if (missing(nt)) {
    #Optimise only r
    f3=function(nt) {re=new_loglik(t,nt,r);if (is.na(re)||is.infinite(re)) return(-1e10) else return(re)}
    o=optim(Ntip(t),f3,method='Brent',lower=ntbounds[1],upper=ntbounds[2],control=list(fnscale=-1))
    return(o$par)
  }
}
