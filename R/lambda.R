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
  while (length(tocoal)>1) {
    l=length(tocoal)
    v=lambda(l,k=2:l)*choose(l,2:l)
    dt=rexp(1,rate=sum(v))
    if (l==2) siz=2 else siz=sample(2:l,1,prob=v)
    if (siz<l) new=nn+1 else new=n+1#this is the root
    w=sample(tocoal,siz)
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
  class(t)<-'phylo'
  return(t)
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
