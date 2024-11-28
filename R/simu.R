#' Simulation of coalescent probabilities using either rejection sampling or multinomial sampling
#'
#' @param rt Reproduction number
#' @param vt Variance of offspring distribution (if equal to R_t then Poisson is used otherwise NegBin)
#' @param nt Population size at time t
#' @param nt1 Population size at time t+1
#' @param type Either inclusive or exclusive
#' @param n Number of lineages observed at time t
#' @param nrep Number of repeats to simulate
#' @param method Which method to use for simulation, options are: rejection or multinomial
#'
#' @return Vector of inclusive coalescent probabilities for number of lines from 1 to n
#' @export
#'
simul_pcoal=function(rt,vt=rt,nt,nt1=round(rt*nt),type='inclusive',n,nrep=1e3,method='rejection')
{
  if (type=='inclusive') n=nt1
  p=rep(0,n)
  for (i in 1:nrep) {
    if (method=='rejection') {
      ks=0
      while (sum(ks)!=nt1)
        if (vt==rt) ks=rpois(nt,rt)
      else ks=rnbinom(nt,rt*rt/(vt-rt),rt/vt)
    }
    if (method=='multinomial') {
      if (vt==rt) ks=rmultinom(1, nt1, rep(1/nt,nt))
      else ks=rmultinom(1, nt1, rgamma(nt,shape=rt*rt/(vt-rt),scale=1))#scale doesn't actually matter
    }
    if (type=='inclusive') p[1]=p[1]+1
    else {
      ks=tabulate(sample(rep(1:nt,ks),n),nbins=nt) #subset from Nt1 down to n
      p[1]=p[1]+sum(ks==1)/n
    }
    for (j in 2:n) {
      v=0:(j-1)
      p0=0
      for (k in ks) p0=p0+prod(k-v)*ifelse(type=='inclusive',k>=j,k==j)
      p[j]=p[j]+p0/prod(n-v)
    }
  }
  p=p/nrep
  return(p)
}

#' Comparison of simulated and calculated coalescent probabilities
#'
#' @param rt Reproduction number
#' @param vt Variance of offspring distribution (if equal to R_t then Poisson is used otherwise NegBin)
#' @param nt Population size at time t
#' @param nt1 Population size at time t+1
#' @param n Number of observed lineages
#' @param type Either 'inclusive' or 'exclusive'
#' @param nrep Number of repeats to simulate
#' @param method Which method to use for simulation, options are: rejection or multinomial
#'
#' @return Vector of inclusive coalescent probabilities for number of lines from 1 to nt1
#' @export
#'
compare_simu=function(rt=1.5,vt=10,nt=10,nt1=20,n=15,type='inclusive',nrep=1e3,method='multinomial')
{
  s=simul_pcoal(rt=rt,vt=vt,nt=nt,nt1=nt1,n=n,type=type,nrep=nrep,method=method)
  if (type=='inclusive') {
    if (rt==vt) p=pois_inclusive(1:nt1,nt)
    else {
      r=rt^2/(vt-rt)
      p=unlist(lapply(1:nt1,function(k) negbin_inclusive(k,nt,r)))
    }
  } else {
    if (rt==vt) p=pois_exclusive(1:n,n,nt)
    else {
      r=rt^2/(vt-rt)
      p=unlist(lapply(1:nt1,function(k) negbin_exclusive(k,n,nt,r)))
    }
  }

  ret=cbind(s,p,s-p)
  colnames(ret)=c('Simulated','Calculated','Difference')
  return(ret)
}

