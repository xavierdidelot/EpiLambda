#' Simulation of inclusive coalescent probabilities using either rejection sampling or multinomial sampling
#'
#' @param rt Reproduction number
#' @param vt Variance of offspring distribution (if equal to R_t then Poisson is used otherwise NegBin)
#' @param nt Population size at time t
#' @param nt1 Population size at time t+1
#' @param nrep Number of repeats to simulate
#' @param method Which method to use for simulation, options are: rejection or multinomial
#'
#' @return Vector of inclusive coalescent probabilities for number of lines from 1 to nt1
#' @export
#'
simul_inclusive=function(rt,vt=rt,nt,nt1=round(rt*nt),nrep=1e3,method='rejection')
{
  p=rep(0,nt1)
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
    p[1]=p[1]+1
    for (j in 2:nt1) {
      v=0:(j-1)
      p0=0
      for (k in ks) p0=p0+prod(k-v)*(k>=j)
      p[j]=p[j]+p0/prod(nt1-v)

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
#' @param nrep Number of repeats to simulate
#' @param method Which method to use for simulation, options are: rejection or multinomial
#'
#' @return Vector of inclusive coalescent probabilities for number of lines from 1 to nt1
#' @export
#'
compare_simu=function(rt=1.5,vt=10,nt=10,nt1=20,nrep=1e3,method='multinomial')
{
  s=simul_inclusive(rt=rt,vt=vt,nt=nt,nt1=nt1,nrep=nrep,method=method)
  if (rt==vt) p=pois_inclusive(1:nt1,nt)
  if (rt!=vt) {
    r=rt^2/(vt-rt)
    p=unlist(lapply(1:nt1,function(k) negbin_inclusive(k,nt,r)))
  }
  ret=cbind(s,p,s-p)
  colnames(ret)=c('Simulated','Calculated','Difference')
  return(ret)
}
