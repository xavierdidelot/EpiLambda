#' Simulation of p_nt using either rejection sampling or multinomial sampling
#'
#' @param rt Reproduction number
#' @param nt Population size at time t
#' @param vt Variance of offspring distribution (if equal to R_t then Poisson is used otherwise NegBin)
#' @param nt1 Population size at time t+1
#' @param nrep Number of repeats to simulate
#' @param max Maximum value of n to consider in p_nt
#' @param method Which method to use for simulation, options are: rejection or multinomial
#'
#' @return Vector of p_nt values for n from 1 to max
#' @export
#'
simul_pnt=function(rt,nt,vt=rt,nt1=round(rt*nt),nrep=1e3,max=10,method='rejection')
{
  p=rep(0,max)
  for (i in 1:nrep) {
    if (method=='rejection') {
      ks=0
      while (sum(ks)!=nt1)
        if (vt==rt) ks=rpois(nt,rt)
          else ks=rnbinom(nt,rt*rt/(vt-rt),rt/vt)
    }
    if (method=='multinomial') {
      if (vt==rt) ks=rmultinom(1, nt1, rep(1/nt,nt))
        else ks=rmultinom(1, nt1, rgamma(nt,shape=rt*rt/(vt-rt),scale=1))#scale doesn't matter
    }
    p[1]=p[1]+1
    for (j in 2:max) {
      v=0:(j-1)
      p0=0
      for (k in ks) p0=p0+prod(k-v)*(k>=j)
      p[j]=p[j]+p0/prod(nt1-v)

    }
  }
  p=p/nrep
  return(p)
}
