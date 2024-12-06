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
simul_pcoal_mixed_pois=function(lambda,q,nt,nt1=round(rt*nt),n=nt1,type='inclusive',nrep=1e3,method='rejection')
{
  n_component <- length(lambda)

  p=rep(0,n)
  for (i in 1:nrep) {
    if (method=='rejection') {
      ks=0
      while (sum(ks)!=nt1){
        m=rmultinom(1, nt, q)
        poi_rates=rep(lambda, m)
        ks=rpois(nt, poi_rates)
      }
    }
    if (method=='multinomial') {
      ks=rmultinom(1, nt1, rep(lambda, rmultinom(1, nt, q)))
    }
    if (type=='inclusive'){
      p[1]=p[1]+1
    } else {
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
compare_simu_mix_pois_2=function(lambda,q,nt=10,nt1=20,n=15,type='inclusive',nrep=1e3,method='multinomial')
{
  # s=simul_pcoal(rt=rt,vt=vt,nt=nt,nt1=nt1,n=n,type=type,nrep=nrep,method=method)
  s=simul_pcoal_mixed_pois(lambda,q,nt,nt1=round(rt*nt),n=nt1,type='inclusive',nrep=1e3,method=method)
  if (type=='inclusive') {
    p=sapply(1:nt1, function(x) mixed_pois_2_inclusive(x, nt, nt1, lambda, q))
  } else {
    p=sapply(1:nt1, function(x) mixed_pois_2_exclusive(x, nt, nt1, lambda, q))
  }

  ret=cbind(s,p,s-p)
  colnames(ret)=c('Simulated','Calculated','Difference')
  return(ret)
}
