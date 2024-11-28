#' Inclusive probability of coalescence of k lineages when the offspring distribution is Poisson
#'
#' @param k Number of lineages to coalesce
#' @param nt Population size at time t
#' @param log Whether to return the log of the probability
#'
#' @return Probability of coalescence
#' @export
#'
pois_inclusive=function(k,nt,log=F) {
  if (log==F)
    return(1/nt^(k-1))
  else
    return(-(k-1)*log(nt))
}

#' Exclusive probability of coalescence of k lineages when the offspring distribution is Poisson
#'
#' @param k Number of lineages to coalesce
#' @param n Number of observed lineages
#' @param nt Population size at time t
#' @param log Whether to return the log of the probability
#'
#' @return Probability of coalescence
#' @export
#'
pois_exclusive=function(k,n,nt,log=F) {
  if (log==F)
    return((nt-1)^{n-k}/nt^(n-1))
  else
    return((n-k)*log(nt-1)-(n-1)*log(nt))
}

#' Inclusive probability of coalescence of k lineages when the offspring distribution is Negative-Binomial
#'
#' @param k Number of lineages to coalesce
#' @param nt Population size at time t
#' @param r Dispersion parameter of Negative-Binomial
#' @param log Whether to return the log of the probability
#'
#' @return Probability of coalescence
#' @export
#'
negbin_inclusive=function(k,nt,r,log=F) {
  if (k==1 && log==F) return(1)
  if (k==1 && log==T) return(0)
  if (log==F)
    return(prod(r+seq(1,k-1))/prod(nt*r+seq(1,k-1)))
  else
    return(sum(log(r+seq(1,k-1)))-sum(log(nt*r+seq(1,k-1))))
}

#' Exclusive probability of coalescence of k lineages when the offspring distribution is Negative-Binomial
#'
#' @param k Number of lineages to coalesce
#' @param n Number of observed lineages
#' @param nt Population size at time t
#' @param r Dispersion parameter of Negative-Binomial
#' @param log Whether to return the log of the probability
#'
#' @return Probability of coalescence
#' @export
#'
negbin_exclusive=function(k,n,nt,r,log=F) {
  # p=nt * prod(r + 0:(k-1)) * prod((nt-1) * r + 0:(n-k-1)) / prod(nt * r + 0:(n-1))
  p=nt * beta(k+r, n-k+(nt-1)*r) / beta(r, (nt-1)*r)

  if (log) return(log(p))
  else return(p)
}
