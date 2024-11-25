#' Probability of coalescence of two lineages when the offspring distribution is Poisson
#'
#' @param nt Population size at time t
#' @param log Whether to return the log of the probability
#'
#' @return Probability of coalescence
#' @export
#'
pois_p2t=function(nt,log=F) {
  if (log==F)
    return(1/nt)
  else
    return(-log(nt))
}

#' Probability of coalescence of n lineages when the offspring distribution is Poisson
#'
#' @param n Number of lineages
#' @param nt Population size at time t
#' @param log Whether to return the log of the probability
#'
#' @return Probability of coalescence
#' @export
#'
pois_pnt=function(n,nt,log=F) {
  if (log==F)
    return(1/nt^(n-1))
  else
    return(-(n-1)*log(nt))
}

#' Probability of coalescence of two lineages when the offspring distribution is Negative-Binomial
#'
#' @param nt Population size at time t
#' @param r Dispersion parameter of Negative-Binomial
#' @param log Whether to return the log of the probability
#'
#' @return Probability of coalescence
#' @export
#'
negbin_p2t=function(nt,r,log=F) {
  if (log==F)
    return((r+1)/(nt*r+1))
  else
    return(log(r+1)-log(nt*r+1))
}

#' Probability of coalescence of n lineages when the offspring distribution is Negative-Binomial
#'
#' @param n Number of lineages
#' @param nt Population size at time t
#' @param r Dispersion parameter of Negative-Binomial
#' @param log Whether to return the log of the probability
#'
#' @return Probability of coalescence
#' @export
#'
negbin_pnt=function(n,nt,r,log=F) {
  if (log==F)
    return(prod(r+seq(1,n-1))/prod(nt*r+seq(1,n-1)))
  else
    return(sum(log(r+seq(1,n-1)))-sum(log(nt*r+seq(1,n-1))))
}
