#' Probability of coalescence of two lineages when the offspring distribution is Poisson
#'
#' @param nt Population size at time t
#'
#' @return Probability of coalescence
#' @export
#'
pois_p2t=function(nt) {
  return(1/nt)
}

#' Probability of coalescence of n lineages when the offspring distribution is Poisson
#'
#' @param n Number of lineages
#' @param nt Population size at time t
#'
#' @return Probability of coalescence
#' @export
#'
pois_pnt=function(n,nt) {
  return(1/nt^(n-1))
}

#' Probability of coalescence of two lineages when the offspring distribution is Negative-Binomial
#'
#' @param nt Population size at time t
#' @param r Dispersion parameter of Negative-Binomial
#'
#' @return Probability of coalescence
#' @export
#'
negbin_p2t=function(nt,r) {
  return((r+1)/(nt*r+1))
}

#' Probability of coalescence of n lineages when the offspring distribution is Negative-Binomial
#'
#' @param n Number of lineages
#' @param nt Population size at time t
#' @param r Dispersion parameter of Negative-Binomial
#'
#' @return Probability of coalescence
#' @export
#'
negbin_pnt=function(n,nt,r) {
  return(prod(r+seq(1,n-1))/prod(nt*r+seq(1,n-1)))
}
