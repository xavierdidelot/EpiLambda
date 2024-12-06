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
  #prod(r+seq(1,k-1))/prod(nt*r+seq(1,k-1))
  if (log==F)
    return(beta(nt*r+1,r+k)/beta(r+1,nt*r+k))
  else
    return(lbeta(nt*r+1,r+k)-lbeta(r+1,nt*r+k))
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
  if (log==F)
    return(nt * beta(k+r, n-k+(nt-1)*r) / beta(r, (nt-1)*r))
  else
    return(log(nt) + lbeta(k+r, n-k+(nt-1)*r) - lbeta(r, (nt-1)*r))
}

#' Inclusive probability of coalescence of k lineages when the offspring distribution is a 2-component Poisson mixture
#'
#' @param k Number of lineages to coalesce
#' @param nt Population size at time t
#' @param nt1 Population size at time t+1
#' @param lambda 2-dimensional vector of Poisson rates
#' @param q Probability of each emmission rate
#' @param log Whether to return the log of the probability
#'
#' @return Probability of coalescence
#' @export
#'
mixed_pois_2_inclusive=function(k,nt,nt1,lambda,q,log=F) {
  sum_poisson_mix <- function(n, N, lambda, q){ #P[X1 + ... + Xn = N]
    sum(sapply(0:n, function(m) choose(n, m) * (q[1] * exp(-lambda[1]))^m * (q[2] * exp(-lambda[2]))^(n-m) * (m * lambda[1] + (n-m) * lambda[2])^N)) / factorial(N)
  }

  conditional_dist <- function(n, N, lambda, q){ #P[X1 = k | X1 + ... + Xn = N] for k=0:N
    con_dist <- numeric(N+1)
    for (k in 0 : N){
      con_dist[k+1] <- (q[1] * dpois(k, lambda[1]) + q[2] * dpois(k, lambda[2])) * sum_poisson_mix(n-1, N-k, lambda, q) / sum_poisson_mix(n, N, lambda, q)
    }

    return(con_dist)
  }

  p <- conditional_dist(nt, nt1, lambda, q)

  if (log==F){
    return(nt * sum(choose(0:nt1, k) * p) / choose(nt1, k))
  } else {
    return(log(nt) + log(sum(choose(0:nt1, k) * p)) - log(choose(nt1, k)))
  }
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
mixed_pois_2_exclusive=function(k,n,nt,log=F) {
  stop('Mixed Poisson exclusive probability not added yet')
}


