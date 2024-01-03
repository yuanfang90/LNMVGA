#' A function for newton update of m
#' This function performs newton update of m within one iteration.
#' @param xold kth element of m, value from last iteration.
#' @param W_i ith observation.
#' @param iSig (generalized) inverse matrices variance parameter for the latent Gaussian variable.
#' @param mu mean parameter for the latent Gaussian variable.
#' @param vsquare diagonal elements of the variational parameter V.
#' @param M sum of counts in K taxa for one observation.
#' @param xi variational parameter xi, vector of n.
#' @param K number of taxa.
#' @keywords newton_up_m
#' @export
#' @examples
#' newton_up_m()

newton_up_m <- function(xold,W_i,iSig,mu,vsquare,M,xi,K){
  xnew <- c(xold-solve(mk_fun.jacobian(xold,iSig,vsquare,M,xi),tol=.Machine$double.neg.eps)%*%mk_fun(xold,W_i,iSig,mu,vsquare,M,xi))
  xnew[(K+1)] <- 0
  return(xnew)
}
