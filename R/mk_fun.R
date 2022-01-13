#' A function for newton update of m
#' This function calculate the first derivative of the elbo with respect to m for it's kth element.
#' @param x kth element of m, variable of this function
#' @param W_i ith observation.
#' @param Sig variance parameter for the latent Gaussian variable.
#' @param mu mean parameter for the latent Gaussian variable.
#' @param vsquare diagonal elements of the variational parameter V.
#' @param M sum of counts in K taxa for one observation.
#' @param xi variational parameter xi, vector of n.
#' @keywords mk
#' @export
#' @examples
#' mk_fun()

mk_fun <- function(x,W_i,Sig,mu,vsquare,M,xi){
  return(W_i-MASS::ginv(Sig)%*%(x-mu)-M*exp(x+vsquare/2)/xi)
}
