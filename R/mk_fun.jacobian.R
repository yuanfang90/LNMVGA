#' A function for newton update of m
#' This function calculate the jacobian matrix of the elbo with respect to m for it's kth element.
#' @param x kth element of m, variable of this function
#' @param Sig variance parameter for the latent Gaussian variable.
#' @param vsquare diagonal elements of the variational parameter V.
#' @param M sum of counts in K taxa for one observation.
#' @param xi variational parameter xi, vector of n.
#' @keywords jacobian
#' @export
#' @examples
#' mk_fun.jacobian()

mk_fun.jacobian <- function(x,Sig,vsquare,M,xi){
  return(-MASS::ginv(Sig)-M*diag(exp(x+vsquare/2))/xi)
}
