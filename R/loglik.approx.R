#' A function to calculate the approximated log-likelihood
#'
#' This function calculate approximated log-likelihood using the elbo.
#' @param W observation data: the count data, n*(K+1) matrix where K is the number of taxa.
#' @param xi variational parameter xi, vector of n.
#' @param m variational parameter m, vector of n.
#' @param V variational parameter V, list of n matrices, each element is a (K+1)*(K+1) matrix.
#' @param Vmat internal parameter, a matrix of n rows as n being number of observations; each row is the diagonal element of the variational parameter V.
#' @param mu mean parameter for the latent Gaussian variable, list of g vectors of (K+1).
#' @param Sig variance parameter for the latent Gaussian variable, list of g (K+1)*(K+1) matrices.
#' @param red_sig internal parameter, a list of G K*K dimensional matrices from the variance parameter.
#' @param pi_g proportion vector of each component.
#' @param G number of component.
#' @param it number of iteration, keep track of iterations for condition checking
#' @keywords loglik
#' @export
#' @examples
#' loglik.approx()

loglik.approx <- function(W,m,V,Vmat,mu,Sig,red_sig,pi_g,G,it){
  n=nrow(W)

  for_elbo <- elbo_fun(W=W,m=m,Vmat=Vmat,V=V,mu=mu,red_sig=red_sig,Sig=Sig,G=G)

  pi_mat<-matrix(pi_g,ncol=G,nrow=n,byrow=TRUE)
  for_tau<-pi_mat*exp(for_elbo)
  for_ll<-rowSums(for_tau)

  if (it<5){
    for_ll[for_ll==0]<-1
  }

  return(sum(log(for_ll)))
}

