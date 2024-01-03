#' A function to calculate posterior clustering probabilities
#'
#' This function calculate z, a matrix with each row representing the probability of the observations come from the corresponding group.
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
#' @keywords z_clus
#' @export
#' @examples
#' z_fun()

z_fun <- function(W,m,V,Vmat,mu,Sig,red_sig,pi_g,G,it){
  n=nrow(W)
  Kfor_elbo <- elbo_fun(W=W,m=m,Vmat=Vmat,V=V,mu=mu,red_sig=red_sig,Sig=Sig,G=G)

  pi_mat<-matrix(pi_g,ncol=G,nrow=n,byrow=TRUE)
  for_tau<-pi_mat*exp(for_elbo)
  for_ll<-rowSums(for_tau)

  if (it<5){
    for_ll[for_ll==0]<-1
    for_z<-pi_mat*exp(for_elbo-apply(for_elbo,1,max))
    tau<-for_z/rowSums(for_z)} else {
      tau<-for_tau/rowSums(for_tau)
    }
  zhat<-tau
  return(zhat)
}

