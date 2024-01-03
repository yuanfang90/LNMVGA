#' An internal function to calculate the expected log sum exponential term for calculation of the elbo
#' This function takes in the parameters for all groups and output a matrix, with the element from the ith row and gth column representing the elbo calculated for the ith observation under the condition that it comes from the gth group.
#' @param W count data
#' @param m variational parameter m
#' @param V variational parameter V, a list of n matrices, n being number of observations.
#' @param Vmat internal parameter, a matrix of n rows as n being number of observations; each row is the diagonal element of the variational parameter V.
#' @param mu mean parameter for the latent Gaussian variable.
#' @param Sig variance parameter for the latent Gaussian variable, a list of G (K+1)*(K+1) matrices, G is the numbe of component and K+1 is the number of taxa (dimension of the count data).
#' @param red_sig internal parameter, a list of G K*K dimensional matrices from the variance parameter.
#' @param G number of components.
#' @keywords elbo
#' @export
#' @examples
#' elbo_fun()

elbo_fun <- function(W,m,Vmat,V,mu,red_sig,Sig,G){
  n=nrow(W)
  temp_ELBO<-matrix(NA,nrow=n,ncol=G)
  for (g in 1:G){
    first_e<-diag(W%*%t(m))
    second_e<-log(rowSums(exp(m+Vmat/2)))*rowSums(W)
    third_e<-0.5*log(det(red_sig[[g]]))
    fourth_e<-0.5*mahalanobis(m,center=mu[g,],cov=ginv(Sig[[g]]),inverted=TRUE)
    fifth_fun<-function(V_hatm){
      sum(diag(ginv(Sig[[g]])%*%V_hatm))
    }
    fifth_e<-0.5*sapply(V,fifth_fun)
    sixth_e<-0.5*rowSums(log(Vmat[,-(K+1)]))
    temp_ELBO[,g]<-first_e-second_e-third_e-fourth_e-fifth_e+sixth_e+K/2+multi_const
  }
  return(temp_ELBO)
}
