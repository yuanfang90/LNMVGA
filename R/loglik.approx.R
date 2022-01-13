#' A function to calculate the approximated log-likelihood
#'
#' This function calculate approximated log-likelihood using the elbo.
#' @param W observation data: the count data, n*(K+1) matrix where K is the number of taxa.
#' @param xi variational parameter xi, vector of n.
#' @param m variational parameter m, vector of n.
#' @param V variational parameter V, list of n matrices, each element is a (K+1)*(K+1) matrix.
#' @param mu mean parameter for the latent Gaussian variable, list of g vectors of (K+1).
#' @param Sig variance parameter for the latent Gaussian variable, list of g (K+1)*(K+1) matrices.
#' @param pi_g proportion vector of each component.
#' @param G number of component.
#' @keywords loglik
#' @export
#' @examples
#' loglik.approx()

loglik.approx <- function(W,xi,m,V,mu,Sig,pi_g,G){
  n=nrow(W)
  K=ncol(W)-1

  forelbo <- matrix(nrow=n,ncol=G)
  for(i in 1:n){
    for(g in 1:G){
      elbo <- elbo_fun(W_i=W[i,],xi_i=xi[i],m_i=m[i,],V_i=V[[i]],mu_g=mu[g,],Sig_g=Sig[[g]])
      # print(elbo)
      if(!is.na(elbo)){
        kk <- pi_g[g]*exp(elbo)
        # print(kk)
        forelbo[i,g]<-kk
      }else{
        forelbo[i,g]<-NA
      }
    }
    # if(sum(forelbo[i,])==0){forelbo[i,]<-pi_g*.Machine$double.neg.eps}
    # if(sum(forelbo[i,])==0){forelbo[i,]<-NA}
  }
  return(sum(log(rowSums(forelbo,na.rm=TRUE))))
}

