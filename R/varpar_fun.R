#' A function to update the variational parameters.
#'
#' This function use the newton update realted function and perform update on the variational parameters in one iteration.
#' @param W observation data: the count data, n*(K+1) matrix where K is the number of taxa.
#' @param m variational parameter m from last iteration.
#' @param V variational parameter Vfrom last iteration.
#' @param mu mean parameter for the latent Gaussian variable, list of g vectors of (K+1).
#' @param iSig (generalized) inverse variance parameter for the latent Gaussian variable, list of g (K+1)*(K+1) matrices.
#' @param K number of taxa.
#' @param z class membership matrix at the current iteration.
#' @keywords varpar
#' @export
#' @examples
#' varpar_fun()

varpar_fun <- function(W,m,V,mu,iSig,K,z){
  n=nrow(W)
  xi_hat <- numeric(n)
  m_hat <- matrix(NA,nrow=n,ncol=K+1)
  V_hat <- vector("list",n)
  Vmat<-matrix(NA,nrow=n,ncol=K+1)
  for(i in 1:n){
    g=which.max(z[i,])
    ## update xi
    xi_hat[i] <- sum(exp(m[i,]+diag(V[[i]])/2))
    ## update m
    m_hat[i,] <- newton_up_m(xold=m[i,],W_i=W[i,],iSig=iSig[[g]],mu=mu[g,],vsquare=diag(V[[i]]),M=sum(W[i,]),xi=xi_hat[i],K=K)
    ## update V
    for_V <- rep(0,(K+1))
    for(k in 1:K){
      for_V[k] <- newton_up_v(xold=sqrt(diag(V[[i]])[k]),iSig_kk=iSig[[g]][k,k],m_k=m[i,k],M=sum(W[i,]),xi=xi_hat[i])
    }
    V_hat[[i]] <- diag((for_V)^2)
    Vmat[i,]<-diag(V_hat[[i]])
  }

  return(list(xi=xi_hat,m=m_hat,V=V_hat,Vmat=Vmat))
}
