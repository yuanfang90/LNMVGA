#' A function to update the variational parameters.
#'
#' This function use the newton update realted function and perform update on the variational parameters in one iteration.
#' @param W observation data: the count data, n*(K+1) matrix where K is the number of taxa.
#' @param m variational parameter m from last iteration.
#' @param V variational parameter Vfrom last iteration.
#' @param mu mean parameter for the latent Gaussian variable, list of g vectors of (K+1).
#' @param Sig variance parameter for the latent Gaussian variable, list of g (K+1)*(K+1) matrices.
#' @param K number of taxa.
#' @param z class membership matrix at the current iteration.
#' @keywords varpar
#' @export
#' @examples
#' varpar_fun()

varpar_fun <- function(W,m,V,mu,Sig,K,z){
  n=nrow(W)
  xi_hat <- numeric(n)
  m_hat <- matrix(NA,nrow=n,ncol=K+1)
  V_hat <- vector("list",n)
  # forp <- colSums(z,na.rm=TRUE)/sum(z,na.rm=TRUE)
  for(i in 1:n){
    # g=ifelse(any(is.na(z[i,])),forp,which.max(z[i,]))
    g=which.max(z[i,])
    ## update xi
    xi_hat[i] <- sum(exp(m[i,]+diag(V[[i]])/2))
    ## update m
    m_hat[i,] <- newton_up_m(xold=m[i,],W_i=W[i,],Sig=Sig[[g]],mu=mu[g,],vsquare=diag(V[[i]]),M=sum(W[i,]),xi=xi_hat[i],K=K)
    ## update V
    for_V <- rep(0,(K+1))
    for(k in 1:K){
      for_V[k] <- newton_up_v(xold=sqrt(diag(V[[i]])[k]),Sig_kk=Sig[[g]][k,k],m_k=m[i,k],M=sum(W[i,]),xi=xi_hat[i])
    }
    V_hat[[i]] <- diag((for_V)^2)
  }
  return(list(xi=xi_hat,m=m_hat,V=V_hat))
}
