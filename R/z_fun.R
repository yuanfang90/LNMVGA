#' A function to calculate posterior clustering probabilities
#'
#' This function calculate z, a matrix with each row representing the probability of the observations come from the corresponding group.
#' @param W observation data: the count data, n*(K+1) matrix where K is the number of taxa.
#' @param xi variational parameter xi, vector of n.
#' @param m variational parameter m, vector of n.
#' @param V variational parameter V, list of n matrices, each element is a (K+1)*(K+1) matrix.
#' @param mu mean parameter for the latent Gaussian variable, list of g vectors of (K+1).
#' @param Sig variance parameter for the latent Gaussian variable, list of g (K+1)*(K+1) matrices.
#' @param pi_g proportion vector of each component.
#' @param G number of component.
#' @keywords z_clus
#' @export
#' @examples
#' z_fun()

z_fun <- function(W,xi,m,V,mu,Sig,pi_g,G){
  n=nrow(W)
  K=ncol(W)-1
  condition <- 1
  W_nmlzd <- W/mean(rowSums(W))

  forz <- matrix(nrow=n,ncol=G)
  z_clus <- forz
  for(i in 1:n){
    for(g in 1:G){
      elbo <- elbo_fun(W_i=W_nmlzd[i,],xi_i=xi[i],m_i=m[i,],V_i=V[[i]],mu_g=mu[g,],Sig_g=Sig[[g]])
      # print(elbo)
      if(!is.na(elbo)){
        kk <- try(pi_g[g]*exp(elbo),silent=TRUE)
        # print(kk)
        # if(is.numeric(kk)&!class(kk) == "try-error"){
        if(!any(class(kk) == "try-error")){
          forz[i,g]<-kk
        }else{
          forz[i,g] <- -Inf
          condition <- 0
        }
      }else{
        forz[i,g] <- NA
        condition <- 0
      }
    }
    # if(sum(forz[i,])==0){forz[i,]<-pi_g*.Machine$double.neg.eps}
    # if(sum(forz[i,])==0){forz[i,]<-NA}
  }
  if (condition == 1) {
    z_clus <- forz/rowSums(forz)
  }
  return(list(z=z_clus,forz=forz,con=condition))
}

