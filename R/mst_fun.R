#' A function to update the model parameters.
#'
#' This function to update the model parameters for the latent Gaussian variables.
#' @param m variational parameter m from last iteration.
#' @param V variational parameter Vfrom last iteration.
#' @param G number of component.
#' @param K number of taxa.
#' @param z class membership matrix at the current iteration.
#' @keywords mst
#' @export
#' @examples
#' mst_fun()
#'
mst_fun <- function(m,V,G,K,z){
  mu <- matrix(NA,nrow=G,ncol=(K+1))
  Sig <- list()
  for(g in 1:G){
    temp <- try(cov.wt(m, wt=(z[,g]), center=TRUE, method="ML"), silent=TRUE)
    if(!class(temp) == "try-error"){
      mu[g,] <- temp$center
      forsig <- Reduce(`+`,Map(`*`,V,(z[,g]/sum(z[,g]))))
      Sig[[g]] <- temp$cov+forsig
    }
  }
  return(list(mu=mu,Sig=Sig))
}
