#' A function for newton update of V
#' This function calculate the first derivative of the elbo with respect to V for it's kth diagonal element.
#' @param x kth diagonal element of V, variable of this function
#' @param iSig_kk k-th diagonal element of the (generalized) inverse variance parameter for the latent Gaussian variable.
#' @param m_k kth element of variational parameter m.
#' @param M sum of counts in K taxa for one observation.
#' @param xi variational parameter xi, vector of n.
#' @keywords vk
#' @export
#' @examples
#' vk_fun()

vk_fun <- function(x,iSig_kk,m_k,M,xi){
  return(1/x-x*iSig_kk-x*M*exp(m_k+x^2/2)/xi)
}
