#' Another function for newton update of V
#' This function calculate the second derivative of the elbo with respect to V for it's kth diagonal element.
#' @param x kth diagonal element of V, variable of this function
#' @param Sig_kk [k,k]th element of the variance parameter for the latent Gaussian variable,.
#' @param m_k kth element of variational parameter m.
#' @param M sum of counts in K taxa for one observation.
#' @param xi variational parameter xi, vector of n.
#' @keywords prime
#' @export
#' @examples
#' vk_fun.prime()

vk_fun.prime <- function(x,Sig_kk,m_k,M,xi){
  return(-1/x^2-1/Sig_kk-(x^2+1)*M*exp(m_k+x^2/2)/xi)
}
