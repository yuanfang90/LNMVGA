#' A function for newton update of V
#' This function performs newton update of V within one iteration.
#' @param xold kth diagonal element of V, value from last iteration.
#' @param Sig_kk [k,k]th element of the variance parameter for the latent Gaussian variable,.
#' @param m_k kth element of variational parameter m.
#' @param M sum of counts in K taxa for one observation.
#' @param xi variational parameter xi, vector of n.
#' @keywords vk
#' @export
#' @examples
#' newton_up_v()

newton_up_v <- function(xold,Sig_kk,m_k,M,xi){
  xnew <- xold-vk_fun(xold,Sig_kk,m_k,M,xi)/vk_fun.prime(xold,Sig_kk,m_k,M,xi)
  return(xnew)
}
