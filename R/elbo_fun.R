#' A function to calculate the expected log sum exponential term for calculation of the elbo
#'
#' This function calculate the expected log sum exponential term which will be used in the calculation of the elbo.
#' @param W_i count data the ith observation.
#' @param xi_i variational parameter xi, for the ith observation.
#' @param m_i variational parameter m, for the ith observation.
#' @param V_i variational parameter V, for the ith observation.
#' @param mu_g mean parameter for the latent Gaussian variable from the gth component.
#' @param Sig_g variance parameter for the latent Gaussian variable from the gth component.
#' @keywords elbo
#' @export
#' @examples
#' elbo_fun()

elbo_fun <- function(W_i,xi_i,m_i,V_i,mu_g,Sig_g){
  K=length(W_i)-1
  funforlb <- function(k){ifelse(W_i[k]==0,return(0),return(sum(log(1:W_i[k]))))}
  const <- sum(log(1:sum(W_i)))-sum(sapply(1:(K+1),funforlb))

  v_vec <- diag(V_i)[1:K]
  expctlogsumexpo <- expected_log_sum_expo(xi_i=xi_i,m_i=m_i,V_i=V_i)
  gamma_integral <- W_i%*%m_i - sum(W_i)*expctlogsumexpo
  entropy <- 0.5*(sum(log(v_vec))+K)
  expctlogpy <- 0.5*(-log(det(Sig_g[1:K,1:K]))-t(m_i-mu_g)%*%MASS::ginv(Sig_g)%*%(m_i-mu_g)-sum(diag(MASS::ginv(Sig_g)%*%V_i)))

  elbo <- gamma_integral+entropy+expctlogpy+const
  if(exp(elbo) != 0){
    return(elbo)
  }else{
    # return(log(.Machine$double.neg.eps))
    return(NA)
  }
}
