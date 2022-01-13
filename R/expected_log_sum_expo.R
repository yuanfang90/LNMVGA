#' A function to calculate the expected log sum exponential term for calculation of the elbo
#'
#' This function calculate the expected log sum exponential term which will be used in the calculation of the elbo.
#' @param xi_i variational parameter xi, for the ith observation.
#' @param m_i variational parameter m, for the ith observation.
#' @param V_i variational parameter V, for the ith observation.
#' @keywords sum_expo
#' @export
#' @examples
#' expected_log_sum_expo()

expected_log_sum_expo <- function(xi_i,m_i,V_i){
  return(sum(exp(m_i+diag(V_i)/2))/xi_i+log(xi_i)-1)
}
