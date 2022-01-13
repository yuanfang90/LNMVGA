#' A Map Function
#'
#' This function use given posterior probabilities of observations belonging to a specific cluster (columns corresponding to clusters) to generate a vector of component labels.
#' @param x Input of the function is the vector of class membership posterior probabilities
#' @keywords which.max
#' @export
#' @examples
#' map()

map <- function(x) {
  return(apply(X=x, MARGIN=1, FUN=which.max))
}
