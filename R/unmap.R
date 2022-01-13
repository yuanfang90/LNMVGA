#' An  Unmap Function
#'
#' This function does the opposite to the map() function. It takes the class membership and unmap it to a matrix with entries that equal 1 if the component membership corresponded to that column. It will not drop empty column/class.
#' @param mapz Vector of class membership
#' @param G Number of component in assumed
#' @keywords mapz
#' @return A matrix with the [i,j]th entry equals 1 when the ith observation belongs to the jth group.
#' @export
#' @examples
#' unmap()

unmap <- function(mapz,G) {
  zunmap <- matrix(data=0, nrow=length(mapz), ncol=G)
  alphabet <- 1:G
  for (u in 1:G) {
    zunmap[which(mapz==alphabet[u]), u] <- 1
  }
  return(zunmap)
}
