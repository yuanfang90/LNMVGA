#' Function for generating data from the mixture of logistic-normal multinomial model
#'
#' This function generates count data from the mixture of logistic-normal multinomial model. It assume the specication of the mean and variance parameter for the latent variable, generate Gaussian latent variable data first, then transform back to compositional data, and then count data.
#' @param num_grp Number of groups that you want to generate the data from.
#' @param num_observation A vector with length as num_grp, specifying the number of observations inside each group.
#' @param K Dimension of the latent variable. Number of "taxa" - 1. For example, if you want to generate counts from 4 taxa, specify K=3.
#' @param true_mu True means of the latent variable. A list contains num_grp vectors. Each vector is of length K+1, with the last element equals 0.
#' @param true_Sig True covariance matrices of the latent variable. A list contains num_grp matrices. Each matrix is a K*K positive definite matrix attached by one row of 0 and one column of 0.
#' @param seed.no set.seed() argument. Setting a seed for retriving the data.
#' @param M A number to generate total observed counts. The total counts are generated as a uniform random numbers in the interval from M/2 to M.
#' @param truelab Logical variable, if is TRUE, the true label of each observation will be returned. Default is TRUE.
#' @return A list cotains both counts data (W) and the latent variable (Y). If truelab is TRUE, true class labels are returned (true_lab).
#' @export
#' @keywords data
#' @examples
#' # generate a two-component 3-dimentional data
#' generate_data(G = 2, num_observation = c(50,50), K = 2, true_mu = list(c(0,1,0),c(-2,-5,0)),true_Sig=list(rbind(cbind(diag(1,2),0),0),rbind(cbind(diag(1,2),0),0)), seed.no = 1234, M = 10000, truelab = TRUE)

generate_data <- function(num_grp,num_observation,K,true_mu,true_Sig,seed.no,M,truelab=TRUE){
  set.seed(seed.no)
  taxa_count_ls <- list()
  Y_ls <- list()
  lab <- list()
  for(g in 1:num_grp){
    Y <- mvtnorm::rmvnorm(num_observation[g],mean=true_mu[[g]],sigma=true_Sig[[g]])
    Y_ls[[g]] <-Y
    eta <- Y
    theta <- exp(eta)/rowSums(exp(eta)) #underlying taxa composition, softmax of Y/eta
    taxa_count_ls[[g]]<-matrix(nrow=num_observation[g],ncol=(K+1))
    for(i in 1:num_observation[g]){
      taxa_count_ls[[g]][i,] <- rmultinom(1,sample(c(M/2):M,1),prob=theta[i,])
    }
    lab[[g]] <- rep(g,num_observation[g])
  }
  taxacount <- Reduce(rbind,taxa_count_ls)
  hidden_y <- Reduce(rbind,Y_ls)
  lab <- Reduce(c,lab)

  if(truelab==TRUE){
    return(list(W=taxacount,Y=hidden_y,true_lab=lab))
  }else{
    return(list(W=taxacount,Y=hidden_y))
  }
}
