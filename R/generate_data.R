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
#' # generate a two-component 21-dimentional data
#' G=2 # number of component
#' K=20 # number of taxa - dimension of the latent variabel
#'
#' n_1 = 250
#' t_mu_1 <- c(runif(20,0.5,1), 0)
#' t_Sigma_1 <- rbind(cbind(genPositiveDefMat(dim=20,covMethod="unifcorrmat",rangeVar=c(0.25,0.75))$Sigma,0),0)
#'
#' n_2 = 250
#' t_mu_2 <- c(runif(20,-1,0.2), 0)
#' t_Sigma_2 <- rbind(cbind(genPositiveDefMat(dim=20,covMethod="unifcorrmat",rangeVar=c(0,0.5))$Sigma,0),0)
#'
#' num_observation=c(n_1,n_2)
#' true_mu<-list(t_mu_1,t_mu_2)
#' true_Sig<-list(t_Sigma_1,t_Sigma_2)
#'
#' true_par <- list(G=G,N=num_observation,K=K,mu=true_mu,Sigma=true_Sig)
#' Data.temp <- generate_data(G=2,num_observation=c(n_1,n_2),K=20,true_mu=list(t_mu_1,t_mu_2),true_Sig=list(t_Sigma_1,t_Sigma_2),seed.no=1234,M=10000,truelab=TRUE)

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
