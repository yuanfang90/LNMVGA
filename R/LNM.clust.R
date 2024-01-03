#' Function for the Clustering Algorithm with specific number of component G_act
#'
#' This function is part of the main clustering function for the proposed algorithm. It fit the proposed LNM_MM model for a specific number of component G_act. If one wants to fit G for 1:5, use the paralleled function or put this function into a loop for G.
#' @param data Input data here. If sim==TRUE, data should be a list of multiple datasets (indexed by "run"), with each dataset as a list of counts W and true_lab (true class label). If no true label, set true_lab as NAs. If not simulation, data should be as the same format as described for each dataset of the simulation.
#' @param run Keep track of run number of datasets. For simulation this could be the index of the simulated data; for other cases, could run several times too with random initialization to pick the highest BIC/ICL. If only want to run 1 time for 1 dataset, specify run=1.
#' @param G_act Input the current actual running number of parameter.
#' @param initial Specify method for initializing z_ig. Possible values could be "kmeans", "random", "small_EM". Default is "kmeans".
#' @param runtime Logical variable, if outputting the running time of the whole procedure or not.
#' @param threshold Threshold for the Atiken's stopping creterion for convergence.
#' @param verb Logical variable, if the key steps of the algortihm and approximated loglikelihood for each iteration are printed.
#' @param maxiter Maximum number of iteration. If specified, algorithm will stop by either below the threshold or maxiter reached. If not specified, algorithm will only be monitored by convergence criterion.
#' @param nrep Default is NA. Only needed if "small_EM" is specified for initial. Number of random starts for the small EM initialization.
#' @param niter Default is NA. Only needed if "small_EM" is specified for initial. Number of iterations for each random start in the small EM initialization.
#' @param sim Indicator of whether this is simulated data. Simulated data input must as a list of multiple datasets (indexed by "run"), with each dataset must be a list of W and true_lab. Default is FALSE.
#' @return  A list contains the parameters when the algorithm converges. pi_g = estimated class size composition; z = soft class membership posterior probability; mu = estimated mean parameter for the latent variable; Sigma = estimated variance paprameter for the latent variable. Others are internal parameters that could used to check model fit and to pass to overall algorithm for model selection.
#' @export
#' @keywords clust
#' @examples
#' # generate data using Data.temp <- generate_data(G = 2, num_observation = c(50,50), K = 2, true_mu = list(c(0,1,0),c(-2,-5,0)),true_Sig=list(rbind(cbind(diag(1,2),0),0),rbind(cbind(diag(1,2),0),0)), seed.no = 1234, M = 10000, truelab = TRUE)
#'
#' LNM.clust(data=Data.temp,run=1,G_act=2,initial="small_EM",runtime=TRUE,threshold=1e-4,verb=TRUE,nrep=30,niter=50,sim=FALSE)

LNM.clust <- function(data,run,G_act,initial="kmeans",runtime=TRUE,threshold,verb=FALSE,maxiter=NA,nrep=NA,niter=NA,sim=FALSE){

  pt<-proc.time()

  if(sim==TRUE){
    W <- data[[run]]$W
    W_init <- W
    if(any(W_init==0)) W_init[which(W_init==0)]<-1

    n = nrow(W)
    K=ncol(W)-1

    true_lab<-data[[run]]$true_lab
  }else{
    W <- data$W
    W_init <- W
    if(any(W_init==0)) W_init[which(W_init==0)]<-1

    n = nrow(W)
    K=ncol(W)-1

    true_lab<-data$true_lab
  }

  theta <- W_init/rowSums(W_init)
  latent <- log(theta/theta[,(K+1)])

  results <- list()
  BIC <- NULL
  ICL <- NULL

  set.seed(run)
  ###Initialization
  if(initial == "kmeans"){
    initialz <- initialize_fun(W=W_init,G=G_act,K=K,n=n,method="kmeans")
  }
  if(initial == "small_EM"){
    nrep=nrep
    niter=niter
    if(verb==TRUE){
      cat("Be patient, small EM initialization is running...", "\n",sep="")
      initialz <- initialize_fun(W=W_init,G=G_act,K=K,n=n,method="small_EM",nrep=nrep,niter=niter,verb=TRUE)
    }else{
      initialz <- initialize_fun(W=W_init,G=G_act,K=K,n=n,method="small_EM",nrep=nrep,niter=niter,verb=FALSE)
    }
  }
  if(initial == "random"){
    initialz <- initialize_fun(W=W_init,G=G_act,K=K,n=n,method="random")
  }
  z <- initialz$z
  par <- initialz$par
  pi_g <- par$pi_g
  xi_hat=par$xi
  m_hat=par$m
  V_hat=par$V
  Vmat<-par$Vmat

  mu_hat=par$mu
  Sig_hat=par$Sig
  red_sig=par$red_sig
  iSig_hat=par$iSig

  #######################EM loop starts######################
  ###########################################################
  aloglik<-NULL
  loglik<-NULL
  aloglik[c(1,2,3,4,5)]<-0
  it_max<-ifelse(is.na(maxiter),500,maxiter)

  it <- 2
  stop <- 0
  zhat<-z
  while(stop < 1){

    ##Updating ELBO
    temp_ELBO<-matrix(NA,nrow=n,ncol=G_act)
    for (g in 1:G_act){
      first_e<-diag(W%*%t(m_hat))
      second_e<-log(rowSums(exp(m_hat+Vmat/2)))*rowSums(W)
      third_e<-0.5*log(det(red_sig[[g]]))
      fourth_e<-0.5*mahalanobis(m_hat,center=mu_hat[g,],cov=ginv(Sig_hat[[g]]),inverted=TRUE)
      fifth_fun<-function(V_hatm){
        sum(diag(ginv(Sig_hat[[g]])%*%V_hatm))
      }
      fifth_e<-0.5*sapply(V_hat,fifth_fun)
      sixth_e<-0.5*rowSums(log(Vmat[,-(K+1)]))
      temp_ELBO[,g]<-first_e-second_e-third_e-fourth_e-fifth_e+sixth_e+K/2+multi_const
    }
    pi_mat<-matrix(pi_g,ncol=G_act,nrow=n,byrow=TRUE)
    for_tau<-pi_mat*exp(temp_ELBO)
    for_ll<-rowSums(for_tau)
    if (it<5){
      for_ll[for_ll==0]<-1
      loglik[it]<-sum(log(for_ll))
      for_z<-pi_mat*exp(temp_ELBO-apply(temp_ELBO,1,max))
      tau<-for_z/rowSums(for_z)} else {
        loglik[it]<-sum(log(for_ll))
        tau<-for_tau/rowSums(for_tau)
      }
    zhat<-tau

    ########Update variational parameters
    varpar_up <- varpar_fun(W=W,m=m_hat,V=V_hat,mu=mu_hat,Sig=Sig_hat,K=K,z=zhat,iSig=iSig_hat)
    xi_hat <- varpar_up$xi
    m_hat <- varpar_up$m
    V_hat <- varpar_up$V
    Vmat<-varpar_up$Vmat

    ####    # M-step
    pi_g_hat <- colSums(zhat)/sum(zhat)
    musig_up <- mst_fun(m=m_hat,V=V_hat,G=G_act,K=K,z=zhat)
    mu_hat <- musig_up$mu
    Sig_hat <- musig_up$Sig

    iSig_hat<-list()
    for(g in 1:G_act){
      iSig_hat[[g]] <- MASS::ginv(Sig_hat[[g]])
    }
    red_sig<-list()
    for(g in 1:G_act){
      red_sig[[g]] <- Sig_hat[[g]][1:K,1:K]
    }

    if (it>5){
      #Aitkaine's stopping criterion
      if ((loglik[it-1]-loglik[it-2])==0) checks<-1 else{
        a<-(loglik[it]-loglik[it-1])/(loglik[it-1]-loglik[it-2])
        add_to<-(1/(1-a)*(loglik[it]-loglik[it-1]))
        # }
        aloglik[it]<-loglik[it-1]+add_to
        if (abs(aloglik[it]-aloglik[it-1])<=threshold) stop<-1 else stop<-stop
      }
    }
    if(verb==TRUE){
      cat(paste("G_act =",G_act))
      cat(paste("iteration",it))
    }

    it<-it+1
    if (it==it_max) stop<-1
    #print(it)
  }
  npar <- (K)*G_act+0.5*(K+1)*(K)*G_act+(G_act-1)
  BIC<-2*loglik[(it-1)]-npar*log(n)
  mapz<-mclust::unmap(mclust::map(zhat))
  forICL<-function(g){sum(log(zhat[which(mapz[,g]==1),g]))}
  ICL<-BIC+2*sum(sapply(1:ncol(mapz),forICL))

  pt2<-proc.time()

  if(runtime == TRUE){
    time_taken<-(pt2-pt)[3] ###In seconds
    return(list(pi_g=pi_g_hat,m=m_hat,V=V_hat,mu=mu_hat,Sig=Sig_hat,z=zhat,loglik=loglik,BIC=BIC,ICL=ICL,true=true_lab,time=time_taken,G=G,dataset=run))
  }else{
    print("No running time recorded. If needed, specify runtim=TRUE. ")
    return(list(pi_g=pi_g_hat,m=m_hat,V=V_hat,mu=mu_hat,Sig=Sig_hat,z=zhat,loglik=loglik,BIC=BIC,ICL=ICL,true=true_lab,G=G_act,dataset=run))
  }
}
