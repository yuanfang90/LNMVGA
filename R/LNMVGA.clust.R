#' Function for the Clustering Algorithm
#'
#' This function is the main clustering function for the proposed algorithm.
#' @param data Matrix of integers. Input data here, data are counts.
#' @param G Input a vector specifying all number of components that the model is fit for. E.g., G = c(1:4) will fit models for G = 1,  2, 3, 4.
#' @param initial Specify method for initializing z_ig. Possible values could be "kmeans", "random", "small_EM".
#' @param runtime Logical variable, if outputting the running time of the whole procedure or not.
#' @param true_lab Vector of true class membership. If not specified, defalt is NA for all observations.
#' @param threshold Threshold for the Atiken's stopping creterion for convergence.
#' @param verb Logical variable, if the key steps of the algortihm and approximated loglikelihood for each iteration are printed.
#' @param maxiter Maximum number of iteration. If specified, algorithm will stop by either below the threshold or maxiter reached. If not specified, algorithm will only be monitored by convergence criterion.
#' @param nrep Default is NA. Only needed if "small_EM" is specified for initial. Number of random starts for the small EM initialization.
#' @param niter Default is NA. Only needed if "small_EM" is specified for initial. Number of iterations for each random start in the small EM initialization.
#' @return  A list contains the models selected by BIC or ICL, classification results, BIC/ICL values for all fitted models, and parameter estimation for the selcted model as a vector. Specicially, result.BIC or result.ICL outputs a list contains results for model selected by BIC or ICL: class=class membership, zhat= soft class membership posterior probability, mu=estimated mean parameter for the latent variable, Sigma=estimated variance paprameter for the latent variable, G=number of component selected by BIC, BIC=BIC for all G fitted. est.vec.BIC or est.vec.ICL output a vector contains estimated vectorized mean and variance parameters for the latent variable, and number of component selected. runtime Only output if "runtime == TRUE" in the input, the time for finishing the clustering algorithm for all observations fitting models 1:G.
#' @export
#' @keywords clust
#' @examples
#' LNMVGA.clust()

LNMVGA.clust <- function(data,G,initial,runtime,true_lab=rep(NA,nrow(data)),threshold,verb,maxiter=NA,nrep=NA,niter=NA){
  ##### data: input data, counts, integer
  ##### G: input vector of all G, specifying all G's that the model is fit for
  ##### initial: possible values could be "kmeans", "random", "small_EM", different way of initialize z_ig
  ##### select: possible values could be "BIC" and "ICL", specifying different model selection creteria
  ##### runtime: logical, if outputting the running time of the whole procedure or not
  ##### true_lab: defaul is NA. if there is any, provide the vector contains the true label, ARI will be calculated and included in the vector estimator.
  W <- data
  W_init <- W
  if(any(W_init==0)) W_init[which(W_init==0)]<-1

  n = nrow(W)
  K=ncol(W)-1

  theta <- W_init/rowSums(W_init)
  latent <- log(theta/theta[,(K+1)])
  # pairs(latent,col=true_lab)
  # latentmean <- matrix(NA,nrow=G,ncol=K+1)
  # latentcov <- list()
  # for(g in 1:G){
  #   temp <- cov.wt(latent,wt=(ztrue[,g]),center=TRUE,method="ML") # weighted mean and cov
  #   latentmean[g, ] <- temp$center
  #   latentcov[[g]] <- temp$cov
  # }

  all_G <- G
  results <- list()
  BIC <- NULL
  ICL <- NULL

  start_time <- Sys.time()
  for(G_act in all_G){
    ###########################################################
    #####Initialization of group membership and parameters#####
    ###########################################################
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
    old_par<-list(pi_g=pi_g,xi=par$xi,m=par$m,V=par$V,mu=par$mu,Sig=par$Sig,z=z)

    ###########################################################
    #######################EM loop starts######################
    ###########################################################
    loglik_approx <- NULL
    loglik_approx[1] <- -Inf
    bohning_asmp_lik <- NULL
    bohning_asmp_lik[1] <- -Inf
    bohning_asmp_lik[2] <- -Inf

    it <- 2
    stop <- 0
    while(stop < 1){
      # E-step update z,alpha,xi,m,V
      tau <- z_fun(W=W,xi=old_par$xi,m=old_par$m,V=old_par$V,mu=old_par$mu,Sig=old_par$Sig,pi_g=old_par$pi_g,G=G_act)
      # if(tau$con == 1){
      #   zhat <- tau$z
      # }else{
      #   if(it==2){zhat <- z}else{zhat <- zhat}
      # }
      if(tau$con == 0){zhat <- old_par$z}else{zhat <- tau$z}
      varpar_up <- varpar_fun(W=W,m=old_par$m,V=old_par$V,mu=old_par$mu,Sig=old_par$Sig,K=K,z=zhat)
      xi_hat <- varpar_up$xi
      m_hat <- varpar_up$m
      V_hat <- varpar_up$V

      # M-step
      pi_g_hat <- colSums(zhat)/sum(zhat)
      musig_up <- mst_fun(m=m_hat,V=V_hat,G=G_act,K=K,z=zhat)
      if(all(!is.na(musig_up$mu))){
        mu_hat <- musig_up$mu
        Sig_hat <- musig_up$Sig
      }else{
        mu_hat <- old_par$mu
        Sig_hat <- old_par$Sig
      }

      # calculate approximated loglikelihood
      loglik_approx[it] <- loglik.approx(W=W,xi=xi_hat,m=m_hat,V=V_hat,mu=mu_hat,Sig=Sig_hat,pi_g=pi_g_hat,G=G_act)

      if(verb==TRUE){
        cat(paste("G_act =",G_act))
        cat(paste("iteration",it))
        print(pi_g_hat)
        # print(mu_hat)
        # print(Sig_hat)
        print(loglik_approx[it])
      }

      if(it > 2){
        a_it<-(loglik_approx[it]-loglik_approx[(it-1)])/(loglik_approx[(it-1)]-loglik_approx[(it-2)])
        l_asmp<-loglik_approx[(it-1)]+(loglik_approx[(it-1)]-loglik_approx[(it-2)])/(1-a_it)
        bohning_asmp_lik[(it-1)] <- l_asmp
      }
      if(verb==TRUE){
        print(abs(bohning_asmp_lik[(it-1)]-bohning_asmp_lik[(it-2)]))
      }

      maxit <- ifelse(is.na(maxiter),1000,maxiter)
      # if((it>10 && abs(bohning_asmp_lik[(it-1)]-bohning_asmp_lik[(it-2)])<=threshold)|it==maxit|loglik_approx[it]<loglik_approx[it-1]){stop = 1}
      # if(loglik_approx[it]<loglik_approx[it-1]){
      #   old_par <- old_par}else{
      #   old_par <- list(pi_g=pi_g_hat,xi=xi_hat,m=m_hat,V=V_hat,mu=mu_hat,Sig=Sig_hat,z=zhat)}

      if((it>10 && abs(bohning_asmp_lik[(it-1)]-bohning_asmp_lik[(it-2)])<=threshold)|it==maxit){stop = 1}
      old_par <- list(pi_g=pi_g_hat,xi=xi_hat,m=m_hat,V=V_hat,mu=mu_hat,Sig=Sig_hat,z=zhat)

      it <- it+1
    }
    # plot(loglik_approx)
    # plot(bohning_asmp_lik)
    # print(old_par$mu)
    # print(old_par$Sig)
    class <- mclust::map(old_par$z)
    # ARI <- mclust::adjustedRandIndex(class,true_lab)
    # pairs(Y,col=mclust::map(zhat),main=paste("ARI=",ARI))
    results[[G_act]] <- list(z=old_par$z,par=old_par)

    # Pamaters mu(K*G),Sig(K*K*G),pi_g(G-1)
    npar <- (K)*G_act+0.5*(K+1)*(K)*G_act+(G_act-1)
    BIC[G_act]<-2*loglik_approx[(it-1)]-npar*log(n)

    forz<-tau$forz
    mapz<-mclust::unmap(mclust::map(forz))
    # mapz<-mclust::unmap(class)
    forICL<-function(g){sum(log(forz[which(mapz[,g]==1),g]))}
    ICL[G_act]<-0.5*BIC[G_act]+sum(sapply(1:ncol(mapz),forICL))
  }
  end_time <- Sys.time()

  G_BIC_select <- which.max(BIC)
  G_ICL_select <- which.max(ICL)

  final.model.BIC <- G_BIC_select

  zhat.BIC <- results[[G_BIC_select]]$z
  pi_g_hat.BIC <- colSums(zhat.BIC)/sum(zhat.BIC)
  class.BIC <- mclust::map(zhat.BIC)

  estimate.BIC <- results[[G_BIC_select]]$par
  muhat.BIC <- estimate.BIC$mu
  Sighat.BIC <- estimate.BIC$Sig

  mu.BIC <- numeric((K+1)*G_BIC_select)
  Sig.BIC <- numeric((K+1)*(K+1)*G_BIC_select)
  for(g in 1:G_BIC_select){
    mu.BIC[(1+(K+1)*(g-1)):((K+1)*g)] <- muhat.BIC[g,]
    Sig.BIC[(1+(K+1)*(K+1)*(g-1)):((K+1)*(K+1)*g)] <- as.vector(Sighat.BIC[[g]])
  }

  final.model.ICL <- G_ICL_select

  zhat.ICL <- results[[G_ICL_select]]$z
  pi_g_hat.ICL <- colSums(zhat.ICL)/sum(zhat.ICL)
  class.ICL <- mclust::map(zhat.ICL)

  estimate.ICL <- results[[G_ICL_select]]$par
  muhat.ICL <- estimate.ICL$mu
  Sighat.ICL <- estimate.ICL$Sig

  mu.ICL <- numeric((K+1)*G_ICL_select)
  Sig.ICL <- numeric((K+1)*(K+1)*G_ICL_select)
  for(g in 1:G_ICL_select){
    mu.ICL[(1+(K+1)*(g-1)):((K+1)*g)] <- muhat.ICL[g,]
    Sig.ICL[(1+(K+1)*(K+1)*(g-1)):((K+1)*(K+1)*g)] <- as.vector(Sighat.ICL[[g]])
  }

  if(all(!is.na(true_lab))){
    ARI.BIC <- mclust::adjustedRandIndex(class.BIC,true_lab)
    pmvec.BIC <- c(mu.BIC,Sig.BIC,pi_g_hat.BIC,final.model.BIC,ARI.BIC)
    BIC.select <- list(class=class.BIC,zhat=zhat.BIC,mu=muhat.BIC,Sigma=Sighat.BIC,G=final.model.BIC,ARI=ARI.BIC,BIC=BIC)

    ARI.ICL <- mclust::adjustedRandIndex(class.ICL,true_lab)
    pmvec.ICL <- c(mu.ICL,Sig.ICL,pi_g_hat.ICL,final.model.ICL,ARI.ICL)
    ICL.select <- list(class=class.ICL,zhat=zhat.ICL,mu=muhat.ICL,Sigma=Sighat.ICL,G=final.model.ICL,ARI=ARI.ICL,ICL=ICL)
  }else{
    pmvec.BIC <- c(mu.BIC,Sig.BIC,pi_g_hat.BIC,final.model.BIC)
    BIC.select <- list(class=class.BIC,zhat=zhat.BIC,mu=muhat.BIC,Sigma=Sighat.BIC,G=final.model.BIC,BIC=BIC)

    pmvec.ICL <- c(mu.ICL,Sig.ICL,pi_g_hat.ICL,final.model.ICL)
    ICL.select <- list(class=class.ICL,zhat=zhat.ICL,mu=muhat.ICL,Sigma=Sighat.ICL,G=final.model.ICL,ICL=ICL)
  }

  if(runtime == TRUE){
    runtime = end_time - start_time
    return(list(result.BIC=BIC.select,result.ICL=ICL.select,est.vec.BIC=pmvec.BIC,est.vec.ICL=pmvec.ICL,runtime=runtime))
  }else{
    print("No running time recorded. If needed, specify runtim=TRUE. ")
    return(list(result.BIC=BIC.select,result.ICL=ICL.select,est.vec.BIC=pmvec.BIC,est.vec.ICL=pmvec.ICL))
  }

  # if(runtime == FALSE){
  #   print("No running time recorded. If needed, specify runtim=TRUE. ")
  #   return(list(final.G=final.model,class=class,zhat=zhat,mu=muhat,Sigma=Sighat,,est.vec=pmvec))
  # }
}
