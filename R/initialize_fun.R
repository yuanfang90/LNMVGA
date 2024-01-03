#' A function for initialization.
#'
#' This function performs initialization of the class membership providing the option of using k-means start, ramdom start, and a small EM or gridsearch type of start: initial randomly on a given number of random starting points, and run for a given number of iterations, then pick the one that ends with highest approximated log-likelihood.
#' @param W observation data: the count data, n*(K+1) matrix where K is the number of taxa.
#' @param G number of component.
#' @param K number of taxa.
#' @param n number of observations.
#' @param method string variable, could be "kmeans", "random", and "small_EM".
#' @param nrep number of random starts. Only useful if method = "small_EM" is specified. default is 50.
#' @param niter number of iterations for each small EM random start. Only useful if method = "small_EM" is specified. default is 10.
#' @param verb logical variable. if is TRUE, each iteration of the small_EM start will be printed.
#' @keywords initialize
#' @export
#' @examples
#' initialize_fun()

initialize_fun <- function(W,G,K,n,method,nrep=NA,niter=NA,verb=FALSE){
  ## kmeans initialization
  if(method == "kmeans"){
    theta <- W/rowSums(W)
    latent <- log(theta/theta[,(K+1)])
    z <- unmap(mapz=kmeans(latent,center=G,nstart=10)$cluster,G=G)

    pi_g <- colSums(z)/sum(z)
    xi_hat <- rep(1,n)
    m_hat <- latent
    V_hat <- vector("list",n)
    Vmat<-matrix(rep(c(rep(0.1,K),0),n),nrow=n,byrow=T)

    for(i in 1:n){V_hat[[i]]<-diag(Vmat[i,])}

    musig_up <- mst_fun(m=m_hat,V=V_hat,G=G,K=K,z=z)
    if(all(!is.na(musig_up$mu))){
      mu_hat <- musig_up$mu
      Sig_hat <- musig_up$Sig
    }else{
      mu_hat <- matrix(rep(colMeans(m_hat),G),nrow=G,ncol=(K+1),byrow=TRUE)
      Sig_hat <- list()
      for(g in 1:G){Sig_hat[[g]] <- rbind(cbind(diag(1,K),0),0)}
    }

    red_sig<-list()
    for(g in 1:G){
      red_sig[[g]] <- Sig_hat[[g]][1:K,1:K]
    }
    iSig_hat<-list()
    for(g in 1:G){
      iSig_hat[[g]] <- MASS::ginv(Sig_hat[[g]])
    }

    par <- list(pi_g=pi_g,xi=xi_hat,m=m_hat,V=V_hat,Vmat=Vmat,mu=mu_hat,Sig=Sig_hat,red_sig=red_sig,iSig=iSig_hat)

    return(list(z=z,par=par))
  }

  if(method == "random"){
    theta <- W/rowSums(W)
    latent <- log(theta/theta[,(K+1)])
    z <- unmap(mapz=sample(1:G,n,replace=TRUE),G=G)

    pi_g <- colSums(z)/sum(z)
    xi_hat <- rep(1,n)
    m_hat <- latent
    V_hat <- vector("list",n)
    Vmat<-matrix(rep(c(rep(0.1,K),0),n),nrow=n,byrow=T)

    for(i in 1:n){V_hat[[i]]<-diag(Vmat[i,])}

    musig_up <- mst_fun(m=m_hat,V=V_hat,G=G,K=K,z=z)
    if(all(!is.na(musig_up$mu))){
      mu_hat <- musig_up$mu
      Sig_hat <- musig_up$Sig
    }else{
      mu_hat <- matrix(rep(colMeans(m_hat),G),nrow=G,ncol=(K+1),byrow=TRUE)
      Sig_hat <- list()
      for(g in 1:G){Sig_hat[[g]] <- rbind(cbind(diag(1,K),0),0)}
    }

    red_sig<-list()
    for(g in 1:G){
      red_sig[[g]] <- Sig_hat[[g]][1:K,1:K]
    }
    iSig_hat<-list()
    for(g in 1:G){
      iSig_hat[[g]] <- MASS::ginv(Sig_hat[[g]])
    }

    par <- list(pi_g=pi_g,xi=xi_hat,m=m_hat,V=V_hat,Vmat=Vmat,mu=mu_hat,Sig=Sig_hat,red_sig=red_sig,iSig=iSig_hat)

    return(list(z=z,par=par))
  }

  ## small EM initialization
  if(method == "small_EM"){
    if(is.na(nrep)&is.na(niter)){S=50;I=10}else{
      if(is.na(nrep)){S=50;I=niter}
      if(is.na(niter)){S=nrep;I=10}
      if(!(is.na(nrep)|is.na(niter))){S=nrep;I=niter}
    }
    # require(foreach)
    # require(doParallel)
    theta <- W/rowSums(W)
    latent <- log(theta/theta[,(K+1)])
    z_store <- list()
    # par_store <- list()
    lik_store <- numeric(S)
    # numCores <- detectCores()
    # registerDoParallel(numCores)
    # foreach(s=1:50) %dopar% {
    for(s in 1:S){
      if(verb){cat("small_EM ",s,"-th start", "\n",sep="")}
      z <- unmap(mapz=sample(1:G,n,replace=TRUE),G=G)
      # random initialization
      pi_g <- colSums(z)/sum(z)
      xi_hat <- rep(1,n)
      m_hat <- latent
      V_hat <- vector("list",n)
      Vmat<-matrix(rep(c(rep(0.1,K),0),n),nrow=n,byrow=T)

      for(i in 1:n){V_hat[[i]]<-diag(Vmat[i,])}

      musig_up <- mst_fun(m=m_hat,V=V_hat,G=G,K=K,z=z)
      if(all(!is.na(musig_up$mu))){
        mu_hat <- musig_up$mu
        Sig_hat <- musig_up$Sig
      }else{
        mu_hat <- matrix(rep(colMeans(m_hat),G),nrow=G,ncol=(K+1),byrow=TRUE)
        Sig_hat <- list()
        for(g in 1:G){Sig_hat[[g]] <- rbind(cbind(diag(1,K),0),0)}
      }

      red_sig<-list()
      for(g in 1:G){
        red_sig[[g]] <- Sig_hat[[g]][1:K,1:K]
      }
      iSig_hat<-list()
      for(g in 1:G){
        iSig_hat[[g]] <- MASS::ginv(Sig_hat[[g]])
      }

      initial_par <- list(pi_g=pi_g,xi=xi_hat,m=m_hat,V=V_hat,Vmat=Vmat,mu=mu_hat,Sig=Sig_hat,red_sig=red_sig,iSig=iSig_hat)

      old_par <- initial_par
      loglik_approx <- NULL
      loglik_approx[1] <- -Inf
      for(it in 2:I){
        zhat <- z_fun(W=W,m=old_par$m,V=old_par$V,Vmat=old_par$Vmat,mu=old_par$mu,Sig=old_par$Sig,red_sig=old_par$red_sig,pi_g=old_par$pi_g,G=G,it=it)

        varpar_up <- varpar_fun(W=W,m=old_par$m,V=old_par$V,mu=old_par$mu,Sig=old_par$Sig,K=K,z=zhat,iSig=old_par$iSig)
        xi_hat <- varpar_up$xi
        m_hat <- varpar_up$m
        V_hat <- varpar_up$V
        Vmat<-varpar_up$Vmat

        # M-step
        pi_g_hat <- colSums(zhat)/sum(zhat)
        musig_up <- mst_fun(m=m_hat,V=V_hat,G=G,K=K,z=zhat)
        mu_hat <- musig_up$mu
        Sig_hat <- musig_up$Sig

        iSig_hat<-list()
        for(g in 1:G){
          iSig_hat[[g]] <- MASS::ginv(Sig_hat[[g]])
        }
        red_sig<-list()
        for(g in 1:G){
          red_sig[[g]] <- Sig_hat[[g]][1:K,1:K]
        }

        # calculate approximated loglikelihood
        loglik_approx[it] <- loglik.approx(W=W,m=m_hat,V=V_hat,Vmat=Vmat,mu=mu_hat,Sig=Sig_hat,red_sig=red_sig,pi_g=pi_g_hat,G=G,it=it)

        old_par <- list(pi_g=pi_g_hat,xi=xi_hat,m=m_hat,V=V_hat,Vmat=Vmat,mu=mu_hat,Sig=Sig_hat,red_sig=red_sig,iSig=iSig_hat)

      }
      z_store[[s]]<-zhat
      # par_store[[s]]<-old_par
      # loglik_approx
      lik_store[s]<-loglik_approx[I]
    }
    maxlik <- which.max(lik_store)
    z_clus <- as.matrix(z_store[[maxlik]])
    pi_g_max <- colSums(z_clus)/sum(z_clus)
    # initial_par <- par_store[[maxlik]]
    initial_par$pi_g <- pi_g_max

    return(list(z=z_clus,par=initial_par))
  }
}
