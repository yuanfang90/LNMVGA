#' Function for the Clustering Algorithm
#'
#' This function is the paralleled running version of the main clustering algorithm
#' @param data Input data here. If sim==TRUE, data should be a list of multiple datasets, with each dataset as a list of counts W and true_lab (true class label). If no true label, set true_lab as NAs. If not simulation, data should be as the same format as described for each dataset of the simulation.
#' @param run Only specify when sim==FALSE. When sim==TRUE, automatically becomes the number of datasets contained in the simulation dataset list.
#' @param Gmax Input the maximum of number of component wants to fit.
#' @param initial Specify method for initializing z_ig. Possible values could be "kmeans", "random", "small_EM". Default is "kmeans".
#' @param runtime Logical variable, if outputting the running time of the whole procedure or not.
#' @param threshold Threshold for the Atiken's stopping creterion for convergence.
#' @param verb Logical variable, if the key steps of the algortihm and approximated loglikelihood for each iteration are printed.
#' @param maxiter Maximum number of iteration. If specified, algorithm will stop by either below the threshold or maxiter reached. If not specified, algorithm will only be monitored by convergence criterion.
#' @param nrep Default is NA. Only needed if "small_EM" is specified for initial. Number of random starts for the small EM initialization.
#' @param niter Default is NA. Only needed if "small_EM" is specified for initial. Number of iterations for each random start in the small EM initialization.
#' @param sim Indicator of whether this is simulated data. Simulated data input must as a list of multiple datasets (indexed by "run"), with each dataset must be a list of W and true_lab. Default is FALSE.
#' @return  A list contains the results for all datasets (runs) and all number of components (from 1 to Gmax). Results include the BIC, ICL, ARI if true labels are not NAs, run time in seconds if runtime==TRUE, as matrices, and the best select G by BIC/ICL for each dataset together with the corresponding ARIs. Results inherited from LNM.clust for each G and each data set are also stored.
#' @keywords clust.parallel.all
#' @examples
#' # generate data using Data.temp <- generate_data(G = 2, num_observation = c(50,50), K = 2, true_mu = list(c(0,1,0),c(-2,-5,0)),true_Sig=list(rbind(cbind(diag(1,2),0),0),rbind(cbind(diag(1,2),0),0)), seed.no = 1234, M = 10000, truelab = TRUE)
#'
#' LNMMM(data=Data.temp,run=1,Gmax=5,initial="kmeans",runtime=TRUE,threshold=1e-4,verb=TRUE,sim=FALSE)

LNMMM <- function(data,run,Gmax,initial="kmeans",runtime=TRUE,threshold,verb,maxiter=NA,nrep=NA,niter=NA,sim=FALSE){
  library("plyr")
  library(mclust)
  library(MASS)


  total_run <- ifelse(sim==FALSE,run,length(data))

  Gmax<-Gmax

  para_mat<-data.frame(run=rep(1:total_run,each=Gmax),G=rep(1:Gmax,total_run))

  if(.Platform$OS.type == "windows") {
    library("doParallel")
    no_cores <- parallel::detectCores() - 4
    registerDoParallel(cores=no_cores)
  }

  if (.Platform$OS.type == "unix") {
    library("doMC")
    registerDoMC(50) #parallel::detectCores() gives the total number of cores available
  }

  parallel_all<-function(pm){
    run_Data<-para_mat[pm,1]
    G_act<-para_mat[pm,2]
    results<-LNM.clust(data,run,G_act,initial,runtime=TRUE,threshold,verb,maxiter=NA,nrep=NA,niter=NA,sim=FALSE)
  }

  output_all <- foreach(pm = 1:nrow(para_mat), .errorhandling = "pass") %dopar% {
    parallel_all(pm)
  }

  nn<-para_mat[which(para_mat[,1]==1),]

  BIC_mat<-matrix(NA,nrow=total_run,ncol=Gmax)
  colnames(BIC_mat)<-paste("G=",nn[,2],sep="")
  rownames(BIC_mat)<-paste("Dataset",1:total_run,sep="=")

  ICL_mat<-matrix(NA,nrow=total_run,ncol=Gmax)
  colnames(ICL_mat)<-paste("G=",nn[,2],sep="")
  rownames(ICL_mat)<-paste("Dataset",1:total_run,sep="=")

  ARI_mat<-matrix(NA,nrow=total_run,ncol=Gmax)
  colnames(ARI_mat)<-paste("G=",nn[,2],sep="")
  rownames(ARI_mat)<-paste("Dataset",1:total_run,sep="=")

  if(runtime == TRUE){
    Time_mat<-matrix(NA,nrow=total_run,ncol=Gmax)
    colnames(Time_mat)<-paste("G=",nn[,2],sep="")
    rownames(Time_mat)<-paste("Dataset",1:total_run,sep="=")
  }

  for (i in 1:total_run){
    each_run<-which(para_mat[,1]==i)
    for (g in 1:Gmax){
      if (length(output_all[[each_run[g]]])>2){
        BIC_mat[i,g]<-output_all[[each_run[g]]]$BIC
        ICL_mat[i,g] <- output_all[[each_run[g]]]$ICL
        ARI_mat[i,g]<-mclust::adjustedRandIndex(mclust::map(output_all[[each_run[g]]]$z),output_all[[each_run[g]]]$true)

        if(runtime == TRUE){
          Time_mat[i,g]<-output_all[[each_run[g]]]$time
        }else{print("No running time recorded. If needed, specify runtim=TRUE. ")}
      }
    }
  }

  BIC_select <- apply(BIC_mat, MARGIN = 1, which.max)
  ICL_select <- apply(ICL_mat, MARGIN = 1, which.max)

  ARI_BIC<-NULL
  for (run in 1:total_run){
    ARI_BIC[run]<-ARI_mat[run,BIC_select[run]]
  }

  ARI_ICL<-NULL
  for (run in 1:total_run){
    ARI_ICL[run]<-ARI_mat[run,ICL_select[run]]
  }

  if(runtime == TRUE&all(!is.nan(ARI_mat))){
    return(list(BIC=BIC_mat,ICL=ICL_mat,all_ARI=ARI_mat,runtime=Time_mat,BIC.select=BIC_select,ICL.select=ICL_select,ARI.BIC=ARI_BIC,ARI.ICL=ARI_ICL,all.results=output_all))
  }
  if(runtime == FALSE&all(!is.nan(ARI_mat))){
    return(list(BIC=BIC_mat,ICL=ICL_mat,all_ARI=ARI_mat,BIC.select=BIC_select,ICL.select=ICL_select,ARI.BIC=ARI_BIC,ARI.ICL=ARI_ICL,all.results=output_all))
  }
  if(runtime == TRUE&any(is.nan(ARI_mat))){
    return(list(BIC=BIC_mat,ICL=ICL_mat,runtime=Time_mat,BIC.select=BIC_select,ICL.select=ICL_select,all.results=output_all))
  }
  if(runtime == FALSE&any(is.nan(ARI_mat))){
    return(list(BIC=BIC_mat,ICL=ICL_mat,BIC.select=BIC_select,ICL.select=ICL_select,all.results=output_all))
  }
}
