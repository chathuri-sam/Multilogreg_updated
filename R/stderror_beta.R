#' @title Estimation of variance for first order parameters
#'
#' @description This is a method to estimate the variance of the first order parameters
#'
#' @param X A multivariate point process. X must be of class ppp.
#' @param covariate Covariates. The covariates must be a matrix. The rows corresponds to the points in the point pattern and the
#' columns indicates the corresponding covariates that are observed at the location of the point pattern.
#' @param Beta Estimated first order parameters. Beta must be a matrix. The number of rows must correspond to the number of covariates.
#' The number of columns must correspond to the number of point types.
#' @param r A vector of distances
#' @param bwd The bandwidth used
#' @param kern The kernel function. Default is Epanechnikov. Alternatively, the kernel function can be Indicator.
#' @param cPCF Non parametric estimates of cross PCF ratios. Default is cPCF = NULL. If cPCF = NULL, then the cross PCF ratios are estimated
#' using the function CrossPCF
#' @param method Method to estimate the variance of the first order estimating function. Default method = "single", alternatively method = "pool".
#' If method = "single", the variance is estimated based on single estimates of cross PCF ratios. If method = "pool", the variance is estimated
#' based on an estimate of the pooled cross PCF.
#' @param cut For refined cross PCF ratio estimator. Default cut = NULL. If cut = NULL, the parameter R^* for the constraint is chosen
#' manually using the function ChooseRange
#' @importFrom spatstat.geom closepairs
#' @importFrom spatstat.geom crosspairs
#' @importFrom stats optim
#' @return Return estimated cross PCF ratios, r distances used, bandwidth used, refined cross PCF ratios and R^*
#'
#' @author Kristian Bjørn Hessellund, Ganggang Xu, Yongtao Guan and Rasmus Waagepetersen.
#' @references Hessellund, K. B., Xu, G., Guan, Y. and Waagepetersen, R. (2020)
#' Second order semi-parametric inference for multivariate log Gaussian Cox processes.
#' @examples
#'
#' ### First order analysis of Washington D.C. street crimes ###
#'
#' nspecies <- length(unique(markedprocess$marks))
#' ncovs <- dim(Z_process)[2]
#'
#' # Initial parameters for street crime data
#' Beta0=matrix(runif(nspecies*ncovs,-0.5,0.5),nrow =nspecies)
#'
#' # Choosing the last type of street crime as the baseline
#' Beta0 <- as.matrix(t(scale(Beta0,center=Beta0[nspecies,],scale=FALSE)))
#'
#' # Parameter estimation
#' Betahat=FirstOrderCCL(X=markedprocess,Beta0 = Beta0,covariate = Z_process)
#'
#' # Parameter estimates
#' colnames(Betahat$betahat)=c("Burglary", "Assault w. weapon", "Motor v. theft",
#' "Theft f. auto", "Robbery", "Other theft (baseline)")
#' rownames(Betahat$betahat)=c("Intercept", "%African", "%Hispanic", "%Males",
#'  "Median income", "%Household","%Bachelor", "Distance to police station")
#' Betahat
#'
#' ### Non-parametric estimates of cross PCF ratios ###
#'
#' rseq = seq(0,3000,length=100)
#' cpcf=CrossPCF(X = markedprocess,covariate = Z_process,Beta = Betahat$betahat,r = rseq
#' ,bwd = 200,cut = NULL)
#'
#' plot(rseq,rep(1,100),type="l",col=1,ylim=c(0,2),xlab="r",ylab="",main="PCFs")
#' for(i in 1:(nspecies-1)){
#'   lines(rseq,cpcf$cPCF.m[i,i,],,col=i)
#' }
#'
#' plot(rseq,rep(1,100),type="l",col=1,ylim=c(0,2),xlab="r",ylab="",main="PCFs")
#' for(i in 1:(nspecies-1)){
#'   for(j in 1:(nspecies-1)){
#'     if(i==j){next}
#'     lines(rseq,cpcf$cPCF.m[i,j,],,col=1)
#'   }
#' }
#'
#' ### Standard errors for Beta ###
#'
#' var.est=Firstorder_var(X=markedprocess,covariate = Z_process,Beta=Betahat$betahat,
#' r = seq(0,3000,length=100),bwd = 200,cPCF = cpcf$cPCF.m,method = "pool")
#'
#' std.error=matrix(sqrt(diag(var.est$COV)),nrow=ncovs)
#' colnames(std.error)=c("Burglary", "Assault w. weapon", "Motor v. theft","Theft f. auto", "Robbery")
#' rownames(std.error)=c("Intercept", "%African", "%Hispanic", "%Males",
#'  "Median income", "%Household","%Bachelor", "Distance to police station")
#' std.error
#'
#' @export
#'

Firstorder_var <- function(X,covariate,Beta,r,bwd,kern="Epanechnikov",cPCF=NULL,method="single",cut=1){

  #########Find the sensitivity matrix using pooled process
  # p <- length(pix)+1
  # nspecies <- length(unique(markedprocess$marks))
  # #Extract covariate information
  # tmp <- NULL
  # for(j in 1:(p-1))
  #   tmp <- cbind(tmp,(pix[[j]])[markedprocess])
  # Z <-   cbind(1,tmp)

  p = dim(covariate)[2]
  Z = covariate

  MU <- exp(Z%*%Beta)
  Prob <- MU/rowSums(MU)
  S <- matrix(NA,p*(nspecies-1),p*(nspecies-1))
  for(j in 1:(nspecies-1)){
    probj <- Prob[,j]
    for(k in j:(nspecies-1)){
      probk <- Prob[,k]
      if(j==k)
        tmp <- t(Z*probj*(1-probk))%*%Z else
          tmp <- t(Z*probj*(-probk))%*%Z
        S[1:p+(j-1)*p,1:p+(k-1)*p] <- tmp
        S[1:p+(k-1)*p,1:p+(j-1)*p] <- t(tmp)
    }
  }

  ####Find the variance matrix using pooled process
  ###Step 1: find the cPCF
  if(is.null(cPCF))  cPCF <- CrossPCF(X,covariate,Beta,r,bwd,kern)$cPCF

  n.r <- length(r)
  R <- max(r)
  nspecies <- dim(Beta)[2]
  marks <- X$marks
  mval <- levels(marks)

  V1 <- matrix(0,p*(nspecies-1),p*(nspecies-1))
  Tpct <-  cPCF[,,-1]*0


  for(i1 in 1:nspecies){
    for(i2 in 1:nspecies){

      if(i1==i2){
        ind <- marks==mval[i1]
        Z_sub <- Z[ind,,drop=FALSE]
        subprocess <- X[ind]
        cpairs <- closepairs(subprocess,rmax=R,twice = TRUE,what = "ijd")
        d <- cpairs$d
        indi <- cpairs$i
        indj <- cpairs$j
        Zi <- Z_sub[indi,,drop=FALSE]
        Zj <- Z_sub[indj,,drop=FALSE]
      }else{
        ind1 <- marks==mval[i1]
        ind2 <- marks==mval[i2]
        Z_sub1 <- Z[ind1,,drop=FALSE]
        Z_sub2 <- Z[ind2,,drop=FALSE]
        subprocess1 <- X[ind1]
        subprocess2 <- X[ind2]
        cpairs <- crosspairs(subprocess1,subprocess2,rmax=R,what = "ijd")
        d <- cpairs$d
        indi <- cpairs$i
        indj <- cpairs$j
        Zi <- Z_sub1[indi,,drop=FALSE]
        Zj <- Z_sub2[indj,,drop=FALSE]
      }

      MUi <- exp(Zi%*%Beta)
      MUj <- exp(Zj%*%Beta)
      Probi <- MUi/rowSums(MUi)
      Probj <- MUj/rowSums(MUj)

      ###pooled PCF function (no longer isotropic)
      gxstar <- 0
      Gtmp <- (cPCF[,,-1])/2+(cPCF[,,-n.r])/2
      idx <- .bincode(d, breaks=r, right = TRUE, include.lowest = TRUE)

      for(j in 1:nspecies){
        for(k in 1:nspecies){
          gtmp <- (Gtmp[j,k,])[idx]
          gxstar <- gxstar + gtmp*Probi[,j]*Probj[,k]
        }
      }
      ###My method
      if(method=="pool")
        wstar <- 1/gxstar else if(method=="single"){
          ####Kristian's method
          wstar <- 1/((Gtmp[i1,i2,])[idx]*Probi[,i1]*Probj[,i2])/(nspecies^2)
        }

      for(j in 1:(nspecies-1)){
        for(k in j:(nspecies-1)){
          gtmp <- Gtmp[j,k,]
          w <- gtmp[idx]
          for(l in 1:nspecies){
            gtmp1 <- Gtmp[j,l,]
            gtmp2 <- Gtmp[k,l,]
            w <- w - Probj[,l]*gtmp1[idx]-Probi[,l]*gtmp2[idx]
          }
          w <- (w+gxstar)*Probi[,j]*Probj[,k]*wstar
          ind <- w<0&(d>(R*cut))
          w[ind] <- 0
          V1[1:p+(j-1)*p,1:p+(k-1)*p] <- V1[1:p+(j-1)*p,1:p+(k-1)*p]+ t(Zi*w)%*%Zj
        }
      }

    }
  }

  ###Symmetrize V1
  for(j in 1:(nspecies-1)){
    for(k in j:(nspecies-1)){
      V1[1:p+(k-1)*p,1:p+(j-1)*p] <- t(V1[1:p+(j-1)*p,1:p+(k-1)*p])
      Tpct[k,j,] <- Tpct[j,k,]
    }
  }

  V <- V1+S
  COV <- solve(S)%*%V%*%solve(S)
  list(S=S,V=V,COV=COV,cPCF=cPCF,V1=V1)
}
