#' @title Non-parametric estimation of cross PCF ratios
#'
#' @description This is a method to estimate cross PCF ratios non-parametrically
#'
#' @param X A multivariate point process. X must be of class ppp.
#' @param covariate Covariates. The covariates must be a matrix. The rows corresponds to the points in the point pattern and the
#' columns indicates the corresponding covariates that are observed at the location of the point pattern.
#' @param Beta Estimated first order parameters. Beta must be a matrix. The number of rows must correspond to the number of covariates.
#' The number of columns must correspond to the number of point types.
#' @param r A vector of distances
#' @param bwd The bandwidth used
#' @param kern The kernel function. Default is Epanechnikov. Alternatively, the kernel function can be Indicator.
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
#' @export
#'
CrossPCF <- function(X,covariate,Beta,r,bwd,kern="Epanechnikov",cut=NULL){
  n.r <- length(r)
  R <- max(r)
  nspecies <- dim(Beta)[2]
  marks <- X$marks
  mval <- levels(marks)
  p <- dim(covariate)[2]
  Z=covariate
  #Extract covariate information
  #tmp <- NULL
  #for(j in 1:(p-1))
  #  tmp <- cbind(tmp,(pix[[j]])[markedprocess])
  #Z <-   cbind(1,tmp)

  cPCF <- array(NA,dim=c(nspecies,nspecies,n.r))

  if(kern=="Indicator"){
    Kern_IND <- function(x){
      (abs(x)<=1)/2
    }
    Kern <- Kern_IND
  }else if(kern=="Epanechnikov"){
    Kern_EP <- function(x){
      (1-x^2)*3/4*(abs(x)<=1)
    }
    Kern <- Kern_EP
  }
  ###Find the cross-PCF ratios relative to the pooled process
  ###Find the cross-PCF ratios relative to the pooled process
  for(j in 1:nspecies){
    for(k in j:nspecies){

      if(j==k){
        ind <- marks==mval[j]
        Z_sub <- Z[ind,,drop=FALSE]
        subprocess <- X[ind]
        cpairs <- closepairs(subprocess,rmax=R,twice = FALSE,what = "ijd")
        d <- cpairs$d
        indi <- cpairs$i
        indj <- cpairs$j
        Zi <- Z_sub[indi,]
        Zj <- Z_sub[indj,]
      }else{
        ind1 <- marks==mval[j]
        ind2 <- marks==mval[k]
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
      MUi <- MUi/rowSums(MUi) ##USE probability for cross-PCFs
      MUj <- MUj/rowSums(MUj)


      for(i in 1:n.r){
        ind <- (d>r[i]-bwd)&(d<=r[i]+bwd)
        n.sub <- sum(ind)

        if(n.sub!=0){
          MUi_sub <- MUi[ind,,drop=FALSE]
          MUj_sub <- MUj[ind,,drop=FALSE]
          d_sub <- d[ind]
          u_sub <- (d_sub-r[i])/bwd

          ##Find double the pairs if j==k
          if(j==k) ratio <- 2 else ratio <- 1
          cPCF[k,j,i]<- cPCF[j,k,i] <- sum(Kern(u_sub)/(MUi_sub[,j]*MUj_sub[,k]))*ratio
        }
      }
    }
  }

  ###Use the last one as the base line
  for(j in 1:nspecies){
    for(k in j:nspecies){
      cPCF[k,j,]<- cPCF[j,k,] <-  cPCF[k,j,]/cPCF[nspecies,nspecies,]
    }
  }

  #######Constraint modification

  if(is.null(cut))  cut <- ChooseRange(X,covariate,Beta,r,bwd,kern,cPCF=cPCF)$R/R


  refined_cPCF <- function(cPCF,r,cut=0){
    n.s <- dim(cPCF)[1]
    n.r <- length(r)
    cPCF.m <- cPCF
    R <- max(r)
    idx <- min(which(r>cut*R))
    for(i in idx:n.r){
      Yhat <- cPCF.m[,,i]
      Ydiag <- diag(1/sqrt(diag(Yhat)))
      Ystd <- Ydiag%*%Yhat%*%Ydiag
      if(sum(Ystd>1)>0){
        ft <- function(par){
          d <- diag(c(par[1:(n.s-1)],1))
          b= matrix(0, n.s, n.s)
          b[upper.tri(b, diag=FALSE)]=par[-(1:(n.s-1))]
          b <- b+t(b)+diag(n.s)
          Y <- d%*%b%*%d
          -mean((Y-Yhat)^2)
        }

        par <- optim(par =rep(1,(n.s-1)*(n.s+2)/2), fn = ft, control = list(fnscale = -1),method ="L-BFGS-B",lower = c(rep(10^(-10),n.s-1),rep(0,n.s*(n.s-1)/2)),upper=c(rep(Inf,n.s-1),rep(1,n.s*(n.s-1)/2)))$par

        d <- diag(c(par[1:(n.s-1)],1))
        b= matrix(0, n.s, n.s)
        b[upper.tri(b, diag=FALSE)]=par[-(1:(n.s-1))]
        b <- b+t(b)+diag(n.s)
        cPCF.m[,,i] <- d%*%b%*%d
      }
    }

    cPCF.m
  }

  cPCF.m <- refined_cPCF(cPCF,r,cut)

  list(cPCF=cPCF,r=r,bwd=bwd,cPCF.m=cPCF.m,cut=cut)
}

