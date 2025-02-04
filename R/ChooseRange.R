
#' @title Choose range for the refined cross PCF ratio estimator
#'
#' @description This is a method to choose the range for the refined cross PCF ratio estimator
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
#' @importFrom spatstat.geom closepairs
#' @importFrom spatstat.geom crosspairs
#' @importFrom stats aggregate
#' @return Return estimated cross PCF ratios, r distances used, bandwidth used, refined cross PCF ratios and R^*
#'
#' @author Kristian Bjørn Hessellund, Ganggang Xu, Yongtao Guan and Rasmus Waagepetersen.
#' @references Hessellund, K. B., Xu, G., Guan, Y. and Waagepetersen, R. (2020)
#' Second order semi-parametric inference for multivariate log Gaussian Cox processes.
#'
#' @export
#'

ChooseRange <- function(X,covariate,Beta,r,bwd,kern="Epanechnikov",cPCF=NULL){

  #########Find the sensitivity matrix using pooled process
  # p <- length(pix)+1
  # #Extract covariate information
  # tmp <- NULL
  # for(j in 1:(p-1))
  #   tmp <- cbind(tmp,(pix[[j]])[markedprocess])
  # Z <-   cbind(1,tmp)

  p = dim(covariate)[2]
  Z = covariate

  ####Find the variance matrix using pooled process
  ###Step 1: find the cPCF
  if(is.null(cPCF))  cPCF <- CrossPCF(X,covariate,Beta,r,bwd,kern)$cPCF

  n.r <- length(r)
  R <- max(r)
  nspecies <- dim(Beta)[2]
  marks <- X$marks
  mval <- levels(marks)

  Count <- Tpct <-  Tpart <-  matrix(0,nspecies,n.r-1)

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
      #idx <- .bincode(d, breaks=r, right = TRUE, include.lowest = FALSE)
      idx <- cut(d,breaks=r,include.lowest=TRUE)
      for(j in 1:nspecies){
        for(k in 1:nspecies){
          gtmp <- (Gtmp[j,k,])[idx]
          gxstar <- gxstar + gtmp*Probi[,j]*Probj[,k]
        }
      }

      for(j in 1:(nspecies)){
        gtmp <- Gtmp[j,j,]
        w <- gtmp[idx]
        for(l in 1:nspecies){
          gtmp1 <- Gtmp[j,l,]
          gtmp2 <- Gtmp[j,l,]
          w <- w - Probj[,l]*gtmp1[idx]-Probi[,l]*gtmp2[idx]
        }
        w <- (w+gxstar)*Probi[,j]*Probj[,j]/gxstar
        datatmp <- data.frame(y=(w<0),y1=w*0+1,y2=w,idx=idx)
        tmp <- aggregate(y~idx,data=datatmp,FUN="sum",drop=FALSE)$y
        ind <- is.na(tmp)
        tmp[ind] <- 0
        Tpct[j,] <- Tpct[j,]+ tmp
        tmp <- aggregate(y1~idx,data=datatmp,FUN="sum",drop=FALSE)$y1
        ind <- is.na(tmp)
        tmp[ind] <- 0
        Count[j,] <- Count[j,]+ tmp

        tmp <- aggregate(y2~idx,data=datatmp,FUN="sum",drop=FALSE)$y2
        ind <- is.na(tmp)
        tmp[ind] <- 0
        Tpart[j,] <- Tpart[j,]+ tmp
      }
    }
  }

  PCT <- Tpct/(Count+0.1)
  Rmax <- numeric()
  for(j in 1:(nspecies)){
    if(sum(PCT[j,]>0.05)!=0)
      Rmax[j] <- min(which(PCT[j,]>0.05)) else Rmax[j] <- n.r
  }

  R <- r[min(Rmax)]

  list(PCT=PCT,r=r,R=R,Tpart=Tpart)
}
