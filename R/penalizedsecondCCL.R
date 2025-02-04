#' @title Second order conditional composite likelihood for multivariate LGCPs
#'
#' @description This is a conditional composite likelihood method for semi-parametric estimation of multivariate log Gaussian Cox processes (LGCP).
#'
#' @param X A multivariate log Gaussian Cox process. Must be of class ppp.
#' @param psi0 Optional. Initial value for psi. If psi0 = NULL, then psi0 is simulated from Unif(-0.25,0.25).
#' @param xi0 Optional. Initial value xi. If xi0 = NULL, then xi0 is simulated using the size of the observation window.
#' @param sigma0 Optional. Initial value sigma. If sigma0 = NULL, then sigma0 is simulated from Unif(0.4,0.6).
#' @param phi0 Optional. Initial value phi. If phi0 = NULL, then phi0 is simulated using the size of the observation window.
#' @param Beta Estimated first order parameters. Beta must be a matrix. The number of rows must correspond to the number of covariates.
#' The number of columns must correspond to the number of point types.
#' @param xibound Optional. Parameter space for xi. Default is xibound = (1e-5,Inf).
#' @param sigmabound Optional. Boundary for sigma. Default is sigmabound = (1e-5,Inf).
#' @param phibound Optional. Boundary for phi. Default is phi.bound = (1e-5,Inf).
#' @param lamb Optional. lasso parameter. Default is lamb = 0. Can be a sequence of lambda-values.
#' @param lat Number of common latent fields in the model. Default is lat = 1.
#' @param Rmax Maximal distance between distinct pairs of points. If Rmax = NULL, then Rmax is determined by the observation window.
#' @param covariate Covariates. Covariates must be a matrix. The rows corresponds to the points in the point pattern and the
#' columns indicates the corresponding covariates that are observed at the location of the point pattern.
#' @param ite Maximal number of iterations. Default is ite = 100.
#' @param tol Convergence parameter for likelihood function. The default is tol = 1e-5.
#' @param tol.alpha Convergence parameter for alpha. The default is tol.alpha = 1e-10.
#' @param tol.eta Convergence parameter for eta. The default is tol.eta = 1e-10.
#' @param mu Contraint parameter. The default is mu = 1.
#' @param prints Optional. Default prints = TRUE. Printing value of likelihood function along with convergence rate for each iteration.
#' @importFrom spatstat.geom closepairs
#' @importFrom spatstat.geom is.im
#' @importFrom spatstat.geom is.ppp
#' @importFrom stats runif
#' @return Return estimated parameters of multivariate LGCP X
#' @details To be updated.
#' @author Kristian Bjørn Hessellund, Ganggang Xu, Yongtao Guan and Rasmus Waagepetersen.
#' @references Hessellund, K. B., Xu, G., Guan, Y. and Waagepetersen, R. (2020)
#' Second order semi-parametric inference for multivariate log Gaussian Cox processes.
#' @examples
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
#' ### Second order analysis of Washington D.C. street crimes ###
#'
#' # Globals
#' len.x=markedprocess$window$xrange[2]-markedprocess$window$xrange[1]
#' len.y=markedprocess$window$yrange[2]-markedprocess$window$yrange[1]
#' lw=(len.x+len.y)/2
#' latent = 1
#' Rmax   = 1000
#' rl     = 100
#' r      = seq(0,Rmax,length=rl)
#'
#' # Parameter spaces (upper,lower)
#' xib = c(lw/10,lw/500)
#' sib = c(10,0.01)
#' phb = c(lw/10,lw/1000)
#'
#' notrun = TRUE
#' if(notrun ==FALSE){
#' # Parameter estimation
#' par.est = PenalizedSecondOrderCCL(X=markedprocess,covariate = Z_process,Beta = Betahat$betahat,
#' Rmax = Rmax,xibound=xib,sigmabound = sib,phibound=phb, lat = latent)
#'
#' # Parameter estimates
#' alp.est = matrix(par.est[[1]]$alpha[,,1])
#' xi.est  = par.est[[1]]$xi
#' sig.est = par.est[[1]]$sigma
#' phi.est = par.est[[1]]$phi
#'
#' ## Pair correlation functions (PCF)
#' pcf.est=matrix(NA,nrow=nspecies,ncol=rl)
#' for(i in 1:nspecies){
#'   tmp=0
#'   for(k in 1:latent){tmp = tmp + alp.est[i,k]*alp.est[i,k]*exp(-r/xi.est[k])}
#'   pcf.est[i,]=exp(tmp+sig.est[i]^2*exp(-r/phi.est[i]))
#' }
#'
#' ## Plot of PCFs
#' lims=c(min(pcf.est,0.9),max(pcf.est,1.1))
#' plot(r,rep(1,rl),type="l",ylim=lims,main="PCFs",ylab="", lwd=2)
#' for(j in 1:nspecies){lines(r,pcf.est[j,],col=j+1,lwd=2)}
#' txt.legend=c(as.expression(bquote(g[11])),as.expression(bquote(g[22])),
#' as.expression(bquote(g[33])),as.expression(bquote(g[44])),
#' as.expression(bquote(g[55])),as.expression(bquote(g[66])))
#' legend("topright",legend =  txt.legend,lty=1,col=c(2,3,4,5,6,7), cex=1,lwd=2)
#'
#' # cross-PCFs
#' cross.pcf.est=matrix(NA,nrow=nspecies*(nspecies-1)/2,ncol=rl)
#' l=0
#' for(i in 1:nspecies){
#'  for(j in i:nspecies){
#'    if(i==j){next}
#'    l=l+1
#'    tmp=0
#'    for(k in 1:latent){tmp = tmp + alp.est[i,k]*alp.est[j,k]*exp(-r/xi.est[k])}
#'    cross.pcf.est[l,]=exp(tmp)
#'   }
#' }
#'
#' # Plot of cross-PCFs
#' ylims=c(min(0.9,cross.pcf.est),max(1.1,cross.pcf.est))
#' ltys=c(rep(1,8),rep(2,8))
#' plot(r,rep(1,rl),type="l",main="cross-PCFs",ylab="",ylim=ylims)
#' for(i in 1:(nspecies*(nspecies-1)/2)){lines(r,cross.pcf.est[i,],col=i+1,lty=ltys[i],lwd=2)}
#' txt.legend=c(as.expression(bquote(g[12])),as.expression(bquote(g[13])),
#' as.expression(bquote(g[14])),as.expression(bquote(g[15])),as.expression(bquote(g[16])),
#' as.expression(bquote(g[23])),as.expression(bquote(g[24])),as.expression(bquote(g[25])),
#' as.expression(bquote(g[26])),as.expression(bquote(g[34])),as.expression(bquote(g[35])),
#' as.expression(bquote(g[36])),as.expression(bquote(g[45])),as.expression(bquote(g[46])),
#' as.expression(bquote(g[56])))
#' legend("topright",legend =  txt.legend,lty=c(rep(1,8),rep(2,8)),col=2:17,cex=0.7,bg="white",lwd=2)
#' }
#'
#' @export
#'

PenalizedSecondOrderCCL <- function(X=NULL,psi0=NULL,xi0=NULL,sigma0=NULL,phi0=NULL,Beta,xibound=c(1e-5,Inf),sigmabound=c(1e-5,Inf),phibound=c(1e-5,Inf),lamb=0,lat=1,Rmax=NULL,covariate=NULL,ite=100,tol=1e-5,tol.alpha=1e-10,tol.eta=1e-10,mu=1, prints = T){

  if(is.ppp(X)==F){
    return("X must be a point pattern")
  }

  Output     = list()
  marks      = X$marks
  mval       = levels(marks)
  len.x      = X$window$xrange[2] - X$window$xrange[1]
  len.y      = X$window$yrange[2] - X$window$yrange[1]
  len.win    = (len.x + len.y)/2
  nlamb      = length(lamb)
  nlat       = length(lat)

  p = length(levels(X$marks))

  if(is.null(Rmax)==T){Rmax=len.win/10}
  # Finding all the closepairs in the multivariate point process
  cp   <- closepairs(X,rmax=Rmax,twice = F,what = "ijd")

  # Indixing the closepairs with respect to the ith and jth process
  indi <- cp$i
  indj <- cp$j

  d <- cp$d

  # Finding the marks with respect to the indices
  marksi <- marks[indi]
  marksj <- marks[indj]


  npairs    = length(d)
  if(is.null(covariate)==T){
    Zi <- cbind(rep(1,npairs))
    Zj <- cbind(rep(1,npairs))
  }
  if(is.null(covariate)==F){

    if(is.matrix(covariate)==T){
      Zi = covariate[indi,]
      Zj = covariate[indj,]
    }
    if(is.im(covariate)==T){
      # Evaluating the covariate wrt. the ith and jth process
      Zi <- cbind(1,covariate[X[indi]])
      Zj <- cbind(1,covariate[X[indj]])
    }
  }


  fi <- exp(Zi%*%Beta)
  fj <- exp(Zj%*%Beta)
  fij <- matrix(NA,npairs,p^2)
  for(j in 1:p){
    for(k in 1:p){
      fij[,(j-1)*p+k] <- fi[,j]*fj[,k]
    }
  }

  B=rbind(diag(p-1),rep(-1,p-1))

  for(n in 1:nlat){

    if(lat[n] == 0){if(is.null(psi0)){psi0 = matrix(0,nrow=p-1)}}
    if(lat[n] >= 1){if(is.null(psi0)){
      psi0 = matrix(runif((p-1)*lat[n],min = -0.25,max=0.25),nrow=p-1)
    }
    }
    alpha0=B%*%psi0

    q = dim(psi0)[2]

    if(lat[n] == 0){if(is.null(xi0)){xi0  = log(runif(1,1*len.win/100,2*len.win/100))}}
    if(lat[n] >= 1){if(is.null(xi0)){xi0  = log(runif(lat[n],1*len.win/100,2*len.win/100))}}
    if(is.null(sigma0)){sigma0 = log(runif(p,0.4,0.6))}
    if(is.null(phi0)){phi0 = log(runif(p,1*len.win/100,2*len.win/100))}

    alpha.final = array(alpha0,c(p,q,nlamb+1))
    xi.final  = matrix(xi0,nrow=(nlamb+1),ncol=q,byrow = T)
    sig.final = matrix(sigma0,nrow=(nlamb+1),ncol=p,byrow = T)
    phi.final = matrix(phi0,nrow=(nlamb+1),ncol=p,byrow = T)

    time   = c()
    iter   = c()
    lambda = c()
    convr  = c()
    result = NULL

    for(k in 1:nlamb){
      if(prints==T){
        print.txt=paste0("latent: ",lat[n],", lambda: ",lamb[k])
        print(print.txt)
      }

      psi    = psi0
      xis    = xi0
      sigmas = sigma0
      phis   = phi0

      lik=c()
      start=Sys.time()
      lik[1]=second_order_likelihood(psi = psi,xis = xis, sigmas = sigmas, phis=phis,R = d,fi = fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=T)
      for(i in 1:ite){
        if(lat[n] >= 1){

          alphas=B%*%psi

          # alpha estimation
          est     = gradient.alpha1(alphas = alphas,xis = xis, sigmas = sigmas, phis=phis,R = d,fi = fi,fj=fj,mi=marksi,mj=marksj,mval=mval)
          sen.alp = sensitivity.alpha1(alphas = alphas,xis = xis, sigmas = sigmas, phis=phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval)

          if(lamb[k]==0){
            est1=NULL
            for(i in 1:q){est1=rbind(est1,t(B)%*%est[(1+p*(i-1)):(p*i)])}
            est=est1

            sen1=matrix(NA,nrow=(p*q-q),ncol=(p*q-q))
            for(i1 in 1:q){
              for(i2 in 1:q){
                sen1[(1+(p-1)*(i1-1)):((p-1)*i1),(1+(p-1)*(i2-1)):((p-1)*i2)]=t(B)%*%sen.alp[(1+p*(i1-1)):(p*i1),(1+p*(i2-1)):(p*i2)]%*%B
              }
            }
            sen.alp=sen1
            dim.a = p*q-q

            S       = eigen(sen.alp)
            # If sensitivity not PSD
            if(min(S$values)<=0){S$values = S$values + 1e-8}

            U       = S$vectors
            D1      = diag(S$values^(1/2),ncol=dim.a,nrow=dim.a)
            D2      = diag(1/S$values,ncol=dim.a,nrow=dim.a)
            X       = U%*%D1%*%t(U)
            inv.sen = U%*%D2%*%t(U)
            alps=psi

            Y       = X%*%(-inv.sen%*%est + c(alps))

            XX=eigen(t(X)%*%X)
            if(min(XX$values)<=0){
              print(c("not PS alpha"))
              XX$values = XX$values + 1e-8
            }
            inv.XX = XX$vectors%*%diag(1/XX$values,ncol=dim.a,nrow=dim.a)%*%t(XX$vectors)
            update=matrix(inv.XX%*%t(X)%*%Y,nrow = p-1,byrow=F)

            stp=1
            lik.old = second_order_likelihood(psi = psi,xis = xis, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=T)#+lamb[k]*sum(abs(alphas))
            for(l in 1:10){
              if(l==10){psi.new  = psi; break}
              psi.new = psi+stp*(update-psi)
              lik.new = second_order_likelihood(psi = psi.new,xis = xis, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=T)#+lamb[k]*sum(abs(alp.new))
              if(is.nan(lik.new)){stp = stp/2; next}
              if(lik.new <= lik.old){break}
              else{stp = stp/2}
            }

            alp.new=matrix(B%*%psi.new,nrow=p)
          }

          if(lamb[k]!=0){
            S       = eigen(sen.alp)
            # If sensitivity not PSD
            if(min(S$values)<=0){S$values = S$values + 1e-8}
            dim.a = p*q

            U       = S$vectors
            D1      = diag(S$values^(1/2),ncol=dim.a,nrow=dim.a)
            D2      = diag(1/S$values,ncol=dim.a,nrow=dim.a)
            X       = U%*%D1%*%t(U)
            inv.sen = U%*%D2%*%t(U)

            Y       = X%*%(-inv.sen%*%est + c(alphas))

            soft.tres=function(z,gm){
              if(z>0 & gm < abs(z)){val=z-gm}
              if(z<0 & gm < abs(z)){val=z+gm}
              if(gm >= abs(z)){val=0}
              val
            }

            update=update.old=c(alphas)
            C=matrix(0,nrow=p*q,ncol=q)
            for(i in 1:q){
              h=((i-1)*p+1):(i*p)
              C[h,i]=1
            }

            eta=numeric(q)
            for(t in 1:10000){
              l=0
              for(m1 in 1:p){
                for(m2 in 1:q){
                  l = l+1
                  tmp=(1/p)*t(X[,l])%*%(Y-as.matrix(X[,-l])%*%update[-l])-mu/2*(update[-l]%*%(C%*%C[l,])[-l,] +t(C[l,])%*%eta)
                  update[l]=soft.tres(c(tmp),lamb[k]/(2))/(t(X[,l])%*%(X[,l])/p+mu*t(C[l,])%*%C[l,])
                }
              }
              eta = eta + t(C)%*%update
              if(sum(abs(update.old-update))<=tol.alpha & sum(abs(t(C)%*%update)) <=tol.eta){break}
              update.old=update
            }
            update=matrix(update,nrow=p,byrow=F)

            stp=1
            lik.old = second_order_likelihood(psi = alphas,xis = xis, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=F)+lamb[k]*sum(abs(alphas))
            for(l in 1:10){
              if(l==10){alp.new  = alphas; break}
              alp.new = alphas+stp*(update-alphas)
              lik.new = second_order_likelihood(psi = alp.new,xis = xis, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=F)+lamb[k]*sum(abs(alp.new))
              if(is.nan(lik.new)){stp = stp/2; next}
              if(lik.new <= lik.old){break}
              else{stp = stp/2}
            }
            psi.new = alp.new[1:(p-1),]
          }

          sen.xi   = sensitivity.xi(alphas = alp.new,xis = xis, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval)
          est      = gradient.xi(alphas = alp.new,xis = xis, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval)
          S        = eigen(sen.xi)

          if(min(S$values)<=0){S$values = S$values + 1e-8}

          U       = S$vectors
          D1      = diag(S$values^(1/2),ncol=q,nrow=q)
          D2      = diag(1/S$values,ncol=q,nrow=q)
          X       = U%*%D1%*%t(U)
          inv.sen = U%*%D2%*%t(U)
          Y       = X%*%(-inv.sen%*%est + c(t(xis)))

          XX=eigen(t(X)%*%X)
          if(min(XX$values)<=0){XX$values = XX$values + 1e-8}

          inv.XX  = XX$vectors%*%diag(1/XX$values,ncol=q,nrow=q)%*%t(XX$vectors)
          update = inv.XX%*%t(X)%*%Y

          lik.old  = second_order_likelihood(psi =  alp.new,xis = xis, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=F)
          stp=1
          for(l in 1:10){
            if(l==10){xis.new  = xis; break}
            xis.new  = xis + stp*(update-xis)
            lik.new = second_order_likelihood(psi = alp.new,xis = xis.new, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=F)
            if(is.nan(lik.new)){stp = stp/2; next}
            if(lik.new <= lik.old){break}
            else{stp = stp/2}
          }
          xis.new[xis.new>=log(xibound[1])]=log(xibound[1])
          xis.new[xis.new<=log(xibound[2])]=log(xibound[2])
        }

        else{
          psi.new = psi =  matrix(rep(0,p-1),nrow=p-1)
          alp.new = alphas = matrix(rep(0,p),nrow=p)
          xis.new = xis
        }

        sen.sig  = sensitivity.sigma2(alphas = alp.new,xis = xis.new, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval)
        est      = gradient.sigma2(alphas = alp.new,xis = xis.new, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval)
        S        = eigen(sen.sig)
        if(min(S$values)<=0){S$values = S$values + 1e-8}

        U       = S$vectors
        D1      = diag(S$values^(1/2),ncol=p,nrow=p)
        D2      = diag(1/S$values,ncol=p,nrow=p)
        X       = U%*%D1%*%t(U)
        inv.sen = U%*%D2%*%t(U)
        Y       = X%*%(-inv.sen%*%est + c(t(sigmas)))

        XX=eigen(t(X)%*%X)
        if(min(XX$values)<=0){XX$values = XX$values + 1e-8}

        inv.XX  = XX$vectors%*%diag(1/XX$values,ncol=p,nrow=p)%*%t(XX$vectors)
        update = inv.XX%*%t(X)%*%Y

        lik.old = second_order_likelihood(psi = alp.new,xis = xis.new, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=F)
        stp=1
        for(l in 1:10){
          if(l==10){sig.new  = sigmas; break}
          sig.new  = sigmas + stp*(update-sigmas)
          lik.new = second_order_likelihood(psi = alp.new,xis = xis.new, sigmas = sig.new,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=F)
          if(is.nan(lik.new)){stp = stp/2; next}
          if(lik.new <= lik.old){break}
          else{stp = stp/2}
        }
        sig.new[sig.new>=log(sigmabound[1])]=log(sigmabound[1])
        sig.new[sig.new<=log(sigmabound[2])]=log(sigmabound[2])

        sen.phi  = sensitivity.phis(alphas = alp.new,xis = xis.new, sigmas = sig.new,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval)
        est      = gradient.phis(alphas = alp.new,xis = xis.new, sigmas = sig.new,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval)
        S        = eigen(sen.phi)
        if(prod(S$values)<=0){S$values = S$values + 1e-8}

        U       = S$vectors
        D1      = diag(S$values^(1/2),ncol=p,nrow=p)
        D2      = diag(1/S$values,ncol=p,nrow=p)
        X       = U%*%D1%*%t(U)
        inv.sen = U%*%D2%*%t(U)
        Y       = X%*%(-inv.sen%*%est + c(t(phis)))

        XX=eigen(t(X)%*%X)
        if(min(XX$values)<=0){
          print(c("not PD phi"))
          XX$values = XX$values + 1e-8
        }

        inv.XX  = XX$vectors%*%diag(1/XX$values,ncol=p,nrow=p)%*%t(XX$vectors)
        update = inv.XX%*%t(X)%*%Y

        stp=1
        lik.old = second_order_likelihood(psi = alp.new,xis = xis.new, sigmas = sig.new,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=F)
        for(l in 1:10){
          if(l==10){phi.new  = phis;break}
          phi.new  = phis + stp*(update-phis)
          lik.new = second_order_likelihood(psi = alp.new,xis = xis.new, sigmas = sig.new,phis = phi.new,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=F)
          if(is.nan(lik.new)){stp = stp/2;next}
          if(lik.new <= lik.old){break}
          else{stp = stp/2}
        }
        phi.new[phi.new>=log(phibound[1])]=log(phibound[1])
        phi.new[phi.new<=log(phibound[2])]=log(phibound[2])

        lik.old=second_order_likelihood(psi = alphas,xis = xis, sigmas = sigmas,phis = phis,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=F)+lamb[k]*sum(abs(alphas))
        lik.new=second_order_likelihood(psi = alp.new,xis = xis.new, sigmas = sig.new,phis = phi.new,R=d,fi=fi,fj=fj,mi=marksi,mj=marksj,mval=mval,stz=F)+lamb[k]*sum(abs(alp.new))

        if(prints==T){
          print.txt=paste0("likelihood: ",round(lik.new,digits = 1),", convergence: ",round(abs((lik.new-lik.old)/lik.old),digits = 6))
          print(print.txt)
        }
        if(abs((lik.new-lik.old)/lik.old)<tol){break}
        psi    = psi.new
        xis    = xis.new
        sigmas = sig.new
        phis   = phi.new
        if(i == (ite-1)){break}
      }

      time[k]   = Sys.time() - start
      iter[k]    = i
      lambda[k] = lamb[k]
      if(i == (ite-1)){convr[k] = 1}
      if(i < (ite-1)){convr[k] = 0}

      alpha.final[,,k] = alp.new
      xi.final[k,]     = exp(xis.new)
      sig.final[k,]    = sqrt(exp(sig.new))
      phi.final[k,]    = exp(phi.new)

      # Warm start
      psi0 = matrix(psi.new,nrow=p-1)
      if(lat[n] >= 1){for(l in 1:q){if(sum(psi0[,l])==0){
        tmp=runif(p-1,min = -0.25,max=0.25)
        psi0[,l]=tmp}
      }}
      xi0    = xis.new
      sigma0 = sig.new
      phi0   = phi.new
    }

    result$alpha  = array(alpha.final[,,1:nlamb],c(p,q,nlamb))
    result$xis    = matrix(xi.final[1:nlamb,],nrow=nlamb)
    result$sigmas = matrix(sig.final[1:nlamb,],nrow=nlamb)
    result$phis   = matrix(phi.final[1:nlamb,],nrow=nlamb)
    result$lambda = lambda
    result$conv   = convr
    result$time   = time
    result$ite    = iter
    result$lik    = lik.new
    Output[[n]]   = result

    # Reset initial parameters
    psi0 = xi0= sigma0 = phi0 = NULL
  }
  return(Output)
}
