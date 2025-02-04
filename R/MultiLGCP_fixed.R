library(Multilogreg)
sim_lgcp_multi_fixed <- function(basecov,covariate,betas,alphas,xis,sigmas,phis, n.window,n.points,beta0s=NULL){
  
  nspecies <- dim(betas)[1]
  dimyx <- ncol(basecov)
  n.latent <- dim(alphas)[2]
  ###find the intercepts such that there are n.points[i] for each species
  logintensity.case <- matrix(0,nspecies,dimyx^2)
  
  ###If all intercepts are given, then use given beta0s; otherwise calculate beta0s so that
  ### average number of points are n.points[i] for process X_i
  
  if(is.null(beta0s)){
    beta0s <- numeric()
    for(i in 1:nspecies){
      logint <- basecov
      for(j in seq_len(ncol(betas))) {
        beta0s[i] <- log(n.points[i]/area.owin(n.window))-mean(c(as.matrix(covariate[[j]]*betas[i,j]+logint)))
        logint <- logint + covariate[[j]]*betas[i,j]
        
      }
      
      logintensity.case[i,] <- c(t(as.matrix(beta0s[i]+logint)))
    }
  }else{
    for(i in 1:nspecies){
      for(j in seq_len(ncol(betas))) {
        logint <- logint + covariate[[j]]*betas[i,j]
      }
      logintensity.case[i,] <- c(t(as.matrix(beta0s[i]+ logint +log(n.points[i]))))
    }
  }
  
  ###Create sqaure window for point processes
  gridcord <- as.data.frame(gridcenters(window = n.window, nx = dim(basecov)[1], ny = dim(basecov)[2]))

  ###The uniq x and y seq
  xx <- unique(gridcord$x)
  yy <- unique(gridcord$y)
  
  ###Generate a Gaussian Random Field for the common latent intensity
  for (j in 1:n.latent){
    #Y1 <- RFsimulate(RMexp(var=1,scale=scales_c[j]), x=xx, y=yy, grid=TRUE)
    Y1 <- RFsimulate(RMexp(var=1,scale=xis[j]), x=xx, y=yy, grid=TRUE)
    for (l in 1:nspecies)
      logintensity.case[l,]=logintensity.case[l,]+alphas[l,j]*c(t(as.matrix(Y1)))
  }
  
  ###Generate a Gaussian Random Field for the individual latent intensity
  for (l in 1:nspecies){
    #U = RFsimulate(RMexp(var=1,scale=scales_res[l]), x=xx, y=yy, grid=TRUE)
    U = RFsimulate(RMexp(var=1,scale=phis[l]), x=xx, y=yy, grid=TRUE)
    logintensity.case[l,]=logintensity.case[l,]+sigmas[l]*c(t(as.matrix(U)))
  }
  
  mu <- process <- list()
  markedprocess <- NULL
  for (l in 1:nspecies){
    ###The intensity of LGCP
    intensity_lgcp <- exp(logintensity.case[l,]-sum(alphas[l,]^2)/2-sigmas[l]^2/2)
    intensity.df <- data.frame("x"= gridcord$x, "y" = gridcord$y, "value" = c(intensity_lgcp))
    mu[[l]] <- as.im(intensity.df)
    
    ###Simulate inhomogeneous poisson process from lgcp intensity
    process[[l]] <- rpoispp(lambda = mu[[l]])
    
  }
  
  
  ##Exact the pixel image for covariates
  
  pix <- as.im(data.frame("x"= gridcord$x, "y" = gridcord$y, "value" = c(t(covariate))))
  
  #gridcord <- data.frame(x=gridcord$x[isin],y=gridcord$y[isin])
  
  ###Aggregate all points into one marked point process
  win <- process[[1]]$window
  tmp <- NULL
  for (l in 1:nspecies){
    tmp <-rbind(tmp,cbind(process[[l]]$x,process[[l]]$y,l))
  }
  
  isin <- inside.owin(x = tmp[,1], y = tmp[,2], w = n.window)
  markedprocess <- ppp(x=tmp[isin,1],y=tmp[isin,2],marks=as.factor(tmp[isin,3]),window=win)
  
  
  list(process = process, mu = mu, pix = pix, beta0s=beta0s, markedprocess=markedprocess)
}

