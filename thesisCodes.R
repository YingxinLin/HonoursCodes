
bic<-function(loglik,n,p){
  return(-2*loglik+p*log(n))
}


aic<-function(loglik,p){
  return(-2*loglik+2*p)
}

icl_bic<-function(loglik,postprob,n,p){
  postprob<-postprob[postprob>0]
  EN=-sum(postprob*log(postprob))
  return(-2*loglik+2*EN+p*log(n))
}

#Gamma-Normal Mixture Model
gammaNormMix = function(data, thresh = 1e-7, maxiter = 10000, 
                        removeZeroes = TRUE, plot = TRUE, hist=TRUE,hist_col="light cyan",
                        verbose = FALSE, forceExponential = FALSE,
                        calculateAreaDifference = FALSE,
                        minDataPoints = 5,onlyAddCurves=FALSE,
                        addContextData = FALSE, contextData = NULL,ylim=c(0,0.6),xlim=c(0,20),xlab=NA,main=NA,breaks=50) {
  
  require(distr)
  require(scales)

  
  #fitting a 2 component normal and gamma mixture model
  
  # add other data to fit the model with as well, but only return 
  # the classification for those we're interested in 
  
  if (addContextData) {
    nOriginal = length(data)
    data <- c(data,contextData)
  }
  
  # assume all values exactly zero already belong to the gamma comp
  # and remove them from the EM algorithm
  
  if (removeZeroes) {
    nonZeroInd = which(data>0)
    x = data[nonZeroInd]
  } else {
    x = data
  }
  
  if(length(x)<minDataPoints) {
    if (verbose) cat("Not enough data points to fit mixture model!")
    return(NA)
  }
  
  # initiate
  n = length(x)
  z = rbinom(n,1,0.5)
  z_iter = z
  mu = -100
  mu_iter = 10
  sig2 = -100
  sig2_iter = 0
  alpha = -100
  alpha_iter = 1
  beta = -100
  beta_iter = 1
  rho = -100
  rho_iter = 0.5
  niter = 0
  
  
  
  while(any(c(
    abs(mu - mu_iter)>thresh ,
    abs(sig2 - sig2_iter)>thresh ,
    abs(alpha - alpha_iter)>thresh ,
    abs(beta - beta_iter)>thresh ,
    abs(rho - rho_iter)>thresh
  )) & (niter<maxiter) ) {
    
    # save old parameters
    mu = mu_iter
    sig2 = sig2_iter
    alpha = alpha_iter
    beta = beta_iter
    rho = rho_iter
    if (forceExponential) alpha_iter = 1 # SHILA
    
    niter = niter + 1
    
    # M step
    mu_iter = sum(z_iter*x)/sum(z_iter)
    sig2_iter = sum(z_iter*(x-mu_iter)*(x-mu_iter))/sum(z_iter)
    if (sig2_iter<=0 | is.na(sig2_iter)) sig2_iter = 1e-11
    beta_iter = alpha_iter*sum(1-z_iter)/sum((1-z_iter)*x)
    if (beta_iter<=0 | is.na(beta_iter)) beta_iter = 3
    if (!forceExponential) {
      alpha_iter = igamma(sum((log(beta_iter) + log(x))*(1-z_iter))/sum(1-z_iter))
    }
    if (alpha_iter>150 | is.na(alpha_iter)) alpha_iter=150
    rho_iter = sum(z_iter)/n
    
    
    # E step
    eta_iter =  
      - 0.5*log(2*pi*sig2_iter) - 
      ((x-mu_iter)*(x-mu_iter))/(2*sig2_iter) - 
      alpha_iter*log(beta_iter) + 
      log(gamma(alpha_iter)) - 
      (alpha_iter-1)*log(x) + 
      beta_iter*x +
      log(rho_iter/(1-rho_iter))
    z_iter = 1/(1+exp(-eta_iter))
    
    if (verbose) cat(niter,mu_iter,sqrt(sig2_iter),alpha_iter,beta_iter,rho_iter,"\n")
  }
  
  
  ll<-sum(log(rho_iter*dnorm(x,mu_iter,sqrt(sig2_iter)) + 
                (1-rho_iter)*dgamma(x,shape=alpha_iter,rate=beta_iter)))
  
  
  xg <- seq(0,max(x)+1,length.out=300)
  c1g <- rho_iter*dnorm(xg,mu_iter,sqrt(sig2_iter))
  c2g <- (1-rho_iter)*dgamma(xg,shape=alpha_iter,rate=beta_iter)
  fg <- rho_iter*dnorm(xg,mu_iter,sqrt(sig2_iter)) + 
    (1-rho_iter)*dgamma(xg,shape=alpha_iter,rate=beta_iter)
  if (plot) {
    if(hist){
      hist(x,probability=TRUE,col=hist_col,breaks=breaks,main=main,xlab=xlab,ylab="Density (zeroes removed)",ylim=ylim,xlim=xlim,lwd=2,cex.lab=1.5)
    }

    lines(xg,c1g,col=alpha("red",0.8),lwd=2)#Normal Lines
    lines(xg,c2g,col=alpha("blue",0.8),lwd=2)#Gamma lines
    lines(xg,fg,col=alpha("black",0.8),lwd=2)#Mixture model line
    if (onlyAddCurves) return(list(xg=xg,c1g=c1g,c2g=c2g,fg=fg))
  }
  if (calculateAreaDifference) {
    f1 <- approxfun(xg, approxfun(density(x,from=0))(xg)-fg)
    # piecewise linear function
    f2 <- function(x) abs(f1(x))
    # take the positive value
    AreaDifference = integrate(f2, min(x[x!=0]), max(x))$value
  } else {
    AreaDifference = NULL
  }
  
  if (removeZeroes) {
    z = rep(0,length(data))
    z[nonZeroInd] <- z_iter
  } else {
    z = z_iter
  }
  
  # force prob expression values above the max to stay the same value
  maxdata = data[which.max(z)]
  z[which(data>maxdata)] <- max(z)
  
  
  
  if (addContextData) {
    z <- z[1:nOriginal]
  }
  if (plot){
    if(addContextData){
      points(data[1:nOriginal],z*0,pch="|",cex=1,col=alpha(rgb(z,0,1-z),0.4))
    }else{
      points(data,z*0,pch="|",cex=1,col=alpha(rgb(z,0,1-z),0.4))
    }
  }
  model_bic<-bic(ll,n,5)
  model_aic<-aic(ll,5)
  model_icl_bic<-icl_bic(ll,z,n,5)
  return(list(probExpressed = z,
              propExpressed = n*rho_iter/length(data),
              numExpressed=length(which(z>0.5)),
              mu = mu_iter,
              sd = sqrt(sig2_iter),
              alpha = alpha_iter,
              beta = beta_iter,
              rho = rho_iter,
              niter = niter,
              loglik=ll,
              BIC=model_bic,
              AIC=model_aic,
              ICL_BIC=model_icl_bic,
              AreaDifference=AreaDifference))
}




#Constrained Gamma-Normal Mixture Model
gammaNormFix = function(data, alpha_fix=NULL,beta_fix=NULL,thresh = 1e-7, maxiter = 10000, 
                        removeZeroes = TRUE, plot = TRUE, hist=TRUE,hist_col="light cyan",
                        verbose = FALSE, forceExponential = FALSE,
                        calculateAreaDifference = FALSE,
                        minDataPoints = 5,onlyAddCurves=FALSE,xlab=NULL,
                        addContextData = FALSE, contextData = NULL,xlim=c(0,15),ylim=c(0,0.6),breaks=50,main=NA) {
  
  require(distr)
  require(scales)
  
  #fitting a 2 component normal and gamma mixture model
  
  # add other data to fit the model with as well, but only return 
  # the classification for those we're interested in 
  
  if (addContextData) {
    nOriginal = length(data)
    data <- c(data,contextData)
  }
  
  # assume all values exactly zero already belong to the gamma comp
  # and remove them from the EM algorithm
  
  if (removeZeroes) {
    nonZeroInd = which(data>0)
    x = data[nonZeroInd]
  } else {
    x = data
  }
  
  if(length(x)<minDataPoints) {
    if (verbose) cat("Not enough data points to fit mixture model!")
    return(NA)
  }
  
  # initiate
  n = length(x)
  z = rbinom(n,1,0.5)
  z_iter = z
  mu = -100
  mu_iter = 10
  sig2 = -100
  sig2_iter = 0
  alpha = alpha_fix
  beta = beta_fix
  rho = -100
  rho_iter = 0.5
  niter = 0
  
  while(any(c(
    abs(mu - mu_iter)>thresh ,
    abs(sig2 - sig2_iter)>thresh ,
    abs(rho - rho_iter)>thresh
  )) & (niter<maxiter) ) {
    
    # save old parameters
    mu = mu_iter
    sig2 = sig2_iter
    rho = rho_iter
    
    niter = niter + 1
    
    # M step
    mu_iter = sum(z_iter*x)/sum(z_iter)
    if (mu_iter<=0 | is.na(mu_iter)) mu_iter=mu
    sig2_iter = sum(z_iter*(x-mu_iter)*(x-mu_iter))/sum(z_iter)
    if (sig2_iter<=0 | is.na(sig2_iter)) sig2_iter = 1e-11
    rho_iter = sum(z_iter)/n
    
    
    # E step
    eta_iter =  
      - 0.5*log(2*pi*sig2_iter) - 
      ((x-mu_iter)*(x-mu_iter))/(2*sig2_iter) - 
      alpha*log(beta) + 
      log(gamma(alpha)) - 
      (alpha-1)*log(x) + 
      beta*x +
      log(rho_iter/(1-rho_iter))
    z_iter = 1/(1+exp(-eta_iter))
    
    if (verbose) cat(niter,mu_iter,sqrt(sig2_iter),alpha,beta,rho_iter,"\n")
  }
  
  ll<-sum(log(rho_iter*dnorm(x,mu_iter,sqrt(sig2_iter)) + 
                (1-rho_iter)*dgamma(x,shape=alpha,rate=beta)))
  
  
  xg <- seq(0,max(x)+1,length.out=300)
  c1g <- rho_iter*dnorm(xg,mu_iter,sqrt(sig2_iter))
  c2g <- (1-rho_iter)*dgamma(xg,shape=alpha,rate=beta)
  fg <- rho_iter*dnorm(xg,mu_iter,sqrt(sig2_iter)) + 
    (1-rho_iter)*dgamma(xg,shape=alpha,rate=beta)
  if (plot) {
    if(hist){
      hist(x,probability=TRUE,col=hist_col,breaks=breaks,main=main,xlab=xlab,ylab="Density (zeroes removed)",ylim=ylim,xlim=xlim,lwd=2,cex.lab=1.5)
    }
    lines(xg,c1g,col=alpha("red",0.8),lwd=2)#Normal Lines
    lines(xg,c2g,col=alpha("blue",0.8),lwd=2)#Gamma lines
    lines(xg,fg,col=alpha("black",0.8),lwd=2)#Mixture model line
    if (onlyAddCurves) return(list(xg=xg,c1g=c1g,c2g=c2g,fg=fg))
  }
  if (calculateAreaDifference) {
    f1 <- approxfun(xg, approxfun(density(x,from=0))(xg)-fg)
    # piecewise linear function
    f2 <- function(x) abs(f1(x))
    # take the positive value
    AreaDifference = integrate(f2, min(x[x!=0]), max(x))$value
  } else {
    AreaDifference = NULL
  }
  
  if (removeZeroes) {
    z = rep(0,length(data))
    z[nonZeroInd] <- z_iter
  } else {
    z = z_iter
  }
  
  # force prob expression values above the max to stay the same value
  maxdata = data[which.max(z)]
  z[which(data>maxdata)] <- max(z)
  
  
  
  if (addContextData) {
    z <- z[1:nOriginal]
  }
  if (plot){
    if(addContextData){
      points(data[1:nOriginal],z*0,pch="|",cex=1,col=alpha(rgb(z,0,1-z),0.4))
    }else{
      points(data,z*0,pch="|",cex=1,col=alpha(rgb(z,0,1-z),0.4))
    }
  }
  
  
  model_bic<-bic(ll,n,3)
  model_aic<-aic(ll,3)
  model_icl_bic<-icl_bic(ll,z,n,3)
  return(list(probExpressed = z,
              propExpressed = n*rho_iter/length(data),
              numExpressed=length(which(z>0.5)),
              mu = mu_iter,
              sd = sqrt(sig2_iter),
              alpha = alpha,
              beta = beta,
              rho = rho_iter,
              niter = niter,
              loglik=ll,
              BIC=model_bic,
              AIC=model_aic,
              ICL_BIC=model_icl_bic,
              AreaDifference=AreaDifference))
}


#Initialisation function for Gamma Mixture Model (Modified from gammamix.init() of the R package mixtools)
gammamix.init <- function(x, lambda = NULL, alpha = NULL, beta = NULL, k = 2,cv=TRUE){
  require(cvTools)
  require(mclust)
  x<-x[x>0]
  n <- length(x)

  if(cv==TRUE){
    if(n>=10){
      CV<-cvFolds(n,K=10)
      km_subset<-c()
      for(i in 1:10){
        km_vec<-rep(NA,n)
        km_vec[CV$which!=i]<-kmeans(x[CV$which!=i],k,nstart=2)$cluster
        km_subset<-cbind(km_subset,km_vec)
      }
      km<-round(apply(km_subset,1,mean,na.rm=TRUE))
    }else{
      CV<-cvFolds(n,K=n)
      km_subset<-c()
      for(i in 1:10){
        km_subset<-cbind(km_subset,kmeans(x[CV$which!=i],k,nstart=2)$cluster)
      }
      km<-round(apply(km_subset,1,mean))
    }
    
  }else{
    km<-kmeans(x,k)$cluster
  }
  
  
  if (is.null(lambda)) {
    lambda = table(km)
    lambda = lambda/sum(lambda)
  } else k = length(lambda)
  
  if(k==1){
    x.bar=mean(x)
    x2.bar=mean(x^2)
  } else{
    x.part=list()
    for(j in 1:k){
      x.part[[j]]=x[which(km==j)]
    }
    x.bar=sapply(x.part,mean)
    x.var=sapply(x.part,var)
  }
  if(is.null(alpha)){
    alpha=x.bar^2/x.var
  }
  
  if(is.null(beta)){
    beta=x.bar/x.var
  }
  
  for(i in 1:k){
    if(length(table(km))==k&table(km)[i]==1){
      alpha[i]=x[which(km==i)]
      beta[i]=1
    }
  }
  
  list(lambda=lambda, alpha=alpha, beta=beta, k=k)
}



#Gamma-Gamma Mixture Model
gammaMix = function(data, thresh = 1e-8, maxiter = 10000, 
                    alpha1_cons=NULL,beta1_cons=NULL,alpha2_cons=NULL,beta2_cons=NULL,
                    removeZeroes = TRUE, plot = TRUE, hist=TRUE,hist_col="light cyan",
                    verbose = FALSE, forceExponential = FALSE,
                    calculateAreaDifference = FALSE,
                    minDataPoints = 5,onlyAddCurves=FALSE,
                    addContextData = FALSE, contextData = NULL,xlim=c(0,12),ylim=c(0,0.6),breaks=50,main=NA) {
  
  require(distr)
  require(scales)
  
  #fitting a 2 component normal and gamma mixture model
  
  # add other data to fit the model with as well, but only return 
  # the classification for those we're interested in 
  
  if (addContextData) {
    nOriginal = length(data)
    data <- c(data,contextData)
  }
  
  # assume all values exactly zero already belong to the gamma comp
  # and remove them from the EM algorithm
  
  if (removeZeroes) {
    nonZeroInd = which(data>0)
    x = data[nonZeroInd]
  } else {
    x = data
  }
  
  if(length(x)<minDataPoints) {
    if (verbose) cat("Not enough data points to fit mixture model!")
    return(NA)
  }
  
  est<-gammamix.init(x)
  while(est$lambda[1]==1){
    est<-gammamix.init(x)
  }
  #print(est)
  est.min<-1
  est.max<-2
  # initiate
  n = length(x)
  
  alpha1 = -100
  alpha1_iter = est$alpha[est.min]
  beta1 = -100
  beta1_iter = est$beta[est.min]
  alpha2 = -100
  alpha2_iter = est$alpha[est.max]
  beta2 = -100
  beta2_iter = est$beta[est.max]
  rho1 = -100
  rho1_iter = est$lambda[est.min]
  rho2 = 100
  rho2_iter = est$lambda[est.max]
  niter = 0
  diff_ll<-thresh+1
  
  z = rbinom(n,1,0.5)
  z1_iter = z
  z2_iter = z
  ll_old<-sum(log(rho2_iter*dgamma(x,shape=alpha2_iter,rate=beta2_iter) + 
                    rho1_iter*dgamma(x,shape=alpha1_iter,rate=beta1_iter)))

  if(is.null(alpha1_cons)){
    if(alpha1_iter>=100|alpha2_iter>=100){
      alpha2_cons=500
      alpha1_cons=100
    }else{
      alpha2_cons=50
      alpha1_cons=20
    }
  }
  
  while(any(c(
    abs(alpha1 - alpha1_iter)>thresh ,
    abs(beta1 - beta1_iter)>thresh ,
    abs(alpha2 - alpha2_iter)>thresh ,
    abs(beta2 - beta2_iter)>thresh ,
    abs(rho1 - rho1_iter)>thresh,
    abs(rho2 - rho2_iter)>thresh
    #diff_ll>thresh
  )) & (niter<maxiter) ) {
    
    # save old parameters
    
    alpha1 = alpha1_iter
    beta1 = beta1_iter
    alpha2 = alpha2_iter
    beta2 = beta2_iter
    rho1 = rho1_iter
    rho2 = rho2_iter
    if (forceExponential) alpha1_iter = 1 
    niter = niter + 1
    
    
    #M step
    beta1_iter = alpha1_iter*sum(1-z1_iter)/sum((1-z1_iter)*x)
    if (beta1_iter<=0 | is.na(beta1_iter)) beta1_iter =3
    if (!forceExponential) {
      alpha1_iter = igamma(sum((log(beta1_iter) + log(x))*(1-z1_iter))/sum(1-z1_iter))
    }
    #if (is.na(alpha1_iter)) alpha1_iter=alpha1
    if (alpha1_iter>alpha1_cons | is.na(alpha1_iter)) alpha1_iter=alpha1
    
    beta2_iter = alpha2_iter*sum(z1_iter)/sum((z1_iter)*x)
    if (beta2_iter<=0 | is.na(beta2_iter)) beta2_iter = 3
    alpha2_iter = igamma(sum((log(beta2_iter) + log(x))*(z1_iter))/sum(z1_iter))
    #if (is.na(alpha2_iter)) alpha2_iter=alpha2
    if (alpha2_iter>alpha2_cons | is.na(alpha2_iter)) alpha2_iter=alpha2
    rho2_iter = sum(z1_iter[is.finite(z1_iter)])/n
    rho1_iter = sum((1-z1_iter[is.finite(z1_iter)]))/n
    dgamma2<-dgamma(x,shape=alpha2_iter,rate=beta2_iter)
    dgamma2[dgamma2==0]<-0.0000001
    dgamma1<-dgamma(x,shape=alpha1_iter,rate=beta1_iter)
    dgamma1[dgamma1==0]<-0.0000001
    ll_new<-sum(log(rho2_iter*dgamma2 + 
                      rho1_iter*dgamma1))
    diff_ll<-ll_new-ll_old
    ll_old<-ll_new
    # # E step
    
    if(alpha1_iter>160){
      logGammaAlpha1<-(alpha1_iter-1)*log(alpha1_iter-1)-(alpha1_iter-1)
    }else{
      logGammaAlpha1<-log(gamma(alpha1_iter))
    }
    if(alpha2_iter>160){
      logGammaAlpha2<-(alpha2_iter-1)*log(alpha2_iter-1)-(alpha2_iter-1)
    }else{
      logGammaAlpha2<-log(gamma(alpha2_iter))
    }
    eta_iter =
      alpha1_iter*log(beta1_iter) -
      alpha2_iter*log(beta2_iter) -
      logGammaAlpha1 +
      logGammaAlpha2 +
      (alpha1_iter-alpha2_iter)*log(x) -
      (beta1_iter-beta2_iter)*x +
      log(rho1_iter/rho2_iter)
    z1_iter = 1/(1+exp(eta_iter))
    
    
    if (verbose) cat(niter,alpha1_iter,beta1_iter,alpha2_iter,beta2_iter,rho1_iter,rho2_iter,"\n")
  }
  
  dgamma2<-dgamma(x,shape=alpha2_iter,rate=beta2_iter)
  dgamma2[dgamma2==0]<-0.0000001
  dgamma1<-dgamma(x,shape=alpha1_iter,rate=beta1_iter)
  dgamma1[dgamma1==0]<-0.0000001
  
  ll<-sum(log(rho2_iter*dgamma2 + 
                rho1_iter*dgamma1))
  if(alpha1_iter/beta1_iter>alpha2_iter/beta2_iter){
    temp<-rho1_iter
    rho1_iter<-rho2_iter
    rho2_iter<-temp
    temp<-alpha1_iter
    alpha1_iter<-alpha2_iter
    alpha2_iter<-temp
    temp<-beta1_iter
    beta1_iter<-beta2_iter
    beta2_iter<-temp
    z1_iter<-1-z1_iter
  }
  if(round(rho1_iter,1)==1){
    temp<-rho1_iter
    rho1_iter<-rho2_iter
    rho2_iter<-temp
    temp<-alpha1_iter
    alpha1_iter<-alpha2_iter
    alpha2_iter<-temp
    temp<-beta1_iter
    beta1_iter<-beta2_iter
    beta2_iter<-temp
    z1_iter<-1-z1_iter
  }
  xg <- seq(0,max(x)+1,length.out=300)
  c1g <- rho2_iter*dgamma(xg,shape=alpha2_iter,rate=beta2_iter) #Gamma2 Lines
  c2g <- rho1_iter*dgamma(xg,shape=alpha1_iter,rate=beta1_iter) #Gamma1 Lines
  fg <- rho2_iter*dgamma(xg,shape=alpha2_iter,rate=beta2_iter) + 
    rho1_iter*dgamma(xg,shape=alpha1_iter,rate=beta1_iter)

  if (plot) {
    if(hist){
      hist(x,probability=TRUE,col=hist_col,breaks=breaks,main=main,xlab=NA,ylab="Density (zeroes removed)",ylim=ylim,xlim=xlim,lwd=2,cex.lab=1.5)
    }

    lines(xg,c1g,col=alpha("red",0.8),lwd=2)#Gamma2 Lines
    lines(xg,c2g,col=alpha("blue",0.8),lwd=2)#Gamma1 lines
    lines(xg,fg,col=alpha("black",0.8),lwd=2)#Mixture model line
    if (onlyAddCurves) return(list(xg=xg,c1g=c1g,c2g=c2g,fg=fg))
  }
  
  if (calculateAreaDifference) {
    f1 <- approxfun(xg, approxfun(density(x,from=0))(xg)-fg)
    # piecewise linear function
    f2 <- function(x) abs(f1(x))
    # take the positive value
    AreaDifference = integrate(f2, min(x[x!=0]), max(x))$value
  } else {
    AreaDifference = NULL
  }
  
  if (removeZeroes) {
    z = rep(0,length(data))
    z[nonZeroInd] <- z1_iter
  } else {
    z = z1_iter
  }
  
  # force prob expression values above the max to stay the same value
  maxdata = data[which.max(z)]
  z[which(data>maxdata)] <- max(z)
  model_bic<-bic(ll,n,5)
  model_aic<-aic(ll,5)
  model_icl_bic<-icl_bic(ll,z,n,5)
  if (addContextData) {
    z <- z[1:nOriginal]
  }
  if (plot){
    if(addContextData){
      points(data[1:nOriginal],z*0,pch="|",cex=1,col=alpha(rgb(z,0,1-z),0.4))
    }else{
      points(data,z*0,pch="|",cex=1,col=alpha(rgb(z,0,1-z),0.4))
    }
  }
  return(list(probExpressed = z,
              propExpressed = n*rho2_iter/length(data),
              numExpressed=length(which(z>0.5)),
              alpha1 = alpha1_iter,
              beta1 = beta1_iter,
              alpha2 = alpha2_iter,
              beta2 = beta2_iter,
              rho1 = rho1_iter,
              rho2 = rho2_iter,
              niter = niter,
              loglik=ll,
              BIC=model_bic,
              AIC=model_aic,
              ICL_BIC=model_icl_bic,
              AreaDifference=AreaDifference))
}



##Constrained Gamma Mixture Model
gammaMixFix = function(data, alpha_fix=NULL,beta_fix=NULL,thresh = 1e-10, maxiter = 10000, 
                       alpha2_cons=NULL,removeZeroes = TRUE, plot = TRUE, hist=TRUE,hist_col="light cyan",
                       verbose = FALSE, forceExponential = FALSE,
                       calculateAreaDifference = FALSE,
                       minDataPoints = 5,onlyAddCurves=FALSE,
                       addContextData = FALSE, contextData = NULL,xlim=c(0,12),ylim=c(0,0.6),breaks=50,main=NA) {
  require(distr)
  require(scales)
  
  if (addContextData) {
    nOriginal = length(data)
    data <- c(data,contextData)
  }
  
  # assume all values exactly zero already belong to the gamma comp
  # and remove them from the EM algorithm
  if (removeZeroes) {
    nonZeroInd = which(data>0)
    x = data[nonZeroInd]
  } else {
    x = data
  }
  if(length(x)<minDataPoints) {
    if (verbose) cat("Not enough data points to fit mixture model!")
    return(NA)
  }
  est<-gammamix.init(x)
  while(est$lambda[1]==1){
    est<-gammamix.init(x)
  }
  est.max<-which.max(est$alpha/est$beta)
  
  # initiate
  n = length(x)
  z = rbinom(n,1,0.5)
  z1_iter = z
  z2_iter = z
  alpha1 = alpha_fix
  beta1 = beta_fix
  alpha2=-100
  alpha2_iter = est$alpha[est.max]
  beta2 = -100
  beta2_iter = est$beta[est.max]
  rho2 = 100
  rho2_iter = est$lambda[est.max]
  rho1 = -100
  rho1_iter = 1-est$lambda[est.max]
  niter = 0
  
  if(alpha2_iter>1000){
    alpha2_iter=1000
  }

  if(is.null(alpha2_cons)){
    if(alpha2_iter>=100){
      alpha2_cons=500
    }else{
      alpha2_cons=50
    }
  }
  
  while(any(c(
    abs(alpha2 - alpha2_iter)>thresh ,
    abs(beta2 - beta2_iter)>thresh ,
    abs(rho1 - rho1_iter)>thresh,
    abs(rho2 - rho2_iter)>thresh
  )) & (niter<maxiter) ) {
    
    # save old parameters
    alpha2 = alpha2_iter
    beta2 = beta2_iter
    rho1 = rho1_iter
    rho2 = rho2_iter
    if (forceExponential) alpha1_iter = 1
    niter = niter + 1
    
    #M step
    beta2_iter = alpha2_iter*sum(z1_iter)/sum((z1_iter)*x)
    if (beta2_iter<=0 | is.na(beta2_iter)) beta2_iter = 3
    alpha2_iter = igamma(sum((log(beta2) + log(x))*(z1_iter))/sum(z1_iter))
    if (alpha2_iter>alpha2_cons| is.na(alpha2_iter)) alpha2_iter=alpha2
    rho2_iter = sum(z1_iter)/n
    rho1_iter = sum((1-z1_iter))/n
  
    ##E step
    if(alpha1>160){
      logGammaAlpha1<-(alpha1-1)*log(alpha1-1)-(alpha1-1)
    }else{
      logGammaAlpha1<-log(gamma(alpha1))
    }
    if(alpha2_iter>160){
      logGammaAlpha2<-(alpha2_iter-1)*log(alpha2_iter-1)-(alpha2_iter-1)
    }else{
      logGammaAlpha2<-log(gamma(alpha2_iter))
    }
    
    eta_iter =
      alpha1*log(beta1) -
      alpha2_iter*log(beta2_iter) -
      logGammaAlpha1 +
      logGammaAlpha2 +
      (alpha1-alpha2_iter)*log(x) -
      (beta1-beta2_iter)*x +
      log(rho1_iter/rho2_iter)
    z1_iter = 1/(1+exp(eta_iter))
    
    
    if (verbose) cat(niter,alpha2_iter,beta2_iter,rho1_iter,rho2_iter,"\n")
  }
  
  dgamma2<-dgamma(x,shape=alpha2_iter,rate=beta2_iter)
  dgamma2[dgamma2==0]<-0.0000001
  dgamma1<-dgamma(x,shape=alpha1,rate=beta1)
  dgamma1[dgamma1==0]<-0.0000001
  
  ll<-sum(log(rho2_iter*dgamma2 + 
                rho1_iter*dgamma1))
  xg <- seq(0,max(x)+1,length.out=300)
  c1g <- rho2_iter*dgamma(xg,shape=alpha2_iter,rate=beta2_iter) #Gamma2 Lines
  c2g <- rho1_iter*dgamma(xg,shape=alpha1,rate=beta1) #Gamma1 Lines
  fg <- rho2_iter*dgamma(xg,shape=alpha2_iter,rate=beta2_iter) + 
    rho1_iter*dgamma(xg,shape=alpha1,rate=beta1)
  if (plot) {
    if(hist){
      hist(x,probability=TRUE,col=hist_col,breaks=breaks,main=main,xlab=NA,ylab="Density (zeroes removed)",ylim=ylim,xlim=xlim)
    }

    lines(xg,c1g,col=alpha("red",0.8),lwd=1.5)#Gamma2 Lines
    lines(xg,c2g,col=alpha("blue",0.8),lwd=1.5)#Gamma1 lines
    lines(xg,fg,col=alpha("black",0.8),lwd=1.5)#Mixture model line
    if (onlyAddCurves) return(list(xg=xg,c1g=c1g,c2g=c2g,fg=fg))
  }
  if (calculateAreaDifference) {
    f1 <- approxfun(xg, approxfun(density(x,from=0))(xg)-fg)
    # piecewise linear function
    f2 <- function(x) abs(f1(x))
    # take the positive value
    AreaDifference = integrate(f2, min(x[x!=0]), max(x))$value
  } else {
    AreaDifference = NULL
  }
  
  if (removeZeroes) {
    z = rep(0,length(data))
    z[nonZeroInd] <- z1_iter
  } else {
    z = z1_iter
  }
  
  # force prob expression values above the max to stay the same value
  maxdata = data[which.max(z)]
  z[which(data>maxdata)] <- max(z)
  model_bic<-bic(ll,n,3)
  model_aic<-aic(ll,3)
  model_icl_bic<-icl_bic(ll,z,n,3)  
  
  if (addContextData) {
    z <- z[1:nOriginal]
  }
  if (plot){
    if(addContextData){
      points(data[1:nOriginal],z*0,pch="|",cex=1,col=alpha(rgb(z,0,1-z),0.4))
    }else{
      points(data,z*0,pch="|",cex=1,col=alpha(rgb(z,0,1-z),0.4))
    }
  }
  
  
  return(list(probExpressed = z,
              propExpressed =n*rho2_iter/length(data),
              numExpressed=length(which(z>0.5)),
              alpha1 = alpha1,
              beta1 = beta1,
              alpha2 = alpha2_iter,
              beta2 = beta2_iter,
              rho1 = rho1_iter,
              rho2 = rho2_iter,
              niter = niter,
              loglik=ll,
              BIC=model_bic,
              AIC=model_aic,
              ICL_BIC=model_icl_bic,
              AreaDifference=AreaDifference))
}





