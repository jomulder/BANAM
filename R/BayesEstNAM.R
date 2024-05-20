

#' @title Bayesian estimation of the network autocorrelation model
#' @description The \code{banam} function can be used for Bayesian estimation of the
#' network autocorrelation model (NAM). In the case of a single weight matrix, a flat prior,
#' the independence Jeffreys prior, and a normal prior can be specified for the network autocorrelation
#' parameter. In the case of multiple weight matrices, a multivariate normal prior can be specified.
#' @param y A numeric vector containing the observations of the outcome variable.
#' @param X The design matrix of the predictor variables. If absent a column of ones
#' is automatically to model the intercept.
#' @param W A weight matrix (in the case of a NAM with a single weight matrix) or
#' a list of weight matrices (in the case of a NAM with multiple weight matrices).
#' @param prior A character string specifying which prior to use in the case of a
#' NAM with a single weight matrix. The options are 'flat', 'IJ', and 'normal',
#' for a flat prior, the independence Jeffreys prior, and a normal prior, respectively.
#' @param prior.mean A scalar (or vector) specifying the prior mean(s) of the network
#' autocorrelation(s). The default prior mean is 0.
#' @param prior.Sigma A scalar (or matrix) specifying the prior variance (or prior covariance
#' matrix) of the network autocorrelation(s). In the univariate case, the default prior variance
#' is 1e6. In the multivariate case, the default prior covariance matrix is the identity matrix
#' times 1e6.
#' @param postdraws An integer specifying the number of posterior draws after burn-in.
#' @param burnin An integer specifying the number of draws for burn-in.
#'
#' @return The output is an object of class \code{banam}. For users of \code{BANAM}, the following
#'         are the useful objects:
#' \itemize{
#' \item \code{rho.draws} Matrix of posterior draws for the network autocorrelation parameter(s).
#' \item \code{beta.draws} Matrix of posterior draws for the coefficients.
#' \item \code{sigma2.draws} Matrix of posterior draws for the error variance.
#' \item \code{summarystats} Table with summary statistics of the posterior.
#' \item \code{residuals} Residuals based on all posterior draws.
#' \item \code{fitted.values} Fitted values based on all posterior draws.
#' }
#' @references Dittrich, D., Leenders, R.Th.A.J., & Mulder, J. (2017).
#' Bayesian estimation of the network autocorrelation model. Social Network, 48, 213–236.
#' \url{https://doi.org/10.1016/j.socnet.2016.09.002}
#' @references Dittrich, D., Leenders, R.Th.A.J., & Mulder, J. (2020). Network
#' Autocorrelation Modeling: Bayesian Techniques for Estimating and Testing Multiple
#' Network Autocorrelations. Sociological Methodology, 50, 168-214.
#' \url{https://doi.org/10.1177/0081175020913899}
#' @examples
#' \donttest{
#' #example analyses
#' #generate example data
#' set.seed(234)
#' n <- 50
#' d1 <- .2
#' Wadj1 <- sna::rgraph(n, tprob=d1, mode="graph")
#' W1 <- sna::make.stochastic(Wadj1, mode="row")
#' d2 <- .4
#' Wadj2 <- sna::rgraph(n, tprob=d2, mode="graph")
#' W2 <- sna::make.stochastic(Wadj2, mode="row")
#' # set rho, beta, sigma2, and generate y
#' rho1 <- .3
#' K <- 3
#' beta <- rnorm(K)
#' sigma2 <- 1
#' X <- matrix(c(rep(1, n), rnorm(n*(K-1))), nrow=n, ncol=K)
#' y <- c((solve(diag(n) - rho1*W1))%*%(X%*%beta + rnorm(n)))
#'
#' #Bayesian estimation of NAM with a single weight matrix using a flat prior for rho
#' best1 <- banam(y,X,W1)
#' print(best1)
#'
#' #Bayesian estimation of NAM with two weight matrices using standard normal priors
#' best2 <- banam(y,X,W=list(W1,W2))
#' print(best2)
#'
#' #Bayes factor testing of equality/order hypotheses using environment of package 'BFpack'
#' BFbest2 <- BF(best2,hypothesis="rho1>rho2>0; rho1=rho2>0; rho1=rho2=0")
#' }
#' @export
#' @rdname banam
banam <- function(y,X,W,prior="flat",prior.mean=NULL,prior.Sigma=NULL,postdraws=5e3,
                    burnin=1e3){

  priormean <- prior.mean
  priorsigma <- prior.Sigma

  #check matrix of covariates X, and if missing add column of ones for intercept
  if(is.null(X)){
    X <- matrix(1,nrow=length(y),ncol=1)
  }else{
    if(is.data.frame(X)){X <- as.matrix(X)}
    if(is.matrix(X)<1){stop("If specified, X must be a matrix.")}
    if(!(mean(X[,1])==1 & sd(X[,1])==0)){
      #add vector of ones for intercept
      X <- cbind(1,X)
    }
  }

  if(is.null(X)){stop("The design maatrix 'X' cannot be empty.")}
  if(is.null(y)){stop("The response vector 'y' cannot be empty.")}
  if(is.null(W)){stop("The response vector 'W' cannot be empty.")}
  if(is.data.frame(X)){X <- as.matrix(X)}
  if(is.data.frame(y)){y <- c(as.matrix(y))}
  if(is.data.frame(W)){W <- as.matrix(W)}
  if(!(is.matrix(W) | is.list(W))){stop("The weight matrix W should be a g x g matrix or a list of weight matrices.")}
  if(!is.vector(y)){stop("The response vector 'y' must be a vector.")}
  if(!is.numeric(y)){stop("The response vector 'y' must be numeric.")}
  if(abs(length(y)-nrow(X))>.5){stop("Number of observations in 'X' must match length of 'y'.")}
  if((postdraws%%1==0)<1 | postdraws<1){stop("The number of desired posterior draws 'N' must be a positive integer.")}
  if((burnin%%1==0)<1 | burnin<0){stop("The 'burnin' must be a non-negative integer.")}

  message(paste0("BANAM: Posterior sampling"))

  if(is.list(W) & length(W) > 1){
    # NAM with multiple rho's
    if(is.null(priormean)){
      priormean <- rep(0,length(W))
    }
    if(is.null(priorsigma)){
      priorsigma <- diag(length(W)) * 1e6
    }
    prior <- "multivariate normal"
    Best1 <- MCMC.multiple(y,X,W,mu.prior=priormean,Sigma.prior=priorsigma,postdraws,burnin)
  }else{
    # NAM with single rho
    g <- length(y)
    if(is.list(W)){W <- W[[1]]}
    if(nrow(W)!=length(y) | ncol(W)!=length(y) ){stop("The weight matrix W should be a g x g matrix")}
    if(!is.numeric(W)){stop("The connectivity matrix must be numeric.")}
    if(length(unique(g,nrow(W),ncol(W)))>1){stop("Order of the connectivity matrix must match length of 'y'.")}

    #set starting values
    out.mle <- lm(y ~ -1 + X)
    startval <- c(0,out.mle$coefficients,sigma(out.mle)**2)
    if(prior == "flat"){
      priormean <- NULL
      priorsigma <- NULL
      Best1 <- MCMC.F(N = postdraws + burnin, y, X, W, startval, burn=burnin)
    }else if(prior == "IJ"){
      priormean <- NULL
      priorsigma <- NULL
      Best1 <- MCMC.IJ(N = postdraws + burnin, y, X, W, startval, burn=burnin)
      prior <- "independence Jeffreys"
    }else if(prior == "normal"){
      if(is.null(priormean)){
        priormean <- 0
      }
      if(is.null(priorsigma)){
        priorsigma <- 1e6 #default prior variance
      }
      if(sum(is.vector(priormean)+length(priormean))<2){stop("The prior mean must be a scalar.")}
      if(!is.numeric(priorsigma)){stop("The prior mean must be numeric.")}
      if(sum(is.vector(priorsigma)+length(priorsigma)+as.numeric(priorsigma>0))<3){
        stop("The prior variance must be a positive scalar.")}
      Best1 <- MCMC.N(N = postdraws + burnin, y, X, W, m = priormean,
                      std2 = priorsigma, startval, burn = burnin)
    }
    W <- list(W)
  }
  cat("\n")

  message(paste0("Sampling finished"))

  o <- list()
  o$y <- y
  o$X <- X
  o$W <- W
  o$priormean <- priormean
  o$priorsigma <- priorsigma
  o$postdraws <- postdraws
  o$burnin <- burnin
  o$rho.draws <- Best1$rho.draws
  o$beta.draws <- Best1$beta.draws
  o$sigma2.draws <- Best1$sigma2.draws
  o$rho.mean <- apply(Best1$rho.draws,2,mean)
  if(ncol(Best1$rho.draws)>1){
    paste0(rep("rho",ncol(Best1$rho.draws)),1:(ncol(Best1$beta.draws)))
  }else{
    names(o$rho.mean) <- "rho"
  }
  o$beta.mean <- apply(Best1$beta.draws,2,mean)
  #names(o$beta.mean) <- paste0(rep("beta",ncol(Best1$beta.draws)),1:ncol(Best1$beta.draws))
  if(ncol(Best1$beta.draws)>1){
    names(o$beta.mean) <- c("(Intercept)",
                            paste0(rep("beta",ncol(Best1$beta.draws)-1),1:(ncol(Best1$beta.draws)-1)))
  }else{
    names(o$beta.mean) <- "(Intercept)"
  }
  o$sigma2.mean <- apply(Best1$sigma2.draws,2,mean)
  names(o$sigma2.mean) <- "sigma2"
  o$rho.sd <- apply(Best1$rho.draws,2,sd)
  names(o$rho.sd) <- "rho"
  o$beta.sd <- apply(Best1$beta.draws,2,sd)
  names(o$beta.sd) <- "beta"
  o$sigma2.sd <- apply(Best1$sigma2.draws,2,sd)
  names(o$sigma2.sd) <- "sigma2"
  summarystats <- rbind(cbind(o$rho.mean,
              o$rho.sd,
              apply(o$rho.draws>0,2,mean),
              apply(o$rho.draws,2,quantile,.025),
              apply(o$rho.draws,2,quantile,.975)),
        cbind(o$beta.mean,
              o$beta.sd,
              apply(o$beta.draws>0,2,mean),
              apply(o$beta.draws,2,quantile,.025),
              apply(o$beta.draws,2,quantile,.975)),
        cbind(o$sigma2.mean,
              o$sigma2.sd,
              1,
              quantile(o$sigma2.draws,.025),
              quantile(o$sigma2.draws,.975))
  )
  colnames(summarystats) <- c("Post.mean","Post.sd","Pr(>0)","2.5%","97.5%")
  o$summarystats <- summarystats
  o$fitted.values <- Best1$fitted.vals
  o$residuals <- Best1$res
  o$prior <- prior
  o$prior.mean <- priormean
  o$prior.Sigma <- priorsigma
  o$call <- match.call()
  class(o) <- c("banam")
  o

  return(o)

}


# Help functions
loglik <- function(rho, EW, g, yMy, yMWy, yWMWy){
  (-1)*(sum(log(1 - rho*EW)) - (g/2)*log(yMy - 2*rho*yMWy + rho**2*yWMWy))}
loglik.gr <- function(rho, EW, g, yMy, yMWy, yWMWy){
  (-1)*((-1)*sum(EW/(1 - rho*EW)) + g*(yMWy - rho*yWMWy)/(yMy - 2*rho*yMWy + rho**2*yWMWy))}

# helpers for MCMC samplers
# for flat prior
mh.up.rhobeta1.F <- function(mu, Sigma, lb, ub, y, Wy, EW, lndet, s2, Xbetatilde, g,
                             SS, cval, Ay) {
  pval <- c(rtmvnorm(1, mean=mu, sigma=Sigma, lower=c(lb, -Inf), upper=c(ub, Inf)))
  Ay_pval <- y-pval[1]*Wy
  lndet_pval <- sum(log(1 - pval[1]*EW))
  if(log(runif(1)) <
     lndet_pval - lndet
     - (1/(2*s2))*(sum(Ay_pval**2) - 2*pval[2]*sum(Ay_pval) - 2*sum(Ay_pval*Xbetatilde)
                   + g*pval[2]**2 + 2*pval[2]*sum(Xbetatilde) + sum(Xbetatilde**2) - SS) +
     dmvnorm(cval, mean=mu, sigma=Sigma, log=T) -
     dmvnorm(pval, mean=mu, sigma=Sigma, log=T)){
    val <- pval
    Ay <- Ay_pval
    lndet  <- lndet_pval
  }
  else{
    val <- cval
    Ay <- Ay
    lndet <- lndet
  }
  outlist = list(val, Ay, lndet)
  names(outlist) = c("value", "Ay", "lndet")
  return(outlist)
}
# for normal prior
mh.up.rhobeta1.N <- function(mu, Sigma, lb, ub, y, Wy, EW, lndet, s2, Xbetatilde, g,
                             SS, cval, Ay, m, std2) {
  pval <- c(rtmvnorm(1, mean=mu, sigma=Sigma, lower=c(lb, -Inf), upper=c(ub, Inf)))
  Ay_pval <- y-pval[1]*Wy
  lndet_pval <- sum(log(1 - pval[1]*EW))
  if(log(runif(1))
     < lndet_pval - lndet
     - (1/(2*s2))*(sum(Ay_pval**2) - 2*pval[2]*sum(Ay_pval) - 2*sum(Ay_pval*Xbetatilde) +
                   g*pval[2]**2 + 2*pval[2]*sum(Xbetatilde) + sum(Xbetatilde**2) - SS)
     - (1/(2*std2))*(((pval[1] - m)**2) - ((cval[1] - m)**2)) +
     dmvnorm(cval, mean=mu, sigma=Sigma, log=T) -
     dmvnorm(pval, mean=mu, sigma=Sigma, log=T)){
    val <- pval
    Ay <- Ay_pval
    lndet <- lndet_pval
  }else{
    val <- cval
    Ay <- Ay
    lndet <- lndet
  }
  outlist <- list(val, Ay, lndet)
  names(outlist) <- c("value", "Ay", "lndet")
  return(outlist)}
# for independence Jeffreys prior
mh.up.rhobeta1.IJ <- function(mu, Sigma, lb, ub, Id, W, betatilde, X, EW, cval, s2, y,
                              SS, g, trBtB, trB2, SSBXbeta_cval, trB){
  pval <- c(rtmvnorm(1, mean=mu, sigma=Sigma, lower=c(lb, -Inf), upper=c(ub, Inf)))
  A_pval <- Id - pval[1]*W
  B_pval <- W%*%solve(A_pval)
  beta_pval <- c(pval[2], betatilde)
  SSBXbeta_pval <- sum(c(B_pval%*%X%*%beta_pval)**2)
  if(runif(1) <=
     exp(sum(1 - pval[1]*EW) - sum(1 - cval[1]*EW)
         - (1/(2*s2))*(sum((c(A_pval%*%y) - c(X%*%beta_pval))**2) - SS)
         + .5*(log(sum(B_pval**2) + sum((EW/(1 - pval[1]*EW))**2) + SSBXbeta_pval/s2 - (2/g)*(tr(B_pval)**2))
               - log(trBtB + trB2 + SSBXbeta_cval/s2 - (2/g)*(trB**2))) +
         dmvnorm(cval, mean=mu, sigma=Sigma, log=T) -
         dmvnorm(pval, mean=mu, sigma=Sigma, log=T))){
    pval
  }else{cval}}
mh.up.betatilde.IJ <- function(XtXiXt, Ay, s2, XtXi, k, beta1, B, X, SS, EW, g, trBtB,
                               trB2, SSBXbeta_cval, trB, cval){
  mu_beta <- c(XtXiXt%*%Ay)
  Sigma_beta <- s2*XtXi
  Sigma_beta11 <- c(Sigma_beta[1:1])
  Sigma_beta12 <- c(Sigma_beta[1, 2:k])
  Sigma_beta22 <- Sigma_beta[2:k, 2:k]
  mu_betatilde <- mu_beta[-1] + Sigma_beta12*(beta1 - mu_beta[1])/Sigma_beta11
  Sigma_betatilde <- Sigma_beta22 - outer(Sigma_beta12, Sigma_beta12)/Sigma_beta11
  pval <- c(rmvnorm(1, mean=mu_betatilde, sigma=Sigma_betatilde))
  beta_pval <- c(beta1, pval)
  SSBXbeta_pval <- sum(c(B%*%X%*%beta_pval)**2)
  if(runif(1)
     <= exp(-(1/(2*s2))*(sum((Ay - c(X%*%beta_pval))**2) - SS)
            + .5*(log(trBtB + trB2 + SSBXbeta_pval/s2 - (2/g)*(trB**2))
                  - log(trBtB + trB2 + SSBXbeta_cval/s2 - (2/g)*trB**2)) +
            dmvnorm(cval, mean=mu_betatilde, sigma=Sigma_betatilde, log=T) -
            dmvnorm(pval, mean=mu_betatilde, sigma=Sigma_betatilde, log=T))){
    pval}
  else{cval}
}
mh.up.s2.IJ <- function(g, SS, cval, trBtB, trB2, SSBXbeta, trB){
  pval <- rinvgamma(1, alpha=(g+1)/2, beta=SS/2)
  if(runif(1) <=
     exp((-(g + 2)/2)*(log(pval) - log(cval)) - .5*SS*((1/pval) - (1/cval))
         + .5*(log(trBtB + trB2 + (1/pval)*SSBXbeta - (2/g)*(trB**2))
               - log(trBtB + trB2 + (1/cval)*SSBXbeta - (2/g)*(trB**2)))
         + log(dinvgamma(cval, alpha=(g+1)/2, beta=SS/2))
         - log(dinvgamma(pval, alpha=(g+1)/2, beta=SS/2)))){
    pval
  }else{cval}
}
# for NAM with multiple rho's
in.space.2.rsd <- function(rho,Wlist,R){
  inout <- 0
  if(sum(rho>=0)>1.5){if(sum(rho)<1){inout=1}else{inout=0}}
  else{
    mf <- c(1,rho[2]/rho[1])
    W <- Wlist[[1]]+mf[2]*Wlist[[2]]
    if(rho[1]>0){
      rho1b <- 1/Re(eigs(W,1,which="LR",opts=list(retvec=FALSE))$values)
      if(identical(rho1b,complex(0))>0){rho1b=1/max(Re(eigen(W,only.values=TRUE)$values))}
    }
    else{
      rho1b=1/Re(eigs(W,1,which="SR",opts=list(retvec=FALSE))$values)
      if(identical(rho1b,complex(0))>0){rho1b=1/min(Re(eigen(W,only.values=TRUE)$values))}
    }
    if(sum(abs(rho)<=abs(mf*rho1b))==2){inout=1}
  }
  return(inout)
}
in.space.2 <- function(rho,Wlist,R){
  inout <- 0
  mf <- c(1,rho[2]/rho[1])
  W <- Wlist[[1]]+mf[2]*Wlist[[2]]
  if(rho[1]>0){
    rho1b <- 1/Re(eigs(W,1,which="LR",opts=list(retvec=FALSE))$values)
    if(identical(rho1b,complex(0))>0){
      rho1b <- 1/max(Re(eigen(W,only.values=TRUE)$values))
    }
  }else{
    rho1b <- 1/Re(eigs(W,1,which="SR",opts=list(retvec=FALSE))$values)
    if(identical(rho1b,complex(0))>0){
      rho1b <- 1/min(Re(eigen(W,only.values=TRUE)$values))
    }
  }
  if(sum(abs(rho)<=abs(mf*rho1b))==2){
    inout <- 1
  }
  return(inout)
}
in.space.R <- function(rho,Wlist,R){
  inout <- 0
  alpha <- NULL
  prod.cos <- NULL
  mf <- NULL
  ss <- NULL
  alpha[1] <- atan2(rho[2],rho[1])
  prod.cos[1:2] <- rep(1,2)
  mf[1:2] <- c(1,tan(alpha[1]))
  ss[1] <- rho[1]**2
  ss[2] <- ss[1]+rho[2]**2
  W <- mf[1]*Wlist[[1]]+mf[2]*Wlist[[2]]
  for(i in 3:R){
    ss[i] <- ss[i-1]+rho[i]**2
    if(rho[i]>0){
      alpha[i-1] <- acos(sqrt(ss[i-1]/ss[i]))
    }else{
      alpha[i-1] <- -acos(sqrt(ss[i-1]/ss[i]))
    }
    prod.cos[i] <- prod.cos[i-1]*cos(alpha[i-2])
    mf[i] <- tan(alpha[i-1])/prod.cos[i]
    W <- W + mf[i]*Wlist[[i]]
  }
  if(rho[1]>0){
    rho1b <- 1/Re(eigs(W,1,which="LR",opts=list(retvec=FALSE))$values)
    if(identical(rho1b,complex(0))>0){
      rho1b=1/max(Re(eigen(W,only.values=TRUE)$values))
    }
  }else{
    rho1b <- 1/Re(eigs(W,1,which="SR",opts=list(retvec=FALSE))$values)
    if(identical(rho1b,complex(0))>0){
      rho1b=1/min(Re(eigen(W,only.values=TRUE)$values))
    }
  }
  if(sum(abs(rho)<abs(mf*rho1b))==R){inout=1}
  return(inout)
}
mh.rho.R.beta0 <- function(mu.cand,Sigma.cand,R,Wlist,Id,y,Wy,lndet.curr,s2.curr,
                           X.beta.tilde.curr,g,rss,val.curr,mu.prior,Sigma.prior,Ay.curr,
                           FUN) {
  inout <- 0
  while(inout<1){
    val.cand <- c(rmvnorm(1,mean=mu.cand,sigma=Sigma.cand))
    rho.cand <- val.cand[-(R+1)]
    beta0.cand <- val.cand[R+1]
    if(FUN(rho.cand,Wlist,R)>0){inout <- 1}
  }
  rho.cand.Wlist <- list()
  for(i in 1:R){
    rho.cand.Wlist[[i]] <- rho.cand[i]*Wlist[[i]]
  }
  A.cand <- Id-Reduce("+", rho.cand.Wlist)
  lndet.cand <- determinant(A.cand)$modulus[1]
  Ay.cand <- y-colSums((rho.cand*t(Wy)))
  if(log(runif(1))<
     lndet.cand-lndet.curr
     -(1/(2*s2.curr))*(sum(Ay.cand**2)-2*beta0.cand*sum(Ay.cand)-2*sum(Ay.cand*X.beta.tilde.curr)+
                       g*beta0.cand**2+2*beta0.cand*sum(X.beta.tilde.curr)+sum(X.beta.tilde.curr**2)
                       -rss) +
     dmvnorm(rho.cand,mean=mu.prior,sigma=Sigma.prior,log=T) -
     dmvnorm(val.curr[-(R+1)],mean=mu.prior,sigma=Sigma.prior,log=T) +
     dmvnorm(val.curr,mean=mu.cand,sigma=Sigma.cand,log=T) -
     dmvnorm(val.cand,mean=mu.cand,sigma=Sigma.cand,log=T)){
    val <- val.cand
    lndet <- lndet.cand
    Ay <- Ay.cand
  }else{
    val <- val.curr
    lndet <- lndet.curr
    Ay <- Ay.curr
  }
  parlist <- list(val,lndet,Ay)
  names(parlist) <- c("value","lndet","Ay")
  return(parlist)
}
# MCMC SAMPLERS for Bayesian NAM with single weight matrix
# flat prior
MCMC.F <- function(N=1e4, y, X, W, startval, burn=0) {
  g <- length(y)
  Id <- diag(g)
  ones <- rep(1,g)
  Wy <- c(W%*%y)
  yWones <- sum(Wy)
  yWWy <- sum(Wy**2)
  EW <- Re(eigen(W)$values)
  tau2 <- sum(EW**2)
  lb <- 1/min(EW)
  ub <- 1/max(EW)
  k <- ncol(X)
  if(k==1){
    Xtilde <- matrix(1,nrow=g,ncol=1)
    noeffects <- TRUE
  }else{
    Xtilde <- as.matrix(X[,-1])
    noeffects <- FALSE
  }
  XtXi <- solve(t(X)%*%X)
  XtXiXt <- XtXi%*%t(X)
  M <- Id-X%*%XtXiXt
  Sigmai12_help <- yWones
  Sigmai21_help <- Sigmai12_help
  Sigmai22_help <- g
  # Set starting values
  if(missing(startval)) {
    stop("set starting values")
  }else{
    rhoS <- startval[1]
    betaS <- startval[2:(1+k)]
    beta1S <- betaS[1]
    if(noeffects){
      betatildeS <- 0
    }else{
      betatildeS <- betaS[-1]
    }
    s2S <- startval[2+k]
    Ay <- y - rhoS*Wy
  }
  lndet <- sum(log(1 - rhoS*EW))
  XbetatildeS <- c(Xtilde%*%betatildeS)
  SS <- sum(Ay**2) - 2*beta1S*sum(Ay) - 2*sum(Ay*XbetatildeS) + g*beta1S**2 +
    2*beta1S*sum(XbetatildeS) + sum(XbetatildeS**2)
  # Create output vectors and matrices
  rho.draws <- rep(NA,N)
  beta.draws <- matrix(NA, nrow=N, ncol=k)
  sigma2.draws <- rep(NA,N)
  rho.draws[1] <- rhoS
  beta.draws[1,] <- betaS
  sigma2.draws[1] <- s2S

  fitted.vals <- res <- matrix(NA, nrow=N, ncol=length(y))
  fitted.vals[1,] <- c(solve(Id-rhoS*W)%*%(X%*%betaS+rnorm(g,mean=0,sd=sqrt(s2S))))
  res[1,] <- y - fitted.vals[1,]

  # Run the MCMC sampler
  pb = txtProgressBar(min = 1, max = N, initial = 1, width = 50, style = 3)

  for (r in 2:N) {
    # Draw (rho,beta1)
    Sigmai11_help <- yWWy + s2S*tau2
    Sigmai_help <- rbind(c(Sigmai11_help, Sigmai12_help), c(Sigmai21_help, Sigmai22_help))
    Sigma_help <- (1/(Sigmai11_help*Sigmai22_help - Sigmai12_help**2))*
      rbind(c(Sigmai22_help, -Sigmai12_help), c(-Sigmai21_help, Sigmai11_help))
    Sigma <- s2S*Sigma_help
    restilde <- y - c(Xtilde%*%betatildeS)
    SStilde <- sum(Wy*restilde)
    mu1 <- sum(restilde*yWones-SStilde)/(yWones**2 - g*Sigmai11_help)
    if(mu1 > ub | mu1 < lb){
      a <- (lb - mu1)/sqrt(Sigma[1,1])
      b <- (ub - mu1)/sqrt(Sigma[1, 1])
      z <- pnorm(b) - pnorm(a)
      mu1 <- mu1 + (sqrt(Sigma[1, 1])/z)*(dnorm(a) - dnorm(b))}
    mu2 <- (SStilde - mu1*Sigmai11_help)/yWones
    mu <- c(mu1, mu2)
    Sigma <- as.matrix(nearPD(Sigma)$mat)
    Draw <- mh.up.rhobeta1.F(mu, Sigma, lb, ub, y, Wy, EW, lndet, s2S, XbetatildeS, g, SS,
                             c(rhoS, beta1S), Ay)
    rhoS <- Draw$value[1]
    beta1S <- Draw$value[2]
    Ay <- Draw$Ay
    lndet <- Draw$lndet
    # Draw betatilde
    if(!noeffects){
      mu_beta <- c(XtXiXt%*%Ay)
      Sigma_beta <- s2S*XtXi
      Sigma_beta11 <- c(Sigma_beta[1:1])
      Sigma_beta12 <- c(Sigma_beta[1, 2:k])
      Sigma_beta22 <- Sigma_beta[2:k, 2:k]
      mu_betatilde <- mu_beta[-1] + Sigma_beta12*(beta1S - mu_beta[1])/Sigma_beta11
      Sigma_betatilde <- Sigma_beta22 - outer(Sigma_beta12, Sigma_beta12)/Sigma_beta11
      betatildeS <- c(rmvnorm(1, mean=mu_betatilde, sigma=Sigma_betatilde))
    }
    XbetatildeS <- c(Xtilde%*%betatildeS)
    SS <- sum(Ay**2) - 2*beta1S*sum(Ay) - 2*sum(Ay*XbetatildeS) + g*beta1S**2 +
      2*beta1S*sum(XbetatildeS) + sum(XbetatildeS**2)
    # Draw sigma2
    s2S <- 1/rgamma(1, shape=g/2, rate=SS/2)
    rho.draws[r] <- rhoS
    if(noeffects){
      beta.draws[r,] <- c(beta1S)
    }else{
      beta.draws[r,] <- c(beta1S, betatildeS)
    }
    sigma2.draws[r] <- s2S

    fitted.vals[r,] <- c(solve(Id-rhoS*W)%*%(X%*%beta.draws[r,]+rnorm(g,mean=0,sd=sqrt(s2S))))
    res[r,] <- y - fitted.vals[r,]

    setTxtProgressBar(pb,r)
  }
  message("\n")
  # Return the posterior draws
  outlist <- list(as.matrix(rho.draws[1:(N-burn)+burn]), as.matrix(beta.draws[1:(N-burn)+burn,]),
                  as.matrix(sigma2.draws[1:(N-burn)+burn]),fitted.vals[1:(N-burn)+burn,],
                  res[1:(N-burn)+burn,])
  names(outlist) <- c("rho.draws", "beta.draws", "sigma2.draws","fitted.vals","res")
  return(outlist)
}
# normal prior (default is the weakly informative prior from Dittrich, Leenders, Mulder (2016)
MCMC.N <- function(N=1e4, y, X, W, m, std2, startval, burn=0){

  g <- length(y)
  Id <- diag(g)
  ones <- rep(1,g)
  Wy <- c(W%*%y)
  yWones <- sum(Wy)
  yWWy <- sum(Wy**2)
  EW <- Re(eigen(W)$values)
  tau2 <- sum(EW**2)
  lb <- 1/min(EW)
  ub <- 1/max(EW)
  k <- ncol(X)
  if(k==1){
    Xtilde <- matrix(1,nrow=g,ncol=1)
    noeffects <- TRUE
  }else{
    Xtilde <- as.matrix(X[,-1])
    noeffects <- FALSE
  }
  XtXi <- solve(t(X)%*%X)
  XtXiXt <- XtXi%*%t(X)
  M <- Id - X%*%XtXiXt
  Sigmai12_help <- yWones
  Sigmai21_help <- Sigmai12_help
  Sigmai22_help <- g
  # Set starting values
  if(missing(startval)){
    scalepar=.999
    lb_ml <- scalepar/min(EW)
    ub_ml <- scalepar/max(EW)
    My <- M%*%y
    yMy <- sum(My**2)
    MWy <- M%*%Wy
    yMWy <- sum(y*MWy)
    yWMWy <- sum(MWy**2)
    rhoS <- 0 #optim(par=0, fn=loglik, gr=loglik.gr, method="L-BFGS-B",
    #EW=EW, g=g, yMy=yMy, yMWy=yMWy, yWMWy=yWMWy, lower=lb_ml, upper=ub_ml)$par
    Ay <- y-rhoS*Wy
    yAMAy <- yMy - 2*rhoS*yMWy + rhoS**2*yWMWy
    betaS <- c(XtXiXt%*%Ay)
    beta1S <- betaS[1]
    betatildeS <- betaS[-1]
    s2S <- yAMAy/g
  }else{
    rhoS <- startval[1]
    betaS <- startval[2:(1+k)]
    beta1S <- betaS[1]
    if(noeffects){
      betatildeS <- 0
    }else{
      betatildeS <- betaS[-1]
    }
    s2S <- startval[2+k]
    Ay <- y - rhoS*Wy
  }
  lndet <- sum(log(1 - rhoS*EW))
  XbetatildeS <- c(Xtilde%*%betatildeS)
  SS <- sum(Ay**2) - 2*beta1S*sum(Ay) - 2*sum(Ay*XbetatildeS) + g*beta1S**2 +
    2*beta1S*sum(XbetatildeS) + sum(XbetatildeS**2)
  # Create output vectors and matrices
  rho.draws <- NULL
  beta.draws <- matrix(NA, nrow=N, ncol=k)
  sigma2.draws <- NULL
  rho.draws[1] <- rhoS
  beta.draws[1,] <-betaS
  sigma2.draws[1] <- s2S

  fitted.vals <- res <- matrix(NA, nrow=N, ncol=length(y))
  fitted.vals[1,] <- c(solve(Id-rhoS*W)%*%(X%*%betaS+rnorm(g,mean=0,sd=sqrt(s2S))))
  res[1,] <- y - fitted.vals[1,]

  # Run the MCMC sampler
  pb = txtProgressBar(min = 1, max = N, initial = 1, width = 50, style = 3)
  for (r in 2:N){
    # Draw (rho,beta1)
    Sigmai11_help <- yWWy + s2S*tau2 + s2S/std2
    Sigmai_help <- rbind(c(Sigmai11_help, Sigmai12_help), c(Sigmai21_help, Sigmai22_help))
    Sigma_help <- (1/(Sigmai11_help*Sigmai22_help - Sigmai12_help**2))*
      rbind(c(Sigmai22_help, -Sigmai12_help), c(-Sigmai21_help, Sigmai11_help))
    Sigma <- s2S*Sigma_help
    restilde <- c(y - Xtilde%*%betatildeS)
    SStilde <- sum(Wy*restilde)
    mu1 <- sum(restilde*yWones - SStilde - m*s2S/std2)/(yWones**2 - g*Sigmai11_help)
    if(mu1 >ub | mu1 < lb){
      a <- (lb - mu1)/sqrt(Sigma[1, 1])
      b <- (ub - mu1)/sqrt(Sigma[1, 1])
      z <- pnorm(b) - pnorm(a)
      mu1 <- mu1 + (sqrt(Sigma[1, 1])/z)*(dnorm(a) - dnorm(b))}
    mu2 <- (SStilde + m*s2S/std2 - mu1*Sigmai11_help)/yWones
    mu <- c(mu1, mu2)
    Sigma <- as.matrix(nearPD(Sigma)$mat)
    Draw <- mh.up.rhobeta1.N(mu, Sigma, lb, ub, y, Wy, EW, lndet, s2S, XbetatildeS, g, SS,
                             cval=c(rhoS, beta1S), Ay, m, std2)
    rhoS <- Draw$value[1]
    beta1S <- Draw$value[2]
    Ay <- Draw$Ay
    lndet <- Draw$lndet
    # Draw betatilde
    mu_beta <- c(XtXiXt%*%Ay)
    Sigma_beta <- s2S*XtXi
    if(!noeffects){
      Sigma_beta11 <- c(Sigma_beta[1:1])
      Sigma_beta12 <- c(Sigma_beta[1, 2:k])
      Sigma_beta22 <- Sigma_beta[2:k, 2:k]
      mu_betatilde <- mu_beta[-1] + Sigma_beta12*(beta1S - mu_beta[1])/Sigma_beta11
      Sigma_betatilde <- Sigma_beta22 - outer(Sigma_beta12, Sigma_beta12)/Sigma_beta11
      betatildeS <- c(rmvnorm(1, mean=mu_betatilde, sigma=Sigma_betatilde))
    }
    XbetatildeS <- c(Xtilde%*%betatildeS)
    SS <- sum(Ay**2) - 2*beta1S*sum(Ay) - 2*sum(Ay*XbetatildeS) + g*beta1S**2 +
      2*beta1S*sum(XbetatildeS) + sum(XbetatildeS**2)
    # Draw sigma2
    s2S <- 1/rgamma(1, shape=g/2, rate=SS/2)
    rho.draws[r] <- rhoS
    if(!noeffects){
      beta.draws[r,] <- c(beta1S, betatildeS)
    }else{
      beta.draws[r,] <- c(beta1S)
    }
    sigma2.draws[r] <- s2S

    fitted.vals[r,] <- c(solve(Id-rhoS*W)%*%(X%*%beta.draws[r,]+rnorm(g,mean=0,sd=sqrt(s2S))))
    res[r,] <- y - fitted.vals[r,]

    setTxtProgressBar(pb,r)
  }
  message("\n")
  # Return the posterior draws
  outlist <- list(as.matrix(rho.draws[1:(N-burn)+burn]), as.matrix(beta.draws[1:(N-burn)+burn,]),
                  as.matrix(sigma2.draws[1:(N-burn)+burn]),fitted.vals[1:(N-burn)+burn,],
                  res[1:(N-burn)+burn,])
  names(outlist) <- c("rho.draws", "beta.draws", "sigma2.draws","fitted.vals","res")
  return(outlist)}
# independence Jeffreys prior
MCMC.IJ <- function(N=1e4, y, X, W, startval, burn=0){
  g <- length(y)
  Id <- diag(g)
  ones <- rep(1,g)
  Wy <- c(W%*%y)
  yWones <- sum(Wy)
  yWWy <- sum(Wy**2)
  EW <- Re(eigen(W)$values)
  tau2 <- sum(EW**2)
  lb <- 1/min(EW)
  ub <- 1/max(EW)
  k <- ncol(X)
  if(k==1){
    Xtilde <- matrix(1,nrow=g,ncol=1)
    noeffects <- TRUE
  }else{
    Xtilde <- as.matrix(X[,-1])
    noeffects <- FALSE
  }
  XtXi <- solve(t(X)%*%X)
  XtXiXt <- XtXi%*%t(X)
  Sigmai12_help <- yWones
  Sigmai21_help <- Sigmai12_help
  Sigmai22_help <- g
  #set starting values
  if(missing(startval)) {
    stop("starting values are needed")
  }else{
    rhoS <- startval[1]
    betaS <- startval[2:(1+k)]
    beta1S <- betaS[1]
    if(noeffects){
      betatildeS <- 0
    }else{
      betatildeS <- betaS[-1]
    }
    s2S <- startval[2+k]
  }
  A <- Id-rhoS*W
  Ay <-c(A%*%y)
  SS <- sum((Ay - c(X%*%betaS))**2)
  B <- W%*%solve(A)
  trBtB <- sum(B**2)
  trB2 <- sum((EW/(1 - rhoS*EW))**2)
  SSBXbeta <- sum(c(B%*%X%*%betaS)**2)
  trB <- tr(B)
  theta.samples <- matrix(NA, nrow=N, ncol=2+k)
  theta.samples[1, ] <- c(rhoS, betaS, s2S)

  fitted.vals <- res <- matrix(NA, nrow=N, ncol=length(y))
  fitted.vals[1,] <- c(solve(Id-rhoS*W)%*%(X%*%betaS+rnorm(g,mean=0,sd=sqrt(s2S))))
  res[1,] <- y - fitted.vals[1,]

  pb = txtProgressBar(min = 1, max = N, initial = 1, width = 50, style = 3)
  for (r in 2:N) {
    # Draw (rho,beta1)
    Sigmai11_help <- yWWy + s2S*tau2
    Sigmai_help <- rbind(c(Sigmai11_help, Sigmai12_help), c(Sigmai21_help, Sigmai22_help))
    Sigma_help <- (1/(Sigmai11_help*Sigmai22_help - Sigmai12_help**2))*
      rbind(c(Sigmai22_help, -Sigmai12_help), c(-Sigmai21_help, Sigmai11_help))
    Sigma <- s2S*Sigma_help
    restilde <- y-c(Xtilde%*%betatildeS)
    SStilde <- sum(Wy*restilde)
    mu1 <- sum(restilde*yWones - SStilde)/(yWones**2 - g*Sigmai11_help)
    if(mu1 > ub | mu1 < lb){
      a <- (lb - mu1)/sqrt(Sigma[1, 1])
      b <- (ub - mu1)/sqrt(Sigma[1, 1])
      z <- pnorm(b) - pnorm(a)
      mu1 <- mu1 + (sqrt(Sigma[1, 1])/z)*(dnorm(a) - dnorm(b))}
    mu2 <- (SStilde - mu1*Sigmai11_help)/yWones
    mu <- c(mu1, mu2)
    Sigma <- as.matrix(nearPD(Sigma)$mat)
    Draw <- mh.up.rhobeta1.IJ(mu, Sigma, lb, ub, Id, W, betatildeS, X, EW, c(rhoS, beta1S),
                              s2S, y, SS, g, trBtB, trB2, SSBXbeta, trB)
    rhoS <- Draw[1]
    beta1S <- Draw[-1]
    betaS[1] <- beta1S
    A <- Id - rhoS*W
    Ay <- c(A%*%y)
    SS <- sum((Ay - c(X%*%betaS))**2)
    B <- W%*%solve(A)
    trBtB <- sum(B**2)
    trB2 <- sum((EW/(1 - rhoS*EW))**2)
    SSBXbeta <- sum(c(B%*%X%*%betaS)**2)
    trB <- tr(B)
    # Draw Beta_rest
    if(!noeffects){
      betatildeS <- mh.up.betatilde.IJ(XtXiXt, Ay, s2S, XtXi, k, beta1S, B, X, SS, EW, g,
                                       trBtB, trB2, SSBXbeta, trB, betatildeS)
      betaS <- c(beta1S, betatildeS)
    }
    SS <- sum((Ay - c(X%*%betaS))**2)
    SSBXbeta <- sum(c(B%*%X%*%betaS)**2)
    # Draw Sigma2
    s2S <- mh.up.s2.IJ(g, SS, s2S, trBtB, trB2, SSBXbeta, trB)
    theta.samples[r, ] <- c(rhoS, betaS, s2S)

    fitted.vals[r,] <- c(solve(Id-rhoS*W)%*%(X%*%betaS+rnorm(g,mean=0,sd=sqrt(s2S))))
    res[r,] <- y - fitted.vals[r,]

    setTxtProgressBar(pb,r)
  }
  message("\n")
  return(list(rho.draws=as.matrix(theta.samples[1:(N-burn)+burn,1]),
              beta.draws=as.matrix(theta.samples[1:(N-burn)+burn,1+1:k]),
              sigma2.draws=as.matrix(theta.samples[1:(N-burn)+burn,k+2]),
              fitted.vals=fitted.vals[1:(N-burn)+burn,],res=res[1:(N-burn)+burn,]))
}
# MCMC SAMPLER for Bayesian NAM with multiple weight matrices
MCMC.multiple <- function(y,X,Wlist,mu.prior,Sigma.prior,N,burnin){

  if(is.null(mu.prior)){
    mu.prior <- rep(0,length(Wlist))
  }
  if(is.null(Sigma.prior)){
    Sigma.prior <- diag(length(Wlist))*100
  }
  g <- length(y)
  R <- length(Wlist)

  if(min(sapply(Wlist,is.numeric))<1){stop("The connectivity matrices W must be numeric.")}
  if(length(unique(g,sapply(Wlist,nrow),sapply(Wlist,ncol)))>1){stop("Order(s) of W must match length of y.")}
  if(!is.vector(mu.prior)){stop("The prior mean must be a vector.")}
  if(!is.numeric(mu.prior)){stop("The prior mean must be numeric.")}
  if(!is.matrix(Sigma.prior)){stop("The prior covariance matrix must be a matrix.")}
  if(!is.numeric(Sigma.prior)){stop("The prior covariance matrix must be numeric.")}
  if(abs(length(mu.prior)-R)>.5){stop("Length of prior mean must match with the number of connectivity matrices W.")}
  if(length(unique(R,nrow(Sigma.prior),ncol(Sigma.prior)))>1){stop("Dimensions of prior covariance matrix  must match
                                                                   number of connectivity matrices W.")}
  if(!is.positive.definite(Sigma.prior)){stop("The prior covariance matrix prior covariance matrix must be
                                              positive-definite.")}
  if(nrow(Sigma.prior)!=length(mu.prior)){stop("The dimensions of the prior covariance matrix should match with
                                               the number of connectivity matrices W")}
  if(R<3){
    if(max(sapply(Wlist,rowSums))==1){
      in.space <- in.space.2.rsd
    }else{
      in.space <- in.space.2
    }
  }else{
    in.space <- in.space.R
  }

  Id <- diag(g)
  ones <- rep(1,g)
  k <- ncol(X)
  if(k==1){
    X.tilde <- matrix(1,nrow=g,ncol=1)
    noeffects <- TRUE
  }else{
    X.tilde <- as.matrix(X[,-1])
    noeffects <- FALSE
  }
  XtXi <- solve(t(X)%*%X)
  XtXiXt <- XtXi%*%t(X)
  M <- Id-X%*%XtXiXt
  yMy <- sum(c(M%*%y)**2)
  Wy <- matrix(NA,g,R)
  tr.WW <- matrix(NA,R,R)
  for(i in 1:R){
    Wy[,i] <- c(Wlist[[i]]%*%y)
  }
  sumWy <- colSums(Wy)
  yWWy <- t(Wy)%*%Wy
  for(i in 1:R){
    for(j in i:R){
      tr.WW[i,j] <- sum(t(Wlist[[i]])*Wlist[[j]])
      tr.WW[j,i] <- tr.WW[i,j]
    }
  }
  Sigmai.prior <- solve(Sigma.prior)
  Sigmaih1 <- rbind(cbind(yWWy,sumWy),c(sumWy,g))
  Sigmaih2 <- rbind(cbind(tr.WW+Sigmai.prior,rep(0,R)),rep(0,R+1))

  vals.start <- lm(y ~ -1 + X)
  rho.curr <- rep(0,R)
  beta0.curr <- vals.start$coefficients[1]
  if(k==1){
    beta.tilde.curr <- 0
  }else{
    beta.tilde.curr <- vals.start$coefficients[-1]
  }
  s2.curr <- sigma(vals.start)**2

  rho.curr.Wlist <- list()
  for(j in 1:R){
    rho.curr.Wlist[[j]] <- rho.curr[j]*Wlist[[j]]
  }
  A.curr <- Id - Reduce("+",rho.curr.Wlist)
  lndet.curr <- determinant(A.curr)$modulus[1]
  Ay.curr <- y-colSums((rho.curr*t(Wy)))
  X.beta.tilde.curr <- c(X.tilde%*%beta.tilde.curr)
  rss <- sum(Ay.curr**2)-2*beta0.curr*sum(Ay.curr)-2*sum(Ay.curr*X.beta.tilde.curr)+
    g*beta0.curr**2+2*beta0.curr*sum(X.beta.tilde.curr)+sum(X.beta.tilde.curr**2)

  rho.samples <- matrix(NA,nrow=burnin+N,ncol=R)
  beta.samples <- matrix(NA,nrow=burnin+N,ncol=k)
  sigma.samples <- rep(NA,burnin+N)
  fitted.vals <- res <- matrix(NA,nrow=burnin+N,ncol=length(y))
  rho.samples[1,] <- rho.curr
  if(noeffects){
    beta.samples[1,] <- c(beta0.curr)
  }else{
    beta.samples[1,] <- c(beta0.curr,beta.tilde.curr)
  }
  sigma.samples[1] <- sqrt(s2.curr)
  fitted.vals[1,] <- c(solve(A.curr)%*%(X.beta.tilde.curr+rnorm(g,mean=0,sd=sqrt(s2.curr))))
  res[1,] <- y - fitted.vals[1,]

  pb = txtProgressBar(min = 1, max = burnin+N, initial = 1, width = 50, style = 3)
  for (i in 2:(burnin+N)){
    # Draw rho
    Sigmai.cand <- Sigmaih1/s2.curr+Sigmaih2
    if(is.positive.definite(Sigmai.cand)<1){
      Sigmai.cand <- nearPD(Sigmai.cand)$mat
    }
    Sigma.cand <- solve(Sigmai.cand)
    r.tilde <- y-c(X.tilde%*%beta.tilde.curr)
    z.rho <- c(t(Wy)%*%r.tilde)/s2.curr+c(Sigmai.prior%*%mu.prior)
    z.beta0 <- sum(r.tilde)/s2.curr
    z <- c(z.rho,z.beta0)
    mu.cand <- c(Sigma.cand%*%z)
    draw <- mh.rho.R.beta0(mu.cand,Sigma.cand,R,Wlist,Id,y,Wy,
                           lndet.curr,s2.curr,X.beta.tilde.curr,g,
                           rss,c(rho.curr,beta0.curr),mu.prior,
                           Sigma.prior,Ay.curr,FUN=in.space)
    rho.curr <- draw$value[-(R+1)]
    beta0.curr <- draw$value[R+1]
    for(j in 1:R) {
      rho.curr.Wlist[[j]] <- rho.curr[j]*Wlist[[j]]
    }
    A.curr <- Id - Reduce("+",rho.curr.Wlist)
    lndet.curr <- draw$lndet
    Ay.curr <- draw$Ay
    rho.samples[i,] <- rho.curr
    beta.samples[i,1] <- beta0.curr
    # Draw beta.tilde
    mu.beta <- c(XtXiXt%*%Ay.curr)
    Sigma.beta <- s2.curr*XtXi
    if(k > 1){
      Sigma.beta11 <- Sigma.beta[1:1]
      Sigma.beta12 <- Sigma.beta[1,2:k]
      Sigma.beta22 <- Sigma.beta[2:k,2:k]
      mu.beta.tilde <- mu.beta[-1]+Sigma.beta12*(beta0.curr-mu.beta[1])/Sigma.beta11
      Sigma.beta.tilde <- Sigma.beta22-outer(Sigma.beta12,Sigma.beta12)/Sigma.beta11
      beta.tilde.curr <- c(rmvnorm(1,mean=mu.beta.tilde,sigma=Sigma.beta.tilde))
      X.beta.tilde.curr <- c(X.tilde%*%beta.tilde.curr)
    }
    rss <- sum(Ay.curr**2)-2*beta0.curr*sum(Ay.curr)-2*sum(Ay.curr*X.beta.tilde.curr)+
      g*beta0.curr**2+2*beta0.curr*sum(X.beta.tilde.curr)+sum(X.beta.tilde.curr**2)
    beta.samples[i,-1] <- beta.tilde.curr
    # Draw Sigma2
    s2.curr <- 1/rgamma(1,shape=g/2,rate=rss/2)
    sigma.samples[i] <- sqrt(s2.curr)
    fitted.vals[i,] <- c(solve(A.curr)%*%(c(X%*%beta.samples[i,]) +
                                            rnorm(g,mean=0,sd=sqrt(s2.curr))))
    res[i,] <- y-fitted.vals[i,]

    setTxtProgressBar(pb,i)
  }
  message("\n")

  if(burnin > 0){
    rho.samples <- as.matrix(rho.samples[-(1:burnin),])
    beta.samples <- as.matrix(beta.samples[-(1:burnin),])
    sigma.samples <- as.matrix(sigma.samples[-(1:burnin)])
    fitted.vals <- fitted.vals[-(1:burnin),]
    res <- res[-(1:burnin),]
  }

  if(is.null(colnames(X))>0){
    colnames(beta.samples)=paste("X",1:k,sep="")
  }else{
    colnames(beta.samples)=colnames(X)
  }
  colnames(rho.samples) <- paste("rho",1:R,sep="")
  theta.samples <- cbind(rho.samples,beta.samples,sigma.samples)
  colnames(theta.samples) <- c(colnames(rho.samples),colnames(beta.samples),"sigma")

  out1 <- list()
  out1$beta.draws <- beta.samples
  out1$rho.draws <- rho.samples
  out1$sigma2.draws <- sigma.samples
  out1$fitted.vals <- fitted.vals
  out1$res <- res
  out1$priorsigma <- Sigma.prior
  out1$priormean <- mu.prior

  return(out1)

}

#' @method print banam
#' @export
print.banam <- function(x, digits = max(3, getOption("digits") - 4), ...){
  coefs.mean <- c(x$rho.mean,x$beta.mean)
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("Bayesian estimation of the network autocorrelation model \n")
  cat(paste0("Prior: ",x$prior),"\n", sep = "")
  cat("\n")
  cat("Parameters:\n")
  print.default(format(coefs.mean, digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
}

#' @method summary banam
#' @export
summary.banam <- function(object, ...){
  cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
  cat("Bayesian estimation of the network autocorrelation model \n")
  cat(paste0("Prior: ",object$prior),"\n", sep = "")
  cat("\n")
  if(ncol(object$rho.draws)==1){
    cat("Network autocorrelation:\n")
    results1 <- t(object$summarystats[1,])
    row.names(results1) <- "rho"
  }else{
    cat("Network autocorrelations:\n")
  }
  numrho <- ncol(object$rho.draws)
  if(numrho==1){
    results1 <- t(object$summarystats[1,])
    row.names(results1) <- "rho"
  }else{
    results1 <- object$summarystats[1:numrho,]
    row.names(results1) <- colnames(object$rho.draws)
  }
  print.default(format(results1, digits = 4), print.gap = 2, quote = FALSE)

  cat("\n")
  cat("Coefficients:")
  cat("\n")
  numbeta <- ncol(object$beta.draws)
  if(numbeta==1){
    results2 <- object$summarystats[numrho+1,]
    row.names(results2) <- "(Intercept)"
  }else{
    results2 <- object$summarystats[numrho+1:numbeta,]
    row.names(results2) <- names(object$beta.mean)
  }
  print.default(format(results2, digits = 4), print.gap = 2, quote = FALSE)
  results3 <- t(object$summarystats[numrho+numbeta+1,])
  row.names(results3) <- "sigma2"
  cat("\n")
  cat("Residual variance:")
  cat("\n")
  print.default(format(results3, digits = 4), print.gap = 2, quote = FALSE)
}



#' @method get_estimates banam
#' @export
get_estimates.banam <- function(x, ...){
  out <- list()
  out$estimate <- c(x$rho.mean,x$beta.mean)
  out$Sigma <- list(cov(cbind(x$rho.draws,x$beta.draws)))
  class(out) <- "model_estimates"
  attr(out, "analysisType") <- "banam"
  out
}



