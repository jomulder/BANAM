
## computes the logarithm of the marginal likelihood ## under a precise hypothesis rho=c and data y, X, W
lnmarglik.p <- function(y, X, W, rho0) {
  g <- length(y)
  k <- ncol(X)
  gminusk2 <- .5*(g - k)
  EW <- eigen(W,only.values=TRUE)$values
  XtX <- t(X) %*% X
  M <- diag(g) - X %*% solve(XtX) %*% t(X)
  yMy <- sum(c(M %*% y)**2)
  Wy <- W %*% y
  MWy <- c(M %*% Wy)
  yMWy <- sum(y*MWy)
  yWMWy <- sum(MWy**2)
  yAMAy <- c(yMy - 2*rho0*yMWy + rho0**2*yWMWy)
  lnml <- (- gminusk2*log(pi) + lgamma(gminusk2) - .5*log(det(XtX))
           + Re(sum(log(1 - rho0*EW))) - gminusk2*log(yAMAy))
  return(lnml)
}
## computes the logarithm of the marginal likelihood under an ## interval hypothesis a1<rho<a2 for a normal prior (mean, sd) ## for rho, data y, X, W, and R grid points
lnmarglik.n <- function(gridpoints=1e3, y, X, W, mean, sd, a1, a2) {
  g <- length(y)
  k <- ncol(X)
  gminusk2 <- .5*(g - k)
  EW <- eigen(W,only.values=TRUE)$values
  XtX <- t(X) %*% X
  M <- diag(g) - X %*% solve(XtX) %*% t(X)
  yMy <- sum(c(M %*% y)**2)
  Wy <- W %*% y
  MWy <- c(M %*% Wy)
  yMWy <- sum(y*MWy)
  yWMWy <- sum(MWy**2)
  scalefac <- .999
  rhoseqh <- seq(scalefac*a1, scalefac*a2, length=100)
  yAMAyh <- c(yMy - 2*rhoseqh*yMWy + rhoseqh**2*yWMWy)
  lognormdensh <- dnorm(rhoseqh, mean=mean, sd=sd, log=T)
  rhoseq <- seq(scalefac*a1, scalefac*a2, length=2*gridpoints - 1)
  yAMAy <- c(yMy - 2*rhoseq*yMWy + rhoseq**2*yWMWy)
  spread <- (rhoseq[3] - rhoseq[1])/6
  weights <- c(1, rep(c(4, 2), .5*(2*gridpoints - 4)), c(4, 1))
  normc <- pnorm(a2, mean=mean, sd=sd) - pnorm(a1, mean=mean, sd=sd)
  lognormdens <- dnorm(rhoseq, mean=mean, sd=sd, log=TRUE)
  logint <- unlist(lapply(1:(2*gridpoints - 1),function(r){
    Re(sum(log(1 - rhoseq[r]*EW))) -
      gminusk2*log(yAMAy[r]) + lognormdens[r]
  }))
  lnml <- - gminusk2*log(pi) + lgamma(gminusk2) - .5*log(det(XtX)) - log(normc) +
           log(spread) + log(sum(weights*exp(logint - max(logint)))) + max(logint)
  return(list(logmarglike=lnml,priorprob=normc))
}

#for multiple rho's
marglikeNAM <- function(S,y,X,Wlist,mu_prior=NULL,Sigma_prior=NULL,orderconstraints=NULL,
                        in.space=in.space.2){
  g <- length(y)
  N <- length(Wlist)
  Id <- diag(g)
  k <- ncol(X)
  gminusk2 <- .5*(g - k)
  XtX <- t(X) %*% X
  M <- Id - X %*% solve(t(X) %*% X) %*% t(X)
  yMy <- sum(c(M %*% y)**2)
  Wy <- matrix(NA, g, N)
  MWy <- matrix(NA, g, N)
  yMWy <- NULL
  yWMWy <- matrix(NA, N, N)
  Wy <- sapply(1:N, function(x) Wlist[[x]] %*% y)
  MWy <- M%*%Wy
  yMWy <- c(y%*%MWy)
  yWMWy <- t(MWy)%*%MWy

  #normal approximation of det(A_rho) (I think...)
  mu1 <- rep(0,N)
  W1 <- unlist(lapply(1:N, function(x) rep(list(t(Wlist[[x]])),N)))
  W2 <- unlist(rep(Wlist,N))
  vec <- W1*W2
  l <- length(vec)/N**2
  veclist <- lapply(1:N**2, function(x) vec[((x-1)*l+1):(x*l)])
  Sigma1i <- matrix(sapply(1:N**2, function(x) sum(veclist[[x]])),N,N)
  Sigma1i <- as.matrix(Matrix::nearPD(Sigma1i)$mat)

  #normal approximation for y'A'MAy part
  mu2 <- c(solve(yWMWy)%*%yMWy)
  yAMAy <- sum(outer(mu2, mu2)*yWMWy) - 2*sum(mu2*yMWy) + yMy
  Sigma2i <- 2*yWMWy/yAMAy*gminusk2

  #combine normal approximations in integrand of likelihood part
  Sigma_int <- solve(Sigma1i+Sigma2i)
  mu_int <- c(Sigma_int%*%(c(Sigma1i%*%mu1)+c(Sigma2i%*%mu2)))

  #combine normal parts of likelihood and prior
  #the caovariance matrix is multiplied with a factor of 2 to make the importance sampler
  #a bit more spread out than the integrand
  SigmaIS <- solve(Sigma1i+Sigma2i+solve(Sigma_prior)) * 2
  muIS <- as.vector(SigmaIS%*%(c((Sigma1i+Sigma2i)%*%mu_int)+c(solve(Sigma_prior)%*%mu_prior)))

  draws <- mvtnorm::rmvnorm(S,mean=muIS,sigma=SigmaIS)
  ISin <- unlist(lapply(1:S,function(s){
    in.space(rho=draws[s,],Wlist=Wlist,R=length(Wlist))
  }))
  IScounter <- sum(ISin)
  drawsin <- draws[ISin==1,]

  yAMAyseq <- sapply(1:IScounter, function(s){
    sum( drawsin[s,]%*%t(drawsin[s,])*yWMWy) - 2*sum(drawsin[s,]*yMWy) + yMy
  })
  logyAMAyseq <- log(yAMAyseq)
  logdensIS <- dmvnorm(drawsin,mean=muIS,sigma=SigmaIS,log=TRUE) - log(mean(ISin)) #Imp sampler is truncated normal
  logpriordens <- dmvnorm(drawsin,mean=mu_prior,sigma=Sigma_prior,log=TRUE) #normalizing constant due to truncated added later

  logint <- unlist(lapply(1:IScounter,function(s){
    rhoW <- Reduce("+",lapply(1:N, function(x) Wlist[[x]] * drawsin[s,x]))
    logpriordens[s] + log(abs(det(Id-rhoW))) - gminusk2*logyAMAyseq[s] - logdensIS[s]
  }))
  lnml <- log(mean(exp(logint - max(logint)))) + max(logint) + lgamma(gminusk2) -
    .5*log(det(XtX)) - gminusk2*log(pi)

  if(!is.null(orderconstraints)){
    #then order constraints are present, and thus we need to compute the posterior and prior probabilities
    postdraws <- MCMC.multiple(y,X,Wlist=Wlist,mu.prior=mu_prior,Sigma.prior=Sigma_prior,
                               N=3e3,burnin=1e3)$rho.draws
    if(nrow(orderconstraints)==1){
      RO <- t(orderconstraints[,-(N+1)])
      rO <- orderconstraints[,N+1]
    }else{
      RO <- orderconstraints[,-(N+1)]
      rO <- orderconstraints[,N+1]
    }
    #posterior probability of order constraints conditional on allowed region based on unconstrained NAM
    postprobCondOrder <- mean(apply((postdraws%*%t(RO) - rep(1,nrow(postdraws))%*%t(rO))>0,1,prod))

    # get normalizing constant for prior
    priordrawsMax <- 1e4
    storePriordraws <- rmvnorm(priordrawsMax,mean=mu_prior,sigma=Sigma_prior)
    checkConstraints <- apply(storePriordraws%*%t(RO) - rep(1,priordrawsMax)%*%t(rO) > 0,1,prod)==1
    checkNAM <- unlist(lapply(1:priordrawsMax,function(s){
      in.space(storePriordraws[s,],Wlist,R=length(Wlist))
    })) == 1
    #prior probability of order constraints and allowed region under unconstrained prior
    priorprob <- mean(checkConstraints*checkNAM)
    #prior probability of order constraints conditional on allowed region based on unconstrained NAM
    priorDrawsIn <- storePriordraws[checkNAM,]
    priorprobCondOrder <- mean(apply(priorDrawsIn%*%t(RO) - rep(1,sum(checkNAM))%*%t(rO) > 0,1,prod)==1)
  }else{
    #prior and post probability of order constraints conditional on allowed region based on unconstrained NAM
    postprobCondOrder <- 1
    priorprobCondOrder <- 1
    # get normalizing constant for prior
    priordrawsMax <- 1e4
    storePriordraws <- rmvnorm(priordrawsMax,mean=mu_prior,sigma=Sigma_prior)
    checkNAM <- unlist(lapply(1:priordrawsMax,function(s){
      in.space(storePriordraws[s,],Wlist,R=length(Wlist))
    })) == 1
    priorprob <- mean(checkNAM)
    priorDrawsIn <- storePriordraws[checkNAM,]
  }

  return(list(logmarglike=lnml - log(priorprob) + log(postprobCondOrder), priorprob=priorprobCondOrder,
              postprob=postprobCondOrder, priordraws=priorDrawsIn, priorprobAll=priorprob))
}



#' @method BF banam
#' @export
BF.banam <- function(x,
                  hypothesis = NULL,
                  prior.hyp = NULL,
                  complement = TRUE,
#                  priortype = NULL,
                  ...){

  #check if constrained hypotheses are formulated that can be tested
  if(!is.null(hypothesis)){

    betaHypo <- grepl("beta", hypothesis)
    rhoHypo <- grepl("rho", hypothesis)
    if(betaHypo & rhoHypo){
      stop(paste0("Constrained hypotheses can only be tested on either ",rhonames," or on ",
                  paste(names(x$beta.mean)[-1],collapse = ', ') ))
    }
  }

  #extract information from object
  y <- x$y
  X <- x$X
  Wlist <- x$W
  numrho <- length(Wlist)
  #rhodraws <- x$rho.draws

  # exploratory testing of rho
  #default normal prior
  prior_mean <- rep(0,numrho)
  #scale default prior scale to upper bound each network autocorrelation
  UBs <- unlist(lapply(1:numrho,function(r){
    EW_r <- Re(eigen(Wlist[[r]])$values)
    1/max(EW_r)
  }))
  if(numrho==1){
    prior_covm <- diag(1) * (.5*UBs)**2
  }else{
    prior_covm <- diag((.5*UBs)**2)
  }
  # }else{
  #   if(!is.list(priortype)){stop("For a 'banam' object 'priortype' must be a list of length 2 with (1) a vector of prior means and (2) a prior covariance matrix.")}
  #   if(length(priortype)!=2){stop("For a 'banam' object 'priortype' must be a list of length 2 with (1) a vector of prior means and (2) a prior covariance matrix.")}
  #   if(length(priortype[[1]])!=numrho){stop("For a 'banam' object the first element of 'priortype' must be a vector of prior means of length equal to the number of weight matrices.")}
  #   if(!is.matrix(priortype)){stop("For a 'banam' object 'priortype' must be a list of length 2 with (1) a vector of prior means and (2) a prior covariance matrix.")}
  #   if(sum(dim(length(priortype[[2]]))==c(numrho,numrho))!=2){stop("For a 'banam' object the first element of 'priortype' must be a vector of prior means of length equal to the number of weight matrices.")}
  #   prior_mean <- priortype[[1]]
  #   prior_covm <- priortype[[2]]
  # }
  # #OF via banam object
  # prior_mean <- x$priormean
  # prior_covm <- x$priorsigma

  #Exploratory BF testing of separate rho(s)
  if(numrho==1){
    #univariate NAM
    rhonames <- "rho"

    message("BANAM: Computing the marginal likelihood under the unconstrained NAM \n")

    marglikeH0 <- lnmarglik.p(y, X, W=Wlist[[1]], rho0=0)
    EW <- Re(eigen(Wlist[[1]],only.values=TRUE)$values)
    lb <- 1/min(EW)
    ub <- 1/max(EW)
    marglikeHu <- lnmarglik.n(gridpoints=1e3, y, X, W=Wlist[[1]], mean=prior_mean,
                              sd=sqrt(prior_covm[1,1]), a1 = lb, a2 = ub)

    #unconstrained posterior draws
    message("       Unconstrained posterior sampling... \n")
    #get posterior probbilities for positive and negative rho's
    starting <- c(x$rho.mean,x$beta.mean,x$sigma2.mean)
    postdrawsHu <- MCMC.N(N=5e3, y, X, W=Wlist[[1]], m=prior_mean, std2=sqrt(prior_covm[1,1]),
                          startval=starting, burn=1e3)
    priordrawsHu <- as.matrix(rtruncnorm(1e4,a=lb,b=ub,mean=prior_mean,sd=sqrt(prior_covm[1,1])))
    message("       ...finished \n")

    postNegative <- mean(postdrawsHu$rho.draws<0)
    postPositive <- mean(postdrawsHu$rho.draws>0)

    BFtu_exploratory <- t(c(marglikeH0 - marglikeHu$logmarglike,log(postNegative),log(postPositive)))

  }else{
    #multivariate NAM
    if(numrho<3){
      if(max(sapply(Wlist,rowSums))==1){
        inspaceUnc <- in.space.2.rsd
      }else{
        inspaceUnc <- in.space.2
      }
    }else{
      inspaceUnc <- in.space.R
    }

    #marginal likelihood under unconstrained model
    message("BANAM: Computing the marginal likelihood under the unconstrained NAM \n")

    marglikeHu <- marglikeNAM(S=1e3,y,X,Wlist,mu_prior=prior_mean,Sigma_prior=prior_covm,
                              in.space=inspaceUnc)
    priordrawsHu <- marglikeHu$priordraws
    #unconstrained posterior draws
    message("       Unconstrained posterior sampling... \n")
    #get posterior probbilities for positive and negative rho's
    postdrawsAllHu <- MCMC.multiple(y,X,Wlist,mu.prior=prior_mean,Sigma.prior=prior_covm,
                                    N=5e3,burnin=1e3)

    postdrawsHu <- postdrawsAllHu$rho.draws
    postPositive <- apply(postdrawsHu>0,2,mean)
    postNegative <- apply(postdrawsHu<0,2,mean)
    message("       ...finished \n")

    rhonames <- paste0(rep("rho",numrho),1:numrho)

    BFtu_exploratory <- matrix(unlist(lapply(1:numrho,function(j){
      if(numrho==2){
        Wdummy <- Wlist[-j]
        EW <- Re(eigen(Wdummy[[1]])$values)
        lb <- 1/min(EW)
        ub <- 1/max(EW)
        marglike_explo <- rep(0,3)
        marglike_j <- lnmarglik.n(gridpoints=1e3, y, X, W=Wdummy[[1]], mean=prior_mean[j],
                                  sd=sqrt(prior_covm[j,j]), a1=lb, a2=ub)
        marglike_explo[1] <- marglike_j$logmarglike - marglikeHu$logmarglike
        marglike_explo[2] <- log(postNegative[j])
        marglike_explo[3] <- log(postPositive[j])
      }else{
        Wdummy <- Wlist[-j]
        if(numrho - 1 < 3){
          if(max(sapply(Wdummy,rowSums))==1){
            inSpaceDummy <- in.space.2.rsd
          }else{
            inSpaceDummy <- in.space.2
          }
        }else{
          inSpaceDummy <- in.space.R
        }
        marglike_explo <- rep(0,3)
        marglike_explo[1] <- marglikeNAM(S=1e3,y,X,Wdummy,mu_prior=prior_mean[-j],Sigma_prior=prior_covm[-j,-j],
                                         in.space=inSpaceDummy)[[1]] - marglikeHu[[1]]
        marglike_explo[2] <- log(postNegative[j])
        marglike_explo[3] <- log(postPositive[j])
      }
      return(marglike_explo)
    })),ncol=3,byrow=TRUE)

  }

  # default BF tests on beta using the AAFBF
  # #compute effective sample size for beta's using estimated rho
  # rhoWsum <- Reduce("+",
  #   lapply(1:length(Wlist),function(w){
  #     x$rho.mean[w] * Wlist[[w]]
  #   })
  # )
  # A_rho <- solve(diag(length(y)) - rhoWsum)
  # # effect sample size based on work of Faes, Molenberghs, Aerts, Verbeke, & Kenward (2009)
  # # and Bayarri et al. (2014).
  # tess <- sum(solve(cov2cor(A_rho %*% t(A_rho))))
  get_est <- get_estimates(x)
  Args <- list()
  Args$x <- get_est$estimate[-(1:numrho)]
  Args$Sigma <- get_est$Sigma[[1]][-(1:numrho),-(1:numrho)]
  Args$n <- length(y)
  Args$hypothesis <- NULL
  Args$prior.hyp <- NULL
  Args$complement <- NULL
  out <- do.call(BF, Args)

  BFtu_exploratory <- exp(BFtu_exploratory)
  BFtu_exploratory <- rbind(BFtu_exploratory,out$BFtu_exploratory)
  colnames(BFtu_exploratory) <- c("Pr(=0)","Pr(<0)","Pr(>0)")
  row.names(BFtu_exploratory) <- names(get_est$estimate)
  PHP_exploratory <- round(BFtu_exploratory / apply(BFtu_exploratory,1,sum),3)

  postestimates <- x$summarystats

  if(!is.null(hypothesis)){

    if(rhoHypo == TRUE){
      # test hypotheses on rho's
      parse_hyp <- parse_hypothesis(rhonames,hypothesis)
      parse_hyp$hyp_mat <- do.call(rbind, parse_hyp$hyp_mat)
      RrList <- make_RrList2(parse_hyp)
      RrE <- RrList[[1]]
      RrO <- RrList[[2]]

      # ADD CHECK WHICH TYPES OF HYPOTHESES CAN BE TESTED. FOR EXAMPLE THERE IS NO FUNCTION YET FOR
      # H: rho1=.2,rho2>0
      # check if rho's are only tested against each other or against zero
      numhyp <- length(RrE)
      for(h in 1:numhyp){
        if(!is.null(RrE[[h]])){
          for(r in 1:nrow(RrE[[h]])){
            row1 <- RrE[[h]][r,]
            if(numrho > 1 & !(sum(abs(row1))==1 || sum(row1)==0) ){
              stop("rho's can only be compared with each other or to zero.")
            }
          }
        }
      }

      output_marglike <- matrix(unlist(lapply(1:numhyp, function(h){

        message(paste0("Evaluate H",as.character(h),": ",parse_hyp$original_hypothesis[h],"; \n"))

        RrE_h <- RrE[[h]]
        RrO_h <- RrO[[h]]

        if(numrho == 1){

          if(is.null(RrE_h)){
            #only one-sided constraint(s)
            #compute posterior probability of constraints

            #add lower and upper bound
            RO_h <- as.matrix(RrO_h[,1])
            rO_h <- RrO_h[,numrho+1]
            RO_h <- rbind(RO_h,as.matrix(c(1,-1)))
            rO_h <- c(rO_h,lb,-ub)

            priordraws <- rnorm(1e5,mean=prior_mean,sd=sqrt(prior_covm[1,1]))
            priorProb <- mean(apply((priordraws%*%t(RO_h) - rep(1,length(priordraws))%*%t(rO_h))>0,1,prod))
            rm(priordraws)
            postProb <- mean(apply((postdrawsHu$rho.draws%*%t(RO_h) - rep(1,nrow(postdrawsHu$rho.draws))
                                    %*%t(rO_h))>0,1,prod))

            margLikeHt <- marglikeHu$logmarglike + log(marglikeHu$priorprob) - log(priorProb) + log(postProb)
          }
          if(is.null(RrO_h)){
            #only equality constraint
            rho0 <- RrE_h[1,2] / RrE_h[1,1]
            margLikeHt <- lnmarglik.p(y, X, Wlist[[1]], rho0)
            priorProb <- postProb <- 0
          }

          return(c(unlist(margLikeHt),postProb,priorProb,ifelse(is.null(RrE[[h]]),1,0)))

        }else{
          # code equal rho's with same integer for marglike2_Hq function
          unique_h <- 1:numrho
          if(!is.null(RrE[[h]])){
            unique_h <- rep(NA,numrho)
            zeroMat <- RrE[[h]][which(apply(RrE[[h]],1,sum)==1),]
            if(!is.matrix(zeroMat)){
              zeroMat <- matrix(zeroMat,nrow=1)
            }
            unique_h[which(apply(zeroMat,2,sum)==1)] <- 0
            teller <- 0
            rowEq <- which(apply(RrE[[h]],1,sum)==0)
            rowEqTeller <- 0
            while(length(rowEq)>0){
              #check which rows already have a parameter that is assigned in unique_h
              checkNA <- rep(NA,nrow(RrE[[h]]))
              for(row0 in rowEq){
                checkNA[row0] <- sum(is.na(unique_h[which(RrE[[h]][row0,1:numrho]!=0)]))
              }
              row1 <- which(checkNA==1)[1]
              row1 <- ifelse(is.na(row1),rowEq[1],row1)

              welk1 <- which(RrE[[h]][row1,]!=0)
              isna_h <- is.na(unique_h[welk1])
              if(sum(isna_h)==2){
                teller <- teller + 1
                unique_h[welk1] <- teller
              }else{ #one is already assigned a unique code
                unique_h[welk1] <- unique_h[welk1[!isna_h]]
              }
              rowEq <- rowEq[-which(rowEq==row1)]
            }
            if(sum(is.na(unique_h))>0){ #unconstrained rho's receive unique code
              unique_h[is.na(unique_h)] <- teller + 1:sum(is.na(unique_h))
              teller <- teller + sum(is.na(unique_h))
            }
          }
          if(!is.null(RrO[[h]])){
            unicum <- unique(unique_h[unique_h!=0])
            inequalities_h <- matrix(0,nrow=nrow(RrO[[h]]),ncol=max(unicum)+1)
            inequalities_h[,max(unicum)+1] <- RrO[[h]][,numrho+1]
            for(u in sort(unicum)){

              welk <- which(unique_h == u)
              if(length(welk) > 1){
                if(nrow(RrO[[h]]) > 1){
                  inequalities_h[,u] <- apply(RrO[[h]][,which(unique_h == u)],1,sum)
                }else{
                  inequalities_h[,u] <- apply(t(RrO[[h]][,which(unique_h == u)]),1,sum)
                }
              }else{ #length is 1
                inequalities_h[,u] <- RrO[[h]][,which(unique_h == u)]
              }
            }
          } else inequalities_h = NULL

          if(is.null(RrE[[h]])){
            #only order constraints; use output from unconstrained analysis
            #use unconstrained prior and posterior draws for computation of posterior probability
            if(is.matrix(RrO_h[,-(numrho+1)])){
              RO_h <- RrO_h[,-(numrho+1)]
            }else{
              RO_h <- t(RrO_h[,-(numrho+1)])
            }
            rO_h <- RrO_h[,numrho+1]
            postProbHt <- mean(apply((postdrawsHu%*% t(RO_h) - rep(1,nrow(postdrawsHu))%*%t(rO_h))>0,1,prod))
            priorProbHt <- mean(apply((priordrawsHu%*% t(RO_h) - rep(1,nrow(priordrawsHu))%*%t(rO_h))>0,1,prod))
            marglike_h <- marglikeHu$logmarglike + log(postProbHt) - log(priorProbHt)

            return(c(marglike=marglike_h,postprob=postProbHt,priorprob=priorProbHt,1))
          }else{
            if(sum(unique_h)==0){
              # all rho's are restricted to zero
              margLikeHt <- lnmarglik.p(y, X, W=Wlist[[1]], rho0=0)
              return(c(margLikeHt,1,1,ifelse(is.null(RrE[[h]]),1,0)))
            }else{
              Wlist_h <- lapply(1:max(unique_h),function(uni){
                Reduce("+",Wlist[which(unique_h==uni)])
              })
              numrho_h <- length(Wlist_h)
              if(numrho_h==1){
                # 1 free rho parameter under constrained model
                EW_h <- Re(eigen(Wlist_h[[1]],only.values=TRUE)$values)
                lb_h <- 1/min(EW_h)
                ub_h <- 1/max(EW_h)
                if(is.null(inequalities_h)){
                  marglikeHu <- lnmarglik.n(gridpoints=1e3, y, X, W=Wlist_h[[1]], mean=prior_mean[which(unique_h==1)[1]],
                                            sd=sqrt(prior_covm[which(unique_h==1)[1],which(unique_h==1)[1]]), a1 = lb_h,
                                            a2 = ub_h)
                  return(c(marglikeHu[[1]],1,1,0))
                }else{
                  drawscheck <- runif(1e4,min=lb_h,max=ub_h)
                  welkcheck <- drawscheck[apply(drawscheck %*% t(inequalities_h[,1]) - rep(1,1e4)%*%t(inequalities_h[,2]) > 0, 1,prod)==1]
                  rm(drawscheck)
                  lb_h2 <- round(min(welkcheck),2)
                  ub_h2 <- round(max(welkcheck),2)
                  marglikeHu <- lnmarglik.n(gridpoints=1e3, y, X, W=Wlist_h[[1]], mean=prior_mean[which(unique_h==1)[1]],
                                            sd=sqrt(prior_covm[which(unique_h==1)[1],which(unique_h==1)[1]]), a1 = lb_h2,
                                            a2 = ub_h2)
                  priorProbHt <- ptruncnorm(ub_h2, a=lb_h, b=ub_h, mean=prior_mean[which(unique_h==1)[1]],
                                            sd = sqrt(prior_covm[which(unique_h==1)[1],which(unique_h==1)[1]])) -
                    ptruncnorm(lb_h2, a=lb_h, b=ub_h, mean = prior_mean[which(unique_h==1)[1]],
                               sd = sqrt(prior_covm[which(unique_h==1)[1],which(unique_h==1)[1]]))
                  start_h <- c(mean(postdrawsAllHu$rho.draws),apply(postdrawsAllHu$beta.draws,2,mean),
                               mean(postdrawsAllHu$sigma2.draws))
                  postdraws_h <- MCMC.N(N=2e3, y, X, W=Wlist_h[[1]], m=prior_mean[which(unique_h==1)[1]],
                                        std2=sqrt(prior_covm[which(unique_h==1)[1],which(unique_h==1)[1]]),
                                        start_h, burn=2e2)
                  postProbHt <- mean(postdraws_h$rho.draws > lb_h2 & postdraws_h$rho.draws < ub_h2)
                  rm(postdraws_h)
                  return(c(marglikeHu[[1]],postProbHt,priorProbHt,0))
                }
              }else{
                # number of free rho's is at least 2
                if(numrho_h<3){
                  if(max(sapply(Wlist,rowSums))==1){
                    inspace_h <- in.space.2.rsd
                  }else{
                    inspace_h <- in.space.2
                  }
                }else{
                  inspace_h <- in.space.R
                }
                #marginal likelihood under unconstrained model
                marglike_h <- marglikeNAM(S=1e3,y,X,Wlist=Wlist_h,mu_prior=prior_mean[1:numrho_h],
                                          Sigma_prior=prior_covm[1:numrho_h,1:numrho_h],
                                          orderconstraints = inequalities_h, in.space=inspace_h)
                return(c(marglike_h[[1]],marglike_h$postprob,marglike_h$priorprob,0))
              }
            }
          }
        }
      })), nrow=numhyp, byrow=TRUE)

      if(complement == TRUE){
        #check if complete parameter space is
        if(sum(output_marglike[,4])==0){
          #complement is the same as unconstrained model
          output_marglike <- rbind(output_marglike,c(unlist(marglikeHu)[1:3],1))
        }else{#completement is the complement space of the order constrained models
          which_order <- which(output_marglike[,4]==1)
          if(length(which_order)==1){
            postProbHc <- 1 - sum(output_marglike[which_order,2])
            priorProbHc <- 1 - sum(output_marglike[which_order,3])
            marglikeHc <- marglikeHu[[1]] + log(marglikeHu[[2]]) - log(priorProbHc) + log(postProbHc)
            output_marglike <- rbind(output_marglike,c(marglikeHc,postProbHc,priorProbHc,1))
          }else{
            # check if order constrained models are overlapping
            # get draws that satisfy the constraints of the separate order constrained hypotheses
            numpriordraws <- nrow(priordrawsHu)
            numpostdraws <- nrow(postdrawsHu)
            checksOC <- lapply(which_order,function(h){
              Rorder <- as.matrix(RrO[[h]][,-(1+numrho)])
              if(ncol(Rorder)==1 & numrho>1){
                Rorder <- t(Rorder)
              }
              rorder <- as.matrix(RrO[[h]][,1+numrho])
              apply(priordrawsHu %*%t(Rorder) > rep(1,numpriordraws)%*%t(rorder),1,prod)
            })
            checkOCplus <- Reduce("+",checksOC)

            if(sum(checkOCplus > 0) < numpriordraws){ #then the joint order constrained hypotheses do not completely cover the parameter space.
              if(sum(checkOCplus > 1) == 0){ # then order constrained spaces are nonoverlapping
                postProbHc <- 1 - sum(output_marglike[which_order,2])
                priorProbHc <- 1 - sum(output_marglike[which_order,3])
                marglikeHc <- marglikeHu[[1]] + log(marglikeHu[[2]]) - log(priorProbHc) + log(postProbHc)
                output_marglike <- rbind(output_marglike,c(marglikeHc,postProbHc,priorProbHc,1))
              }else{# the order constrained subspaces at least partly overlap
                #compute prior probability
                checksOCprior <- lapply(which_order,function(h){
                  Rorder <- as.matrix(RrO[[h]][,-(1+numrho)])
                  if(ncol(Rorder)==1){
                    Rorder <- t(Rorder)
                  }
                  rorder <- as.matrix(RrO[[h]][,1+numrho])
                  apply(priordrawsHu%*%t(Rorder) > rep(1,numpriordraws)%*%t(rorder),1,prod)
                })
                priorProbHc <- sum(Reduce("+",checksOCprior) == 0) / numpriordraws

                #compute posterior probability
                checksOCpost <- lapply(which_order,function(h){
                  Rorder <- as.matrix(RrO[[h]][,-(1+numrho)])
                  if(ncol(Rorder)==1){
                    Rorder <- t(Rorder)
                  }
                  rorder <- as.matrix(RrO[[h]][,1+numrho])
                  apply(postdrawsHu%*%t(Rorder) > rep(1,numpostdraws)%*%t(rorder),1,prod)
                })
                postProbHc <- sum(Reduce("+",checksOCpost) == 0) / numpostdraws
                marglikeHc <- marglikeHu[[1]] + log(marglikeHu[[2]]) - log(priorProbHc) + log(postProbHc)
                output_marglike <- rbind(output_marglike,c(marglikeHc,postProbHc,priorProbHc,1))
              }
            }
          }
        }
        row.names(output_marglike) <- c(parse_hyp$original_hypothesis,"complement")
      }else{
        row.names(output_marglike) <- c(parse_hyp$original_hypothesis)
      }
      colnames(output_marglike) <- c("margLikeHt","postProbHt","priorProbHt","onlyOrder")
      #create output
      relcomp <- matrix(c(rep(NA,nrow(output_marglike)),output_marglike[,3]),ncol=2)
      relfit <- matrix(c(rep(NA,nrow(output_marglike)),output_marglike[,2]),ncol=2)
      #compute log marginal likelihood for H* without order constraints
      BF_E <- exp(output_marglike[,1] - log(output_marglike[,2]) + log(output_marglike[,3]) - marglikeHu[[1]])
      BFtu_confirmatory <- exp(output_marglike[,1] - marglikeHu[[1]])
      #compute BFmatrix and PHPs
      numhyp <- length(RrE)
      logBFmatrix <- matrix(rep(output_marglike[,1],numhyp+complement),nrow=numhyp+complement) -
        matrix(rep(output_marglike[,1],each=numhyp+complement),nrow=numhyp+complement)
      if(complement == TRUE){
        row.names(logBFmatrix) <- colnames(logBFmatrix) <- c(parse_hyp$original_hypothesis,"complement")
      }else{
        row.names(logBFmatrix) <- colnames(logBFmatrix) <- c(parse_hyp$original_hypothesis)
      }
      BFmatrix_confirmatory <- round(exp(logBFmatrix),3)
      BFta_confirmatory <- exp(output_marglike[,1] - max(output_marglike[,1]))
      # Change prior probs in case of default setting
      if(is.null(prior.hyp)){
        priorprobs <- rep(1/length(BFta_confirmatory),length(BFta_confirmatory))
      }else{
        if(!is.numeric(prior.hyp) || length(prior.hyp)!=length(BFta_confirmatory)){
          warning(paste0("Argument 'prior.hyp' should be numeric and of length ",as.character(length(BFta_confirmatory)),
                         ". Equal prior probabilities are used."))
          priorprobs <- rep(1/length(BFta_confirmatory),length(BFta_confirmatory))
        }else{
          priorprobs <- prior.hyp
        }
      }
      PHP_confirmatory <- priorprobs*BFta_confirmatory / sum(priorprobs*BFta_confirmatory)
      BFtable <- cbind(relcomp,relfit,BF_E,relfit[,2]/relcomp[,2],
                       BF_E*relfit[,2]/relcomp[,2],PHP_confirmatory)
      row.names(BFtable) <- names(PHP_confirmatory)
      colnames(BFtable) <- c("comp_E","comp_O","fit_E","fit_O","BF_E","BF_O","BF","PHP")
      hypotheses <- names(PHP_confirmatory)

    }else{
      #confirmatory test on beta's using AAFBF
      Args$hypothesis <- hypothesis
      Args$prior.hyp <- prior.hyp
      Args$complement <- complement
      out <- do.call(BF, Args)
      BFmatrix_confirmatory <- out$BFmatrix_confirmatory
      BFtu_exploratory <- out$BFtu_exploratory
      BFtu_confirmatory <- out$BFtu_confirmatory
      PHP_confirmatory <- out$PHP_confirmatory
      BFtable <- out$BFtable
      priorprobs <- out$priorprobs
      hypotheses <- out$hypotheses
    }
  }else{
    BFmatrix_confirmatory <- PHP_confirmatory <- BFtu_confirmatory <-
      hypotheses <- BFtable <- priorprobs <- NULL
  }

  BF_out <- list(
    BFtu_exploratory=BFtu_exploratory,
    PHP_exploratory=PHP_exploratory,
    BFtu_confirmatory=BFtu_confirmatory,
    PHP_confirmatory=PHP_confirmatory,
    BFmatrix_confirmatory=BFmatrix_confirmatory,
    BFtable_confirmatory=BFtable,
    prior=priorprobs,
    hypotheses=hypotheses,
    estimates=postestimates,
    model=x,
    bayesfactor="Bayes factors based on normal priors & AAFBF",
    parameter="network autocorrelations & regression coefficients",
    call=match.call())

  class(BF_out) <- "BF"

  return(BF_out)

}




