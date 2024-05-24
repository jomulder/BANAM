
#' @method BF banam
#' @export
BF.banam <- function(x,
                    hypothesis = NULL,
                    #prior.hyp.explo = NULL,
                    #prior.hyp.conf = NULL,
                    prior.hyp = NULL,
                    complement = TRUE,
                    #log = FALSE,
                    ...){

  numrho <- length(x$W)
  n <- length(x$y)
  if(is.null(x$prior.mean)){x$prior.mean <- rep(0,numrho)}
  if(is.null(x$prior.Sigma)){x$prior.Sigma <- diag(numrho)*1e6}

  if(numrho == 1){
    if(x$prior == "flat" | (x$prior == "normal" & x$prior.mean==0 & x$prior.Sigma>=1e6)){
      #use draws from input object
      post.draws <- cbind(x$rho.draws,x$beta.draws)
      colnames(post.draws) <- names(get_estimates(x)$estimate)
    }else{
      banam1 <- banam(y = x$y, X = x$X, W = x$W, prior = "flat",...)
      post.draws <- cbind(banam1$rho.draws,banam1$beta.draws)
      colnames(post.draws) <- names(get_estimates(banam1)$estimate)
    }
  }else{
    if(sum(abs(x$prior.mean))==0 & sum(abs(x$prior.Sigma - diag(numrho)*1e6))==0){
      #use draws from input object
      post.draws <- cbind(x$rho.draws,x$beta.draws)
      colnames(post.draws) <- names(get_estimates(x)$estimate)
    }else{
      banam1 <- banam(y = x$y, X = x$X, W = x$W,...)
      post.draws <- cbind(banam1$rho.draws,banam1$beta.draws)
      colnames(post.draws) <- names(get_estimates(banam1)$estimate)
    }
  }
  post.mean <- apply(post.draws,2,mean)
  post.covm <- cov(post.draws)
  BF.BANAM <- BF(x = post.mean,
     Sigma = post.covm,
     n = n,
     #log = log,
     hypothesis = hypothesis,
     #prior.hyp.explo = prior.hyp.explo,
     #prior.hyp.conf = prior.hyp.conf,
     prior.hyp = prior.hyp
     )

  BF.BANAM$parameter <- "NAM parameters"

  return(BF.BANAM)
}


