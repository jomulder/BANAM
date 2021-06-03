#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom extraDistr rinvgamma dinvgamma
#' @importFrom Matrix nearPD
#' @importFrom stats rgamma cov2cor cov
#' @importFrom truncnorm ptruncnorm rtruncnorm
#' @importFrom matrixcalc is.positive.definite
#' @importFrom rARPACK eigs
#' @importFrom stats rnorm pnorm dnorm lm sigma sd runif quantile
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom psych tr
#' @importFrom sna lnam
#' @importFrom BFpack BF
#' @importFrom bain get_estimates
parse_hypothesis <- getFromNamespace("parse_hypothesis", "bain")
make_RrList2 <- getFromNamespace("make_RrList2", "BFpack")





