# BANAM

<img src="man/figures/logo_BANAM.png" width = 450 />

R package for Bayesian Analysis of the Network Autocorrelation Model. The methodology
is based on:

Dittrich, D., Leenders, R.Th.A.J. Leenders, & Mulder, J. (2017). Bayesian estimation of the network autocorrelation model. Social Network, 48, 213–236. https://doi.org/10.1016/j.socnet.2016.09.002

Dittrich, D., Leenders, R.Th.A.J. Leenders, & Mulder, J. (2019). Network autocorrelation modeling: A Bayes factor approach for testing (multiple) precise and interval hypotheses. Sociological Methods & Research. https://doi.org/10.1177/0049124117729712

Dittrich, D., Leenders, R.Th.A.J. Leenders, & Mulder, J. (2020). Network Autocorrelation Modeling: Bayesian Techniques for Estimating and Testing Multiple Network Autocorrelations. Sociological Methodology. https://doi.org/10.1177/0081175020913899

Last Modified 09/24/2020

Licensed under the GNU General Public License version 2 (June, 1991)

Basic example
-------------

``` r
library(BANAM)
library(sna)
library(BFpack)

# Generate data

set.seed(3)
n <- 50
d1 <- .1
Wadj1 <- sna::rgraph(n, tprob=d1, mode="graph")
W1 <- sna::make.stochastic(Wadj1, mode="row")
d2 <- .3
Wadj2 <- sna::rgraph(n, tprob=d2, mode="graph")
W2 <- sna::make.stochastic(Wadj2, mode="row")
# set rho, beta, sigma2, and generate y
rho1 <- 0
rho2 <- .4
K <- 3
beta <- rnorm(K)
sigma2 <- 1
X <- matrix(c(rep(1, n), rnorm(n*(K-1))), nrow=n, ncol=K)
y <- c(solve(diag(n) - rho1*W1 - rho2*W2)%*%(X%*%beta + rnorm(n)))

# Bayesian estimation of a NAM with a single weight matrix
best1 <- banam(y,X,W1)

# Bayesian estimation of a NAM with two weight matrices
best2 <- banam(y,X,W=list(W1,W2))

# Bayesian hypothesis testing of equality/order constraints on network
# autocorrelation parameters
BFbest2 <- BFbanam(best2,hypothesis="rho1>rho2>0; rho1=rho2>0; rho1=rho2=0")

```

Installation
------------

You can install BANAM from github with:

``` r
# install.packages("devtools")
devtools::install_github("jomulder/BANAM")
