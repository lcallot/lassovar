library(lassovar)

# A simple VAR(1) parameter matrix
A <- matrix(0,3,3)
diag(A) <- 0.4
A[2,1] <- A[1,3] <- -0.4
# Simulating data (VAR(1))
nobs <- 20
simdata <- matrix(rnorm(3),nobs+1,3)
simdata[1,] <- rnorm(3)
for(t in 1:nobs) simdata[t+1,] <- simdata[t,]%*%A + rnorm(3)
simdata <- tail(simdata,-1)

colnames(simdata) <- c('Cyrus','Cambyses','Darius')

# TEST lassovar

# basic
expect_that(lassovar(simdata,lags = 1),is_a('lassovar'))
# trend
expect_that(lassovar(simdata,lags = 1,trend=TRUE),is_a('lassovar'))
# AIC
expect_that(lassovar(simdata,lags = 1,ic='AIC'),is_a('lassovar'))
# multiple lags
expect_that(lassovar(simdata,lags = 2),is_a('lassovar'))
# exo variables
expect_that(lassovar(simdata,lags = 1,exo=c(1:nobs)),is_a('lassovar'))
# adaptive lasso
expect_that(lassovar(simdata,lags = 1,adaptive='ols'),is_a('lassovar'))


