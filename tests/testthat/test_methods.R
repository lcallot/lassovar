library(lassovar)

# A simple VAR(1) parameter matrix
A <- matrix(0,3,3)
diag(A) <- 0.4
A[2,1] <- A[1,3] <- -0.4
# Simulating data
nobs <- 500
simdata <- matrix(rnorm(3),nobs+1,3)
simdata[1,] <- rnorm(3)
for(t in 1:nobs) simdata[t+1,] <- simdata[t,]%*%A + rnorm(3)
simdata <- tail(simdata,-1)

colnames(simdata) <- c('Cyrus','Cambyses','Darius')

# TEST METHODS
# Estimating a model
mod <- lassovar(simdata,lags = 1,exo=data.frame('exovar'=1:500))
# testing summary
expect_that(summary(mod),is_a('matrix'))
# some fake data for forecasting
nwd <- cbind(matrix(rnorm(6),2,3),nobs+1:2)
# testing predict	
expect_that(predict(mod,nwd),is_a('matrix'))
# get residuals
expect_that(residuals(mod),is_a('data.frame'))
expect_that(resid(mod),is_a('data.frame'))

