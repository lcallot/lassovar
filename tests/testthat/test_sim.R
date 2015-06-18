library(lassovar)

# A simple VAR(1) parameter matrix
A <- matrix(0,3,3)
diag(A) <- 0.4
A[2,1] <- A[1,3] <- -0.4
# Simulating data
nobs <- 50
simdata <- matrix(rnorm(3),nobs+1,3)
simdata[1,] <- rnorm(3)
for(t in 1:nobs) simdata[t+1,] <- simdata[t,]%*%A + rnorm(3)
simdata <- tail(simdata,-1)

colnames(simdata) <- c('Cyrus','Cambyses','Darius')

# TEST lassovar

# basic
expect_that(lassovar(simdata,lags = 1),is_a('lassovar'))
# AIC
expect_that(lassovar(simdata,lags = 1,ic='AIC'),is_a('lassovar'))
# multiple lags
expect_that(lassovar(simdata,lags = 2),is_a('lassovar'))
# exo variables
expect_that(lassovar(simdata,lags = 1,exo=c(1:nobs)),is_a('lassovar'))
# adaptive lasso
expect_that(lassovar(simdata,lags = 1,adaptive='ols'),is_a('lassovar'))


# TEST METHODS

# Estimating a model
mod <- lassovar(simdata,lags = 1)
# testing summary
expect_that(summary(mod),is_a('matrix'))
# some fake data for forecasting
nwd <- matrix(rnorm(6),2,3)
# testing predict	
expect_that(predict(mod,nwd),is_a('matrix'))
# get residuals
expect_that(residuals(mod),is_a('data.frame'))
expect_that(resid(mod),is_a('data.frame'))



# TEST FORECASTS

# basic
expect_that(forecast.lassovar(dat = simdata,ntrain = 48,horizon = 1,fc.window = 'fix',fc.type = 'recursive'),is_a('list'))
# expanding
expect_that(forecast.lassovar(dat = simdata,ntrain = 48,horizon = 1,fc.window = 'expanding',fc.type = 'recursive'),is_a('list'))
# horizon > 1 
expect_that(forecast.lassovar(dat = simdata,ntrain = 47,horizon = 2,fc.window = 'fix',fc.type = 'recursive'),is_a('list'))
# direct forecasts
expect_that(forecast.lassovar(dat = simdata,ntrain = 47,horizon = 2,fc.window = 'fix',fc.type = 'direct'),is_a('list'))

