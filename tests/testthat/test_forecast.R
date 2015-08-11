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


# TEST FORECASTS

# basic
expect_that(forecast.lassovar(dat = simdata,fc.train = 48,horizon = 1,fc.window = 'fix',fc.type = 'recursive'),is_a('list'))
# exo
expect_that(forecast.lassovar(dat = simdata,exo = matrix(rnorm(nobs),ncol=1),fc.train = 48,horizon = 1,fc.window = 'fix',fc.type = 'recursive'),is_a('list'))
# expanding
expect_that(forecast.lassovar(dat = simdata,fc.train = 48,horizon = 1,fc.window = 'expanding',fc.type = 'recursive',silent = TRUE),is_a('list'))
# horizon > 1 
expect_that(forecast.lassovar(dat = simdata,fc.train = 47,horizon = 2,fc.window = 'fix',fc.type = 'recursive',silent = TRUE),is_a('list'))
# direct forecasts
expect_that(forecast.lassovar(dat = simdata,fc.train = 47,horizon = 2,fc.window = 'fix',fc.type = 'direct',silent = TRUE),is_a('list'))

