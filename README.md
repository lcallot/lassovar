lassovar
========

Estimation and forecasting of VAR model with the Lasso.

[![Build Status](https://travis-ci.org/lcallot/lassovar.svg?branch=master)](https://travis-ci.org/lcallot/lassovar)

This package is a wrapper for the _glmnet_ package aimed at facilitating estimation and forecasting with VAR models.
The package is used in:
- [*Oracle inequalities for high dimensional vector autoregressions*](http://lcallot.github.io/pub/oracle-var)
- [*Oracle Efficient estimation and Forecasting with the Adaptive Lasso and the Adaptive Group Lasso in Vector Autoregressions.*](http://lcallot.github.io/pub/oracle-forecasting)
- [*Estimation and Forecasting of Large Realized Covariance Matrices and Portfolio Choice.*](http://lcallot.github.io/wp/rcv-fc/)



Disclaimer
---------
This package is a work in progress.


Usage
--------

* The function `lassovar` provides for the estimation of Vector Autoregressions with the Lasso,
or adaptive Lasso using either the Lasso, OLS, or ridge regressions as the initial estimator.
The penalty parameter is chosen using an information criterion (BIC or AIC). A post-Lasso OLS can also be estimated.  

* `forecast.lassovar` is used to forecast (!) either directly or recursively. 

* `summary`, `residuals`, and `predict` methods are provided. 
