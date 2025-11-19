lassovar (Python)
==================

Estimation and forecasting of VAR model with the Lasso.

This is a Python port of the original R package [lassovar](https://github.com/lcallot/lassovar) by Laurent Callot.

The package provides tools for estimating Vector Autoregressions with Lasso-type estimators and is used in:
- [*Oracle inequalities for high dimensional vector autoregressions*](http://lcallot.github.io/pub/oracle-var)
- [*Oracle Efficient estimation and Forecasting with the Adaptive Lasso and the Adaptive Group Lasso in Vector Autoregressions.*](http://lcallot.github.io/pub/oracle-forecasting)
- [*Estimation and Forecasting of Large Realized Covariance Matrices and Portfolio Choice.*](http://lcallot.github.io/wp/rcv-fc/)

## Features

* Estimation of Vector Autoregressions with the Lasso or adaptive Lasso using either the Lasso, OLS, or ridge regressions as the initial estimator
* Penalty parameter selection using information criteria (BIC or AIC)
* Post-Lasso OLS estimation
* Forecasting support (direct and recursive)
* Parallel computation support
* Summary, residuals, and prediction methods

## Installation

```bash
pip install -e .
```

## Usage

### Basic VAR estimation with Lasso

```python
import numpy as np
import pandas as pd
from lassovar import LassoVAR

# Create sample data
data = pd.DataFrame(np.random.randn(100, 5), columns=['V1', 'V2', 'V3', 'V4', 'V5'])

# Fit a VAR(1) model with Lasso
model = LassoVAR(data, lags=1)
model.fit()

# Summary statistics
print(model.summary())

# Get predictions
predictions = model.predict(model.x[-1:])

# Get residuals
residuals = model.residuals()
```

### Adaptive Lasso

```python
# Fit with adaptive Lasso using OLS as initial estimator
model = LassoVAR(data, lags=2, adaptive='ols')
model.fit()
```

### Forecasting

```python
from lassovar import forecast_lassovar

# Perform pseudo out-of-sample forecasting experiment
fc_results = forecast_lassovar(
    data, 
    fc_train=80,  # Use first 80 observations for training
    horizon=1,     # 1-step ahead forecast
    lags=1,
    fc_window='expanding',  # Expanding window
    fc_type='recursive'     # Recursive forecast
)

# Access forecast errors
print(fc_results['err'])
```

### With Exogenous Variables

```python
# Create exogenous variables
exo = pd.DataFrame(np.random.randn(100, 2), columns=['Exo1', 'Exo2'])

# Fit model with exogenous variables
model = LassoVAR(data, lags=1, exo=exo)
model.fit()
```

## API Reference

### LassoVAR Class

**Parameters:**
- `dat`: DataFrame containing the time series data
- `exo`: Optional DataFrame with exogenous variables (not lagged)
- `lags`: Number of lags to include (default: 1)
- `ic`: Information criterion for penalty selection ('BIC' or 'AIC', default: 'BIC')
- `adaptive`: Initial estimator for adaptive Lasso ('none', 'ols', 'lasso', 'ridge', default: 'none')
- `post`: Whether to compute post-Lasso OLS (default: False)
- `mc`: Enable parallel processing (default: False)
- `n_jobs`: Number of parallel jobs (default: -1, uses all cores)
- `dfmax`: Maximum number of variables excluding intercept (default: None)
- `horizon`: Horizon for h-step ahead estimation (default: 1)
- `trend`: Include linear trend (default: False)
- `lambda_values`: User-defined lambda values (default: None)

**Methods:**
- `fit()`: Estimate the VAR model
- `predict(newdata)`: Make predictions with new data
- `residuals()`: Extract residuals
- `summary()`: Print summary statistics
- `coef()`: Get coefficients matrix

### forecast_lassovar Function

Perform pseudo out-of-sample forecasting experiments.

**Parameters:**
- `dat`: DataFrame containing the time series
- `exo`: Optional exogenous variables
- `fc_train`: Number of training observations
- `horizon`: Forecast horizon (default: 1)
- `lags`: Number of lags (default: 1)
- `fc_window`: 'fix' or 'expanding' (default: 'fix')
- `fc_type`: 'recursive' or 'direct' (default: 'recursive')
- `ic`: Information criterion ('BIC' or 'AIC', default: 'BIC')
- `adaptive`: Adaptive Lasso type (default: 'none')
- `mc`: Enable parallel processing (default: False)
- `n_jobs`: Number of parallel jobs (default: -1)
- `silent`: Suppress output (default: False)
- `trend`: Include trend (default: False)
- `post`: Post-Lasso OLS (default: False)

**Returns:** Dictionary with forecast results including:
- `err`: Forecast errors
- `pred`: Predictions
- `coefficients`: List of coefficient matrices for each forecast
- `lambda`: Lambda values used
- `spectest`: Specification tests

## Requirements

- Python >= 3.7
- numpy >= 1.19.0
- pandas >= 1.1.0
- scikit-learn >= 0.23.0
- scipy >= 1.5.0
- statsmodels >= 0.12.0
- joblib >= 0.16.0

## License

MIT License

Copyright (c) 2015 Laurent Callot

## Citation

If you use this package in your research, please cite the original R package and relevant papers.

## Disclaimer

This package is a work in progress. The Python port aims to maintain feature parity with the original R implementation.
