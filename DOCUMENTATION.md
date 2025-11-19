# LassoVAR Python Package - Complete Documentation

## Overview

This is a complete Python port of the R package `lassovar` by Laurent Callot. The package provides tools for estimating and forecasting Vector Autoregression (VAR) models using Lasso and adaptive Lasso estimators.

## Installation

```bash
# Install dependencies
pip install numpy pandas scikit-learn scipy statsmodels joblib

# Install the package
pip install -e .
```

## Package Structure

```
lassovar-python/
├── lassovar/                 # Main package directory
│   ├── __init__.py          # Package initialization
│   ├── lassovar.py          # Main LassoVAR class
│   ├── estimation.py        # Core estimation functions
│   ├── adaptive.py          # Adaptive weights computation
│   ├── forecasting.py       # Forecasting functionality
│   └── utils.py             # Utility functions
├── tests/                    # Test suite
│   ├── test_sim.py          # Basic estimation tests
│   ├── test_methods.py      # Method tests
│   └── test_forecast.py     # Forecasting tests
├── examples/                 # Usage examples
│   └── basic_usage.py       # Comprehensive examples
├── setup.py                  # Package setup
├── requirements.txt          # Dependencies
├── README.md                 # User documentation
├── LICENSE                   # MIT License
└── DOCUMENTATION.md          # This file
```

## Core Components

### 1. LassoVAR Class (`lassovar.py`)

The main class for fitting VAR models with Lasso or adaptive Lasso.

**Constructor Parameters:**
- `dat`: DataFrame with time series data (T x N)
- `exo`: Optional exogenous variables (T x K)
- `lags`: Number of lags (default: 1)
- `ic`: Information criterion ('BIC' or 'AIC', default: 'BIC')
- `adaptive`: Type of adaptive Lasso ('none', 'ols', 'lasso', 'ridge', 'group')
- `post`: Estimate post-Lasso OLS (default: False)
- `mc`: Enable parallel processing (default: False)
- `n_jobs`: Number of parallel jobs (default: -1)
- `dfmax`: Maximum degrees of freedom (default: None)
- `horizon`: h-step ahead horizon (default: 1)
- `trend`: Include linear trend (default: False)
- `lambda_values`: User-defined lambda values (default: None)

**Methods:**
- `fit()`: Estimate the model
- `predict(newdata)`: Make predictions
- `residuals()`: Extract residuals
- `summary(short=False)`: Print summary statistics
- `coef()`: Get coefficient matrix

**Example:**
```python
from lassovar import LassoVAR
import pandas as pd
import numpy as np

# Create data
data = pd.DataFrame(np.random.randn(100, 3), columns=['A', 'B', 'C'])

# Fit model
model = LassoVAR(data, lags=1, ic='BIC')
model.fit()

# Get results
model.summary()
predictions = model.predict(model.x.iloc[-1:])
residuals = model.residuals()
```

### 2. Forecasting (`forecasting.py`)

**Function: `forecast_lassovar`**

Performs pseudo out-of-sample forecasting experiments.

**Parameters:**
- `dat`: Time series data
- `exo`: Optional exogenous variables
- `fc_train`: Number of training observations (required)
- `horizon`: Forecast horizon (default: 1)
- `lags`: Number of lags (default: 1)
- `fc_window`: 'fix' or 'expanding' (default: 'fix')
- `fc_type`: 'recursive' or 'direct' (default: 'recursive')
- `ic`: Information criterion (default: 'BIC')
- `adaptive`: Adaptive type (default: 'none')
- `mc`: Parallel processing (default: False)
- `n_jobs`: Number of jobs (default: -1)
- `silent`: Suppress output (default: False)
- `trend`: Include trend (default: False)
- `post`: Post-Lasso OLS (default: False)

**Returns:**
Dictionary with:
- `err`: Forecast errors DataFrame
- `pred`: Predictions DataFrame
- `coefficients`: List of coefficient matrices
- `lambda`: Lambda values array
- `spectest`: Specification tests
- `call`: Call parameters

**Example:**
```python
from lassovar import forecast_lassovar

fc_results = forecast_lassovar(
    data,
    fc_train=80,
    horizon=1,
    lags=1,
    fc_window='expanding',
    silent=True
)

# Access results
print(fc_results['err'].mean())
```

### 3. Utility Functions (`utils.py`)

**`make_var_data(data_var, lags, horizon=1, exo=None, trend=False)`**
- Prepares data for VAR estimation by creating lagged variables

**`specification_tests(residuals, y_original, lags=None, fitdf=1)`**
- Performs Ljung-Box, Shapiro-Wilk, and R² tests

**`compute_bic(rss, n_obs, df)` / `compute_aic(rss, n_obs, df)`**
- Compute information criteria

**`ridge_degrees_of_freedom(x, lambda_values)`**
- Calculate effective degrees of freedom for ridge regression

### 4. Adaptive Weights (`adaptive.py`)

Functions for computing adaptive weights using different initial estimators:

**`compute_ols_weights(y, x, mc=False, n_jobs=-1, gamma=1)`**
- OLS-based adaptive weights

**`compute_lasso_weights(y, x, ic='BIC', mc=False, n_jobs=-1, ...)`**
- Lasso-based adaptive weights

**`compute_ridge_weights(y, x, ic='BIC', mc=False, n_jobs=-1, ...)`**
- Ridge-based adaptive weights

**`compute_group_weights(y, x, trend=False, gamma=1)`**
- Group Lasso adaptive weights (simplified implementation)

**`get_adaptive_weights(y, x, adaptive_type, ...)`**
- Main function to get adaptive weights based on type

### 5. Estimation (`estimation.py`)

**`lassovar_equation(y, x, ada_w, ic='BIC', ...)`**
- Core estimation function for VAR equations
- Implements equation-by-equation Lasso/Ridge estimation
- Handles adaptive weights and information criterion selection

**`post_ols_estimation(y, x, selected_params, mc=False, n_jobs=-1)`**
- Post-Lasso OLS on selected variables

**`fit_glmnet_equation(i, y, x, ada_w, ic, alpha, dfmax, trend, lambda_values)`**
- Fits a single equation using elastic net

## Implementation Details

### Differences from R Implementation

1. **Dependencies:**
   - R: glmnet, grpreg, Matrix, biglm, parallel
   - Python: scikit-learn, numpy, pandas, scipy, statsmodels, joblib

2. **Parallel Processing:**
   - R: Uses `mclapply` from parallel package
   - Python: Uses `joblib.Parallel`

3. **Sparse Matrices:**
   - R: Uses Matrix package
   - Python: Uses numpy arrays (sparse matrix support can be added)

4. **Lasso Implementation:**
   - R: Uses glmnet package
   - Python: Uses sklearn.linear_model.Lasso with custom lambda path

5. **Adaptive Lasso:**
   - R: Direct per-feature penalties in glmnet
   - Python: Feature scaling approach (equivalent mathematically)

### Key Algorithmic Components

#### 1. VAR Data Preparation
```python
# Creates lagged variables
for l in range(1, lags + 1):
    lagged = data.shift(l + horizon - 1)
    x_list.append(lagged)
```

#### 2. Information Criterion Selection
```python
# BIC
ic_values = np.log(rss / n_obs) + df * np.log(n_obs) / n_obs

# AIC
ic_values = np.log(rss / n_obs) + df / n_obs

# Select best
best_idx = np.argmin(ic_values)
```

#### 3. Adaptive Weights
```python
# Compute initial estimator
initial_coef = fit_initial_estimator(y, x)

# Compute adaptive weights
weights = np.abs(initial_coef) ** (-gamma)
```

#### 4. Specification Tests
```python
# Ljung-Box for autocorrelation
lb_result = acorr_ljungbox(residuals, lags=[lags])

# Shapiro-Wilk for normality
_, p_value = stats.shapiro(residuals)

# R²
r2 = 1 - (n * var(residuals)) / (var(y) * len(y))
```

## Testing

The package includes comprehensive tests matching the original R package:

### Running Tests

```bash
# Install test dependencies
pip install pytest pytest-cov

# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_sim.py -v

# Run with coverage
pytest tests/ --cov=lassovar --cov-report=html
```

### Test Files

1. **`test_sim.py`**: Basic estimation tests
   - Basic Lasso VAR
   - VAR with trend
   - VAR with AIC criterion
   - Multiple lags
   - Exogenous variables
   - Adaptive Lasso

2. **`test_methods.py`**: Method tests
   - Summary method
   - Predict method
   - Residuals method
   - Coefficient extraction
   - Error handling

3. **`test_forecast.py`**: Forecasting tests
   - Basic forecasting
   - Exogenous variables in forecasting
   - Expanding window
   - Multi-step ahead
   - Direct vs. recursive

## Usage Examples

### Example 1: Basic Lasso VAR

```python
import numpy as np
import pandas as pd
from lassovar import LassoVAR

# Generate data
np.random.seed(42)
data = pd.DataFrame(np.random.randn(100, 3), columns=['X1', 'X2', 'X3'])

# Fit model
model = LassoVAR(data, lags=1)
model.fit()

# Summary
model.summary()
```

### Example 2: Adaptive Lasso

```python
# Fit with OLS initial estimator
model = LassoVAR(data, lags=2, adaptive='ols', ic='BIC')
model.fit()

# Get coefficients
coef = model.coef()
print(f"Coefficient shape: {coef.shape}")
```

### Example 3: With Exogenous Variables

```python
# Create exogenous variables
nobs = 100
exo = pd.DataFrame({
    'trend': np.arange(nobs),
    'seasonal': np.sin(2 * np.pi * np.arange(nobs) / 12)
})

# Fit model
model = LassoVAR(data, lags=1, exo=exo, trend=True)
model.fit()
```

### Example 4: Forecasting

```python
from lassovar import forecast_lassovar

# Pseudo out-of-sample forecast
fc_results = forecast_lassovar(
    data,
    fc_train=80,
    horizon=1,
    lags=1,
    fc_window='expanding',
    fc_type='recursive',
    silent=True
)

# Analyze errors
mae = fc_results['err'].abs().mean()
rmse = np.sqrt((fc_results['err']**2).mean())
print(f"MAE: {mae}")
print(f"RMSE: {rmse}")
```

### Example 5: Parallel Processing

```python
# Enable parallel processing
model = LassoVAR(data, lags=2, mc=True, n_jobs=-1)
model.fit()

# Parallel forecasting
fc_results = forecast_lassovar(
    data,
    fc_train=80,
    horizon=1,
    lags=1,
    mc=True,
    n_jobs=4,
    silent=True
)
```

## Performance Considerations

1. **Memory Usage:**
   - Large VARs can be memory-intensive
   - Use `dfmax` parameter to limit model size

2. **Computational Speed:**
   - Enable parallel processing with `mc=True`
   - Adjust `n_jobs` based on available cores
   - For forecasting, parallelize at forecast level (not equation level)

3. **Numerical Stability:**
   - Feature standardization is automatic in Lasso
   - Ridge regression (`alpha=0`) for ill-conditioned problems
   - Adaptive Lasso helps with variable selection

## Limitations and Future Work

### Current Limitations

1. **Group Lasso:**
   - Simplified implementation
   - Full group structure not yet implemented
   - Would require additional dependencies

2. **Sparse Matrices:**
   - Not fully optimized for sparse storage
   - Could use scipy.sparse for large models

3. **Direct Forecasting:**
   - Implemented but may need further testing
   - Recursive is more thoroughly tested

### Potential Enhancements

1. Add support for:
   - Time-varying parameters
   - Structural breaks
   - Heteroskedasticity
   - Bootstrap confidence intervals

2. Performance improvements:
   - Cython extensions for critical loops
   - Better sparse matrix support
   - GPU acceleration for large models

3. Additional features:
   - Impulse response functions
   - Forecast error variance decomposition
   - Granger causality tests
   - Model diagnostics plots

## References

1. Callot, L. (2015). lassovar: Estimation and forecasting with VAR models using the (adaptive) Lasso. R package.

2. Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society: Series B, 58(1), 267-288.

3. Zou, H. (2006). The adaptive lasso and its oracle properties. Journal of the American Statistical Association, 101(476), 1418-1429.

4. Callot, L. A., & Kock, A. B. (2014). Oracle inequalities for high dimensional vector autoregressions. Journal of Econometrics, 186(2), 325-344.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## License

MIT License - see LICENSE file

## Support

For issues, questions, or contributions:
- GitHub Issues: [github.com/lcallot/lassovar](https://github.com/lcallot/lassovar)
- Original R package documentation

## Acknowledgments

- Original R package by Laurent Callot
- Contributors to scikit-learn, numpy, pandas, and other dependencies
- Research papers on Lasso VAR methods
