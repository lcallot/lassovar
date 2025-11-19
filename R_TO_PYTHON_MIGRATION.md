# Migration Guide: R lassovar to Python lassovar

This guide helps users migrate from the R version of lassovar to the Python version.

## Installation

**R:**
```r
# install.packages("devtools")
devtools::install_github("lcallot/lassovar")
library(lassovar)
```

**Python:**
```bash
pip install numpy pandas scikit-learn scipy statsmodels joblib
pip install -e .  # Install from source
```

```python
from lassovar import LassoVAR, forecast_lassovar
```

## Basic Usage Comparison

### 1. Basic Lasso VAR Estimation

**R:**
```r
library(lassovar)

# Create data
dat <- data.frame(matrix(rnorm(300), ncol=3))
colnames(dat) <- c('V1', 'V2', 'V3')

# Fit model
lv.mod <- lassovar(dat, lags=1)

# Summary
summary(lv.mod)

# Predictions
predict(lv.mod, newdata)

# Residuals
residuals(lv.mod)
```

**Python:**
```python
import numpy as np
import pandas as pd
from lassovar import LassoVAR

# Create data
dat = pd.DataFrame(
    np.random.randn(100, 3),
    columns=['V1', 'V2', 'V3']
)

# Fit model
model = LassoVAR(dat, lags=1)
model.fit()

# Summary
model.summary()

# Predictions
model.predict(newdata)

# Residuals
model.residuals()
```

### 2. Adaptive Lasso

**R:**
```r
# With OLS initial estimator
lv.mod.ada <- lassovar(dat, lags=2, adaptive='ols')
```

**Python:**
```python
# With OLS initial estimator
model_ada = LassoVAR(dat, lags=2, adaptive='ols')
model_ada.fit()
```

### 3. With Exogenous Variables

**R:**
```r
# Create exogenous variables
exo <- data.frame(exo1=rnorm(100), exo2=1:100)

# Fit model
lv.mod <- lassovar(dat, lags=1, exo=exo)
```

**Python:**
```python
# Create exogenous variables
exo = pd.DataFrame({
    'exo1': np.random.randn(100),
    'exo2': np.arange(100)
})

# Fit model
model = LassoVAR(dat, lags=1, exo=exo)
model.fit()
```

### 4. Forecasting

**R:**
```r
# Pseudo out-of-sample forecasting
fc.lv <- forecast.lassovar(
    dat, 
    fc.train=80,
    horizon=1,
    lags=1,
    fc.window='expanding',
    fc.type='recursive',
    silent=TRUE
)

# Access forecast errors
fc.lv$err
```

**Python:**
```python
# Pseudo out-of-sample forecasting
fc_results = forecast_lassovar(
    dat,
    fc_train=80,
    horizon=1,
    lags=1,
    fc_window='expanding',
    fc_type='recursive',
    silent=True
)

# Access forecast errors
fc_results['err']
```

## Parameter Mapping

| R Parameter | Python Parameter | Notes |
|-------------|------------------|-------|
| `dat` | `dat` | Both accept DataFrames |
| `exo` | `exo` | Optional exogenous variables |
| `lags` | `lags` | Number of lags |
| `ic` | `ic` | 'BIC' or 'AIC' |
| `adaptive` | `adaptive` | 'none', 'ols', 'lasso', 'ridge', 'group' |
| `mc` | `mc` | Parallel processing flag |
| `ncores` | `n_jobs` | Number of cores (-1 for all) |
| `dfmax` | `dfmax` | Max degrees of freedom |
| `post` | `post` | Post-Lasso OLS |
| `horizon` | `horizon` | Forecast horizon |
| `trend` | `trend` | Include linear trend |
| `lambda` | `lambda_values` | User-defined lambda values |

## Return Value Comparison

### LassoVAR Object

**R (list components):**
- `$call` - The call
- `$var.names` - Variable names
- `$ada.w` - Adaptive weights
- `$x` - Right-hand side
- `$y` - Left-hand side
- `$coefficients` - Coefficients
- `$RSS` - Residual sum of squares
- `$lambda` - Lambda values
- `$spectest` - Specification tests
- `$estimator` - Estimator name
- `$ic` - Information criterion
- `$nbreq` - Number of equations

**Python (object attributes):**
- `.call` - Call dictionary
- `.var_names` - Variable names list
- `.ada_w` - Adaptive weights dict
- `.x` - Right-hand side DataFrame
- `.y` - Left-hand side DataFrame
- `.coefficients` - Coefficients array
- `.RSS` - Residual sum of squares array
- `.lambda_` - Lambda values array
- `.spectest` - Specification tests DataFrame
- `.estimator` - Estimator string
- `.ic` - Information criterion string
- `.nbreq` - Number of equations

### Forecast Results

**R:**
```r
fc.lv$err         # Forecast errors
fc.lv$pred        # Predictions
fc.lv$coefficients # List of coefficient matrices
fc.lv$lambda      # Lambda matrix
fc.lv$spectest    # 3D array
```

**Python:**
```python
fc_results['err']          # DataFrame
fc_results['pred']         # DataFrame
fc_results['coefficients'] # List of arrays
fc_results['lambda']       # Array
fc_results['spectest']     # 3D array
```

## Method Differences

### Accessing Components

**R:**
```r
# Access coefficients
coef(lv.mod)  # or lv.mod$coefficients

# Get residuals
residuals(lv.mod)  # or resid(lv.mod)

# Predictions
predict(lv.mod, newdata)

# Summary
summary(lv.mod)
summary(lv.mod, short=TRUE)
```

**Python:**
```python
# Access coefficients
model.coef()  # or model.coefficients

# Get residuals
model.residuals()

# Predictions
model.predict(newdata)

# Summary
model.summary()
model.summary(short=True)
```

## Data Structure Differences

### Input Data

**R:** Accepts data.frame, matrix, or tibble
```r
dat <- data.frame(x1=rnorm(100), x2=rnorm(100))
# or
dat <- matrix(rnorm(200), ncol=2)
```

**Python:** Accepts DataFrame, array, or array-like
```python
dat = pd.DataFrame({'x1': np.random.randn(100), 'x2': np.random.randn(100)})
# or
dat = np.random.randn(100, 2)
```

### Output Data

**R:** Returns lists with matrices and data.frames
**Python:** Returns numpy arrays and pandas DataFrames

## Parallel Processing

**R:**
```r
# Uses mclapply from parallel package
lv.mod <- lassovar(dat, lags=1, mc=TRUE, ncores=4)
```

**Python:**
```python
# Uses joblib
model = LassoVAR(dat, lags=1, mc=True, n_jobs=4)
model.fit()
```

## Common Migration Issues

### 1. Explicit fit() call

**Issue:** In Python, you must call `fit()` after creating the model.

**R:**
```r
lv.mod <- lassovar(dat, lags=1)  # Automatically fitted
```

**Python:**
```python
model = LassoVAR(dat, lags=1)  # Not yet fitted
model.fit()  # Explicit fit required
```

### 2. 1-based vs 0-based indexing

**R:** 1-based indexing
```r
first_coef <- lv.mod$coefficients[1, ]
```

**Python:** 0-based indexing
```python
first_coef = model.coefficients[0, :]
```

### 3. $ vs . for attribute access

**R:**
```r
lv.mod$coefficients
lv.mod$RSS
```

**Python:**
```python
model.coefficients
model.RSS
```

### 4. Function vs Method calls

**R:**
```r
summary(lv.mod)
predict(lv.mod, newdata)
residuals(lv.mod)
```

**Python:**
```python
model.summary()
model.predict(newdata)
model.residuals()
```

### 5. Column name handling

**R:** Automatically handles missing column names
**Python:** Best to provide column names explicitly

```python
# Good practice
dat = pd.DataFrame(data, columns=['V1', 'V2', 'V3'])
```

## Performance Comparison

### Similarities:
- Both use regularized regression (Lasso/Ridge)
- Both select lambda via information criteria
- Both support parallel processing
- Both have similar computational complexity

### Differences:

**R (glmnet):**
- Highly optimized Fortran backend
- Efficient sparse matrix support
- Coordinate descent algorithm

**Python (scikit-learn):**
- Optimized C/Cython backend
- Good performance for dense matrices
- Coordinate descent algorithm

**General guideline:** Performance should be comparable for most use cases. R may be slightly faster for very large sparse problems due to glmnet's optimizations.

## Feature Parity Matrix

| Feature | R | Python | Notes |
|---------|---|--------|-------|
| Basic Lasso VAR | ✓ | ✓ | Full support |
| Adaptive Lasso (OLS) | ✓ | ✓ | Full support |
| Adaptive Lasso (Lasso) | ✓ | ✓ | Full support |
| Adaptive Lasso (Ridge) | ✓ | ✓ | Full support |
| Adaptive Lasso (Group) | ✓ | ~ | Simplified in Python |
| Information Criteria (BIC/AIC) | ✓ | ✓ | Full support |
| Post-Lasso OLS | ✓ | ✓ | Full support |
| Exogenous variables | ✓ | ✓ | Full support |
| Trend | ✓ | ✓ | Full support |
| Multiple lags | ✓ | ✓ | Full support |
| Forecasting | ✓ | ✓ | Full support |
| Fixed window | ✓ | ✓ | Full support |
| Expanding window | ✓ | ✓ | Full support |
| Recursive forecasts | ✓ | ✓ | Full support |
| Direct forecasts | ✓ | ✓ | Full support |
| Parallel processing | ✓ | ✓ | Different backends |
| Specification tests | ✓ | ✓ | Full support |
| Summary method | ✓ | ✓ | Full support |
| Predict method | ✓ | ✓ | Full support |
| Residuals method | ✓ | ✓ | Full support |

Legend: ✓ = Full support, ~ = Partial support, ✗ = Not supported

## Example: Complete Migration

**Original R Code:**
```r
library(lassovar)

# Load data
dat <- read.csv("data.csv")

# Fit adaptive Lasso VAR
lv.mod <- lassovar(
    dat,
    lags=2,
    ic='BIC',
    adaptive='ols',
    post=TRUE,
    mc=TRUE,
    ncores=4
)

# Summary
summary(lv.mod)

# Forecasting
fc.lv <- forecast.lassovar(
    dat,
    fc.train=80,
    horizon=1,
    lags=2,
    fc.window='expanding',
    adaptive='ols',
    silent=FALSE
)

# Analyze errors
mean(abs(fc.lv$err))
```

**Equivalent Python Code:**
```python
import pandas as pd
import numpy as np
from lassovar import LassoVAR, forecast_lassovar

# Load data
dat = pd.read_csv("data.csv")

# Fit adaptive Lasso VAR
model = LassoVAR(
    dat,
    lags=2,
    ic='BIC',
    adaptive='ols',
    post=True,
    mc=True,
    n_jobs=4
)
model.fit()

# Summary
model.summary()

# Forecasting
fc_results = forecast_lassovar(
    dat,
    fc_train=80,
    horizon=1,
    lags=2,
    fc_window='expanding',
    adaptive='ols',
    silent=False
)

# Analyze errors
fc_results['err'].abs().mean()
```

## Tips for Smooth Migration

1. **Start Simple:** Begin with basic Lasso VAR before adding complexity
2. **Check Dimensions:** Verify data shapes match expectations
3. **Use DataFrames:** Provide column names for better tracking
4. **Test Incrementally:** Migrate and test one feature at a time
5. **Compare Results:** Check that coefficients are similar between R and Python
6. **Read Documentation:** Both packages have comprehensive documentation
7. **Use Examples:** Refer to example scripts in both packages

## Getting Help

- **Python package:** See `DOCUMENTATION.md` and `examples/`
- **R package:** See R documentation with `?lassovar`
- **Issues:** Report differences or bugs on GitHub
- **Community:** Ask questions on relevant forums/mailing lists

## Known Differences

1. **Group Lasso:** Python version uses simplified implementation
2. **Numerical Precision:** Minor differences due to different optimization routines
3. **Random Seeds:** Different RNG implementations may give different results
4. **Memory Usage:** Python may use more memory for very large models
5. **Error Messages:** Different phrasing but similar information

## Validation

To verify migration:

1. **Run same data through both:**
```r
# R
lv.mod <- lassovar(dat, lags=1, lambda=0.1)
coef.r <- lv.mod$coefficients
```

```python
# Python
model = LassoVAR(dat, lags=1, lambda_values=[0.1])
model.fit()
coef.py = model.coefficients
```

2. **Compare coefficients:**
```python
# Check correlation
correlation = np.corrcoef(coef.r.flatten(), coef.py.flatten())[0, 1]
print(f"Coefficient correlation: {correlation}")  # Should be > 0.95
```

3. **Compare predictions:**
Compare predictions on same test set to verify similar forecasting performance.
