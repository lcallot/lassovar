"""
Utility functions for lassovar package.

This module contains helper functions for data preparation,
specification tests, and other utilities.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.diagnostic import acorr_ljungbox


def make_var_data(data_var, lags, horizon=1, exo=None, trend=False):
    """
    Prepare data for VAR estimation by creating lagged variables.
    
    Parameters
    ----------
    data_var : pd.DataFrame
        Input data with time series
    lags : int
        Number of lags to include
    horizon : int, optional
        Forecast horizon for h-step ahead models (default: 1)
    exo : pd.DataFrame, optional
        Exogenous variables (default: None)
    trend : bool, optional
        Include linear trend (default: False)
        
    Returns
    -------
    dict
        Dictionary with 'y' (dependent variables), 'x' (independent variables),
        'lags', 'horizon', 'nbrser', and 'trend'
    """
    nbrser = data_var.shape[1]
    
    # The dependent variable
    y = data_var.copy()
    if data_var.columns is not None:
        varny = list(data_var.columns)
    else:
        varny = [f'Eq_{i}' for i in range(nbrser)]
    y.columns = varny
    
    # Create lagged variables
    x_list = []
    varnx = []
    
    for l in range(1, lags + 1):
        # Create lag by shifting and padding with NaN
        lagged = data_var.shift(l + horizon - 1)
        x_list.append(lagged)
        varnx.extend([f'{l}L_{v}' for v in varny])
    
    x = pd.concat(x_list, axis=1)
    
    # Add exogenous variables
    if exo is not None:
        if exo.columns is not None:
            vnexo = list(exo.columns)
        else:
            vnexo = [f'Exo_{i}' for i in range(exo.shape[1])]
        varnx.extend(vnexo)
        x = pd.concat([x, exo], axis=1)
    
    # Trim to remove NaN values
    trim_len = lags + horizon - 1
    y = y.iloc[trim_len:]
    x = x.iloc[trim_len:]
    
    # Add trend if requested
    if trend:
        x['trend'] = np.arange(1, len(x) + 1)
        varnx.append('trend')
    
    x.columns = varnx
    
    # Reset indices
    y = y.reset_index(drop=True)
    x = x.reset_index(drop=True)
    
    return {
        'y': y,
        'x': x,
        'lags': lags,
        'horizon': horizon,
        'nbrser': nbrser,
        'trend': trend
    }


def specification_tests(residuals, y_original, lags=None, fitdf=1):
    """
    Perform specification tests on residuals.
    
    Parameters
    ----------
    residuals : np.ndarray or pd.Series
        Model residuals
    y_original : np.ndarray or pd.Series
        Original dependent variable
    lags : int, optional
        Number of lags for Ljung-Box test (default: min(floor(3*fitdf), n/2))
    fitdf : int, optional
        Degrees of freedom for the fit (default: 1)
        
    Returns
    -------
    dict
        Dictionary with test results: 'Ljung-Box' p-value, 
        'Shapiro' p-value (if n <= 5000), and 'R2'
    """
    residuals = np.asarray(residuals).flatten()
    y_original = np.asarray(y_original).flatten()
    
    n = len(residuals)
    
    if lags is None:
        lags = min(int(np.floor(3 * fitdf)), n // 2)
    
    sptest = {}
    
    # Ljung-Box test for autocorrelation
    try:
        lb_result = acorr_ljungbox(residuals, lags=[lags], return_df=False)
        sptest['Ljung-Box'] = float(lb_result[1][0])  # p-value
    except Exception:
        sptest['Ljung-Box'] = np.nan
    
    # Shapiro-Wilk test for normality (only if n <= 5000)
    if n <= 5000:
        try:
            _, p_value = stats.shapiro(residuals)
            sptest['Shapiro'] = p_value
        except Exception:
            sptest['Shapiro'] = np.nan
    
    # R-squared
    try:
        var_res = np.var(residuals, ddof=1)
        var_y = np.var(y_original, ddof=1)
        r2 = 1 - (n * var_res) / (var_y * len(y_original))
        sptest['R2'] = r2
    except Exception:
        sptest['R2'] = np.nan
    
    return sptest


def compute_bic(rss, n_obs, df):
    """
    Compute Bayesian Information Criterion.
    
    Parameters
    ----------
    rss : float or np.ndarray
        Residual sum of squares
    n_obs : int
        Number of observations
    df : int or np.ndarray
        Degrees of freedom
        
    Returns
    -------
    float or np.ndarray
        BIC value(s)
    """
    return np.log(rss / n_obs) + df * np.log(n_obs) / n_obs


def compute_aic(rss, n_obs, df):
    """
    Compute Akaike Information Criterion.
    
    Parameters
    ----------
    rss : float or np.ndarray
        Residual sum of squares
    n_obs : int
        Number of observations
    df : int or np.ndarray
        Degrees of freedom
        
    Returns
    -------
    float or np.ndarray
        AIC value(s)
    """
    return np.log(rss / n_obs) + df / n_obs


def ridge_degrees_of_freedom(x, lambda_values):
    """
    Compute degrees of freedom for ridge regression.
    
    Parameters
    ----------
    x : np.ndarray
        Design matrix
    lambda_values : float or array-like
        Ridge penalty parameter(s)
        
    Returns
    -------
    float or np.ndarray
        Degrees of freedom for each lambda value
    """
    # Add intercept column
    x_with_intercept = np.column_stack([np.ones(x.shape[0]), x])
    
    if np.isscalar(lambda_values):
        lambda_values = [lambda_values]
    
    df_ridge = []
    for lam in lambda_values:
        try:
            xtx = x_with_intercept.T @ x_with_intercept
            ridge_matrix = np.linalg.solve(
                xtx + lam * np.eye(xtx.shape[0]),
                xtx
            )
            hat_matrix = x_with_intercept @ ridge_matrix @ x_with_intercept.T
            df = np.trace(hat_matrix)
            df_ridge.append(df)
        except np.linalg.LinAlgError:
            df_ridge.append(np.nan)
    
    return np.array(df_ridge) if len(df_ridge) > 1 else df_ridge[0]


def coerce_to_dataframe(data):
    """
    Coerce input to pandas DataFrame.
    
    Parameters
    ----------
    data : various types
        Input data (DataFrame, array, list, etc.)
        
    Returns
    -------
    pd.DataFrame
        Data as DataFrame
    """
    if isinstance(data, pd.DataFrame):
        return data
    elif isinstance(data, pd.Series):
        return data.to_frame()
    elif isinstance(data, np.ndarray):
        return pd.DataFrame(data)
    elif isinstance(data, (list, tuple)):
        return pd.DataFrame(data)
    else:
        try:
            return pd.DataFrame(data)
        except Exception:
            raise ValueError(f"Cannot convert {type(data)} to DataFrame")


def validate_lags(lags):
    """
    Validate lags parameter.
    
    Parameters
    ----------
    lags : int
        Number of lags
        
    Returns
    -------
    int
        Validated lags value
        
    Raises
    ------
    ValueError
        If lags is not a positive integer
    """
    if not isinstance(lags, (int, np.integer)) or lags <= 0:
        raise ValueError("'lags' must be a positive integer")
    return int(lags)
