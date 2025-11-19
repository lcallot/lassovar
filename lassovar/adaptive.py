"""
Adaptive weights computation for adaptive Lasso.

This module contains functions to compute adaptive weights using
various initial estimators (OLS, Lasso, Ridge, Group).
"""

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, LassoCV, RidgeCV
from joblib import Parallel, delayed
import logging

logger = logging.getLogger(__name__)


def compute_ols_weights(y, x, mc=False, n_jobs=-1, gamma=1):
    """
    Compute adaptive weights using OLS as initial estimator.
    
    Parameters
    ----------
    y : pd.DataFrame
        Dependent variables (T x N)
    x : pd.DataFrame or np.ndarray
        Independent variables (T x P)
    mc : bool, optional
        Use parallel computation (default: False)
    n_jobs : int, optional
        Number of parallel jobs (default: -1, all cores)
    gamma : float, optional
        Exponent for adaptive weights (default: 1)
        
    Returns
    -------
    dict
        Dictionary with 'ada' (estimator name), 'b' (coefficients),
        'w' (weights), and 'gamma'
    """
    n_eq = y.shape[1]
    x_array = x.values if isinstance(x, pd.DataFrame) else x
    
    def fit_ols_equation(i):
        """Fit OLS for equation i."""
        yi = y.iloc[:, i].values
        model = LinearRegression(fit_intercept=True)
        model.fit(x_array, yi)
        # Return intercept and coefficients
        return np.concatenate([[model.intercept_], model.coef_])
    
    if mc and n_jobs != 1:
        ols_coefs = Parallel(n_jobs=n_jobs)(
            delayed(fit_ols_equation)(i) for i in range(n_eq)
        )
    else:
        ols_coefs = [fit_ols_equation(i) for i in range(n_eq)]
    
    # Stack coefficients
    ols_b = np.column_stack(ols_coefs)
    
    # Compute adaptive weights (exclude intercept)
    ols_w = np.abs(ols_b[1:, :]) ** (-gamma)
    
    return {
        'ada': 'ols',
        'b': ols_b,
        'w': ols_w,
        'gamma': gamma
    }


def compute_lasso_weights(y, x, ic='BIC', mc=False, n_jobs=-1, 
                          dfmax=None, trend=False, gamma=1):
    """
    Compute adaptive weights using Lasso as initial estimator.
    
    Parameters
    ----------
    y : pd.DataFrame
        Dependent variables (T x N)
    x : pd.DataFrame or np.ndarray
        Independent variables (T x P)
    ic : str, optional
        Information criterion ('BIC' or 'AIC', default: 'BIC')
    mc : bool, optional
        Use parallel computation (default: False)
    n_jobs : int, optional
        Number of parallel jobs (default: -1)
    dfmax : int, optional
        Maximum degrees of freedom (default: None)
    trend : bool, optional
        Whether trend is included (default: False)
    gamma : float, optional
        Exponent for adaptive weights (default: 1)
        
    Returns
    -------
    dict
        Dictionary with adaptive weights
    """
    # Import here to avoid circular dependency
    from .estimation import lassovar_equation
    
    # Estimate Lasso VAR
    lv_las = lassovar_equation(
        y, x, ada_w=None, ic=ic, mc=mc, n_jobs=n_jobs,
        alpha=1.0, dfmax=dfmax, trend=trend, lambda_values=None
    )
    
    # Compute adaptive weights from Lasso coefficients
    ada_w = {
        'ada': 'lasso',
        'b': lv_las['coefficients'],
        'w': np.abs(lv_las['coefficients'][1:, :]) ** (-gamma),
        'gamma': gamma
    }
    
    return ada_w


def compute_ridge_weights(y, x, ic='BIC', mc=False, n_jobs=-1,
                          dfmax=None, trend=False, gamma=1):
    """
    Compute adaptive weights using Ridge regression as initial estimator.
    
    Parameters
    ----------
    y : pd.DataFrame
        Dependent variables (T x N)
    x : pd.DataFrame or np.ndarray
        Independent variables (T x P)
    ic : str, optional
        Information criterion ('BIC' or 'AIC', default: 'BIC')
    mc : bool, optional
        Use parallel computation (default: False)
    n_jobs : int, optional
        Number of parallel jobs (default: -1)
    dfmax : int, optional
        Maximum degrees of freedom (default: None)
    trend : bool, optional
        Whether trend is included (default: False)
    gamma : float, optional
        Exponent for adaptive weights (default: 1)
        
    Returns
    -------
    dict
        Dictionary with adaptive weights
    """
    # Import here to avoid circular dependency
    from .estimation import lassovar_equation
    
    # Estimate Ridge VAR (alpha=0)
    lv_ridge = lassovar_equation(
        y, x, ada_w=None, ic=ic, mc=mc, n_jobs=n_jobs,
        alpha=0.0, dfmax=dfmax, trend=trend, lambda_values=None
    )
    
    # Compute adaptive weights from Ridge coefficients
    ada_w = {
        'ada': 'ridge',
        'b': lv_ridge['coefficients'],
        'w': np.abs(lv_ridge['coefficients'][1:, :]) ** (-gamma),
        'gamma': gamma
    }
    
    return ada_w


def compute_group_weights(y, x, trend=False, gamma=1):
    """
    Compute adaptive weights using Group Lasso as initial estimator.
    
    Note: This is a simplified implementation. Full group Lasso
    would require additional dependencies.
    
    Parameters
    ----------
    y : pd.DataFrame
        Dependent variables (T x N)
    x : pd.DataFrame or np.ndarray
        Independent variables (T x P)
    trend : bool, optional
        Whether trend is included (default: False)
    gamma : float, optional
        Exponent for adaptive weights (default: 1)
        
    Returns
    -------
    dict
        Dictionary with adaptive weights
    """
    logger.info("Group Lasso computation")
    logger.info("Note: Simplified implementation without group structure")
    
    # For now, use standard Lasso as approximation
    # A full implementation would require group Lasso optimization
    n_eq = y.shape[1]
    x_array = x.values if isinstance(x, pd.DataFrame) else x
    
    greg_coef = []
    
    for i in range(n_eq):
        logger.debug(f"Processing equation: {i}")
        yi = y.iloc[:, i].values
        
        # Use LassoCV as approximation
        model = LassoCV(cv=5, fit_intercept=True, max_iter=500)
        model.fit(x_array, yi)
        
        # Get coefficients
        coef = np.concatenate([[model.intercept_], model.coef_])
        greg_coef.append(coef[1:])  # Exclude intercept for weights
    
    logger.info("Group lasso computation complete")
    
    greg_coef = np.column_stack(greg_coef)
    
    grp_w = {
        'ada': 'group',
        'b': np.vstack([np.zeros(n_eq), greg_coef]),  # Add intercept row
        'w': np.abs(greg_coef) ** (-gamma),
        'gamma': gamma
    }
    
    return grp_w


def get_adaptive_weights(y, x, adaptive_type, ic='BIC', mc=False, 
                        n_jobs=-1, dfmax=None, trend=False):
    """
    Compute adaptive weights based on specified estimator type.
    
    Parameters
    ----------
    y : pd.DataFrame
        Dependent variables (T x N)
    x : pd.DataFrame or np.ndarray
        Independent variables (T x P)
    adaptive_type : str
        Type of initial estimator ('ols', 'lasso', 'ridge', 'group')
    ic : str, optional
        Information criterion (default: 'BIC')
    mc : bool, optional
        Use parallel computation (default: False)
    n_jobs : int, optional
        Number of parallel jobs (default: -1)
    dfmax : int, optional
        Maximum degrees of freedom (default: None)
    trend : bool, optional
        Whether trend is included (default: False)
        
    Returns
    -------
    dict or None
        Dictionary with adaptive weights, or None if adaptive_type is 'none'
    """
    if adaptive_type == 'none':
        return None
    
    logger.info(f"Initial estimator for the adaptive lasso: {adaptive_type}")
    
    if adaptive_type == 'ols':
        return compute_ols_weights(y, x, mc=mc, n_jobs=n_jobs)
    elif adaptive_type == 'lasso':
        return compute_lasso_weights(y, x, ic=ic, mc=mc, n_jobs=n_jobs,
                                     dfmax=dfmax, trend=trend)
    elif adaptive_type == 'ridge':
        return compute_ridge_weights(y, x, ic=ic, mc=mc, n_jobs=n_jobs,
                                      dfmax=dfmax, trend=trend)
    elif adaptive_type == 'group':
        return compute_group_weights(y, x, trend=trend)
    else:
        raise ValueError(f"Unknown adaptive type: {adaptive_type}")
