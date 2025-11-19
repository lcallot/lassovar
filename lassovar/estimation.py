"""
Core estimation functions for lassovar.

This module contains the main estimation routines for
Lasso and adaptive Lasso VAR models.
"""

import numpy as np
import pandas as pd
from sklearn.linear_model import Lasso, Ridge, ElasticNet, LinearRegression
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
import warnings

from .utils import (
    compute_bic, compute_aic, specification_tests,
    ridge_degrees_of_freedom
)


def fit_glmnet_equation(i, y, x, ada_w, ic, alpha, dfmax, trend, lambda_values):
    """
    Fit a single equation using elastic net (Lasso/Ridge).
    
    Parameters
    ----------
    i : int
        Equation index
    y : pd.DataFrame
        Dependent variables
    x : pd.DataFrame or np.ndarray
        Independent variables
    ada_w : dict or None
        Adaptive weights
    ic : str
        Information criterion ('BIC' or 'AIC')
    alpha : float
        ElasticNet mixing parameter (1.0 = Lasso, 0.0 = Ridge)
    dfmax : int or None
        Maximum degrees of freedom
    trend : bool
        Whether trend is included
    lambda_values : array-like or None
        User-specified lambda values
        
    Returns
    -------
    dict
        Results for this equation
    """
    x_array = x.values if isinstance(x, pd.DataFrame) else x
    yi = y.iloc[:, i].values
    n_obs = len(yi)
    n_features = x_array.shape[1]
    
    # Handle adaptive weights
    if ada_w is None:
        # Plain Lasso/Ridge
        penalty_weights = np.ones(n_features)
        if trend:
            penalty_weights[-1] = 0  # No penalty for trend
        excluded_vars = []
    else:
        # Adaptive Lasso/Ridge
        weights_i = ada_w['w'][:, i]
        
        # Check if all weights are infinite (empty model)
        if np.all(np.isinf(weights_i)):
            warnings.warn('Adaptive weights all equal to Inf. Returning empty model.')
            return {
                'coefficients': np.concatenate([[yi.mean()], np.zeros(n_features)]),
                'rss': np.sum((yi - yi.mean())**2),
                'lambda': 0.0,
                'spectest': specification_tests(yi - yi.mean(), yi),
                'ic': ic,
                'ic_all': None
            }
        
        penalty_weights = weights_i.copy()
        # Variables with infinite weights are excluded
        excluded_vars = np.where(np.isinf(weights_i))[0]
        penalty_weights[np.isinf(penalty_weights)] = 0
        
        if trend:
            penalty_weights[-1] = 0  # No penalty for trend
    
    # Generate lambda path if not provided
    if lambda_values is None:
        # Compute lambda_max
        if ada_w is None or alpha == 1.0:
            # For Lasso, standardize features
            scaler = StandardScaler(with_mean=True, with_std=True)
            x_scaled = scaler.fit_transform(x_array)
            y_scaled = yi - yi.mean()
            
            # Compute lambda_max
            lambda_max = np.max(np.abs(x_scaled.T @ y_scaled)) / (n_obs * alpha) if alpha > 0 else 1.0
        else:
            lambda_max = 1.0
        
        # Create lambda path
        n_lambdas = 100
        lambda_min_ratio = 0.0001
        lambda_values = np.logspace(
            np.log10(lambda_max),
            np.log10(lambda_max * lambda_min_ratio),
            n_lambdas
        )
    
    # Fit models for different lambda values
    coefficients_list = []
    intercepts_list = []
    lambda_used = []
    
    for lam in lambda_values:
        if alpha == 1.0:
            # Lasso
            if ada_w is None:
                model = Lasso(alpha=lam, fit_intercept=True, max_iter=10000, 
                             tol=1e-4, warm_start=False)
            else:
                # Adaptive Lasso: modify penalty
                # sklearn doesn't support per-feature penalties directly,
                # so we scale features
                x_weighted = x_array.copy()
                for j in range(n_features):
                    if j not in excluded_vars and penalty_weights[j] > 0:
                        x_weighted[:, j] = x_weighted[:, j] / penalty_weights[j]
                
                model = Lasso(alpha=lam, fit_intercept=True, max_iter=10000,
                             tol=1e-4, warm_start=False)
                try:
                    model.fit(x_weighted, yi)
                    # Rescale coefficients
                    coef = model.coef_ / np.where(penalty_weights > 0, penalty_weights, 1)
                    coef[excluded_vars] = 0
                    coefficients_list.append(coef)
                    intercepts_list.append(model.intercept_)
                    lambda_used.append(lam)
                    continue
                except:
                    continue
        elif alpha == 0.0:
            # Ridge
            model = Ridge(alpha=lam, fit_intercept=True, max_iter=10000, tol=1e-4)
        else:
            # Elastic Net
            model = ElasticNet(alpha=lam, l1_ratio=alpha, fit_intercept=True,
                              max_iter=10000, tol=1e-4, warm_start=False)
        
        try:
            if alpha != 1.0 or ada_w is None:
                model.fit(x_array, yi)
                coefficients_list.append(model.coef_)
                intercepts_list.append(model.intercept_)
                lambda_used.append(lam)
        except:
            continue
        
        # Check dfmax constraint
        if dfmax is not None:
            n_nonzero = np.sum(model.coef_ != 0)
            if n_nonzero > dfmax:
                break
    
    if len(coefficients_list) == 0:
        # If no models fitted, return empty model
        return {
            'coefficients': np.concatenate([[yi.mean()], np.zeros(n_features)]),
            'rss': np.sum((yi - yi.mean())**2),
            'lambda': 0.0,
            'spectest': specification_tests(yi - yi.mean(), yi),
            'ic': ic,
            'ic_all': None
        }
    
    # Stack results
    coefficients = np.array(coefficients_list)  # shape: (n_lambdas, n_features)
    intercepts = np.array(intercepts_list)      # shape: (n_lambdas,)
    lambda_used = np.array(lambda_used)
    
    # Compute predictions and residuals
    predictions = intercepts.reshape(-1, 1) + x_array @ coefficients.T
    residuals = yi.reshape(-1, 1) - predictions
    
    # Compute RSS
    rss = np.sum(residuals**2, axis=0)
    
    # Compute degrees of freedom
    if alpha == 1.0:
        # For Lasso, df = number of non-zero coefficients
        df = np.sum(coefficients != 0, axis=1) + 1  # +1 for intercept
    elif alpha == 0.0:
        # For Ridge, compute effective df
        df = ridge_degrees_of_freedom(x_array, lambda_used) + 1
    else:
        # For Elastic Net, approximate df
        df = np.sum(coefficients != 0, axis=1) + 1
    
    # Compute information criterion
    if ic == 'BIC':
        ic_values = compute_bic(rss, n_obs, df)
    else:  # AIC
        ic_values = compute_aic(rss, n_obs, df)
    
    # Select best model
    best_idx = np.argmin(ic_values)
    best_coef = np.concatenate([[intercepts[best_idx]], coefficients[best_idx]])
    best_rss = rss[best_idx]
    best_lambda = lambda_used[best_idx]
    
    # Specification tests
    best_residuals = residuals[:, best_idx]
    spectest = specification_tests(best_residuals, yi)
    
    return {
        'coefficients': best_coef,
        'rss': best_rss,
        'lambda': best_lambda,
        'spectest': spectest,
        'ic': ic,
        'ic_all': ic_values
    }


def lassovar_equation(y, x, ada_w, ic='BIC', mc=False, n_jobs=-1,
                      alpha=1.0, dfmax=None, trend=False, lambda_values=None):
    """
    Estimate VAR equations using Lasso or Ridge.
    
    Parameters
    ----------
    y : pd.DataFrame
        Dependent variables (T x N)
    x : pd.DataFrame or np.ndarray
        Independent variables (T x P)
    ada_w : dict or None
        Adaptive weights from initial estimator
    ic : str, optional
        Information criterion ('BIC' or 'AIC', default: 'BIC')
    mc : bool, optional
        Use parallel computation (default: False)
    n_jobs : int, optional
        Number of parallel jobs (default: -1)
    alpha : float, optional
        ElasticNet mixing parameter (default: 1.0 for Lasso)
    dfmax : int, optional
        Maximum degrees of freedom (default: None)
    trend : bool, optional
        Whether trend is included (default: False)
    lambda_values : array-like, optional
        User-specified lambda values (default: None)
        
    Returns
    -------
    dict
        Estimation results with coefficients, RSS, lambda, etc.
    """
    n_eq = y.shape[1]
    var_names = list(y.columns)
    
    # Estimate equation by equation
    if mc and n_jobs != 1:
        results = Parallel(n_jobs=n_jobs)(
            delayed(fit_glmnet_equation)(
                i, y, x, ada_w, ic, alpha, dfmax, trend, lambda_values
            ) for i in range(n_eq)
        )
    else:
        results = [
            fit_glmnet_equation(
                i, y, x, ada_w, ic, alpha, dfmax, trend, lambda_values
            ) for i in range(n_eq)
        ]
    
    # Aggregate results
    coefficients = np.column_stack([r['coefficients'] for r in results])
    rss = np.array([r['rss'] for r in results])
    lambdas = np.array([r['lambda'] for r in results])
    
    # Specification tests
    spectest_names = list(results[0]['spectest'].keys())
    spectest = np.zeros((len(spectest_names), n_eq))
    for i, r in enumerate(results):
        for j, name in enumerate(spectest_names):
            spectest[j, i] = r['spectest'][name]
    
    spectest_df = pd.DataFrame(
        spectest,
        index=spectest_names,
        columns=var_names
    )
    
    # Determine estimator name
    if ada_w is None:
        estimator = 'Lasso' if alpha == 1.0 else 'Ridge' if alpha == 0.0 else 'ElasticNet'
    else:
        estimator = 'Adaptive Lasso' if alpha == 1.0 else 'Adaptive Ridge'
    
    return {
        'var_names': var_names,
        'ada_w': ada_w,
        'x': x,
        'y': y,
        'coefficients': coefficients,
        'RSS': rss,
        'lambda': lambdas,
        'spectest': spectest_df,
        'estimator': estimator,
        'ic': ic,
        'nbreq': n_eq,
        'trend': trend
    }


def post_ols_estimation(y, x, selected_params, mc=False, n_jobs=-1):
    """
    Perform post-Lasso OLS estimation on selected variables.
    
    Parameters
    ----------
    y : pd.DataFrame
        Dependent variables (T x N)
    x : pd.DataFrame or np.ndarray
        Independent variables (T x P)
    selected_params : np.ndarray
        Boolean array indicating selected parameters (excluding intercept)
    mc : bool, optional
        Use parallel computation (default: False)
    n_jobs : int, optional
        Number of parallel jobs (default: -1)
        
    Returns
    -------
    np.ndarray
        Post-Lasso OLS coefficients
    """
    n_eq = y.shape[1]
    n_features = x.shape[1]
    x_array = x.values if isinstance(x, pd.DataFrame) else x
    
    def post_ols_eq(i):
        """Estimate post-OLS for equation i."""
        yi = y.iloc[:, i].values
        selected_i = selected_params[1:, i]  # Exclude intercept
        
        # Initialize with zeros
        coef = np.zeros(n_features + 1)
        
        if np.sum(selected_i) > 0:
            # Fit OLS on selected variables
            x_selected = x_array[:, selected_i]
            model = LinearRegression(fit_intercept=True)
            model.fit(x_selected, yi)
            
            # Place coefficients
            coef[0] = model.intercept_
            coef[1:][selected_i] = model.coef_
        else:
            # No variables selected, use mean
            coef[0] = yi.mean()
        
        return coef
    
    if mc and n_jobs != 1:
        post_coefs = Parallel(n_jobs=n_jobs)(
            delayed(post_ols_eq)(i) for i in range(n_eq)
        )
    else:
        post_coefs = [post_ols_eq(i) for i in range(n_eq)]
    
    return np.column_stack(post_coefs)
