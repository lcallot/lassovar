"""
Main LassoVAR class for Vector Autoregression with Lasso.

This module provides the main interface for fitting VAR models
using Lasso and adaptive Lasso estimators.
"""

import numpy as np
import pandas as pd
import warnings

from .utils import (
    make_var_data, coerce_to_dataframe, validate_lags
)
from .adaptive import get_adaptive_weights
from .estimation import lassovar_equation, post_ols_estimation


class LassoVAR:
    """
    Fits a Vector Autoregressive model by L1 penalization.
    
    This class provides functionality for estimating VAR models using
    Lasso or adaptive Lasso, with penalty parameter selection via
    information criteria (BIC or AIC).
    
    Parameters
    ----------
    dat : pd.DataFrame or array-like
        Data containing the time series (T x N)
    exo : pd.DataFrame or array-like, optional
        Exogenous variables (T x K). Not lagged. (default: None)
    lags : int, optional
        Number of lags to include (default: 1)
    ic : str, optional
        Information criterion ('BIC' or 'AIC', default: 'BIC')
    adaptive : str, optional
        Initial estimator for adaptive Lasso: 'none', 'ols', 'lasso',
        'ridge', 'group' (default: 'none' for plain Lasso)
    post : bool, optional
        Estimate post-Lasso OLS (default: False)
    mc : bool, optional
        Enable parallel processing (default: False)
    n_jobs : int, optional
        Number of parallel jobs, -1 uses all cores (default: -1)
    dfmax : int, optional
        Maximum number of variables excluding intercept (default: None)
    horizon : int, optional
        Horizon for h-step ahead estimation (default: 1)
    trend : bool, optional
        Include linear trend (default: False)
    lambda_values : array-like, optional
        User-defined lambda values (default: None)
        
    Attributes
    ----------
    var_names : list
        Names of endogenous variables
    ada_w : dict or None
        Adaptive weights (if adaptive != 'none')
    x : pd.DataFrame
        Right-hand side variables
    y : pd.DataFrame
        Left-hand side variables
    coefficients : np.ndarray
        Estimated coefficients (P+1 x N)
    RSS : np.ndarray
        Residual sum of squares for each equation
    lambda_ : np.ndarray
        Selected lambda values for each equation
    spectest : pd.DataFrame
        Specification tests results
    estimator : str
        Name of estimator used
    ic : str
        Information criterion used
    nbreq : int
        Number of equations
    post : np.ndarray or None
        Post-Lasso OLS coefficients (if post=True)
    fitted_ : bool
        Whether model has been fitted
        
    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> from lassovar import LassoVAR
    >>> 
    >>> # Generate sample data
    >>> data = pd.DataFrame(np.random.randn(100, 3), 
    ...                     columns=['V1', 'V2', 'V3'])
    >>> 
    >>> # Fit basic Lasso VAR
    >>> model = LassoVAR(data, lags=1)
    >>> model.fit()
    >>> 
    >>> # Fit adaptive Lasso VAR
    >>> model_ada = LassoVAR(data, lags=2, adaptive='ols')
    >>> model_ada.fit()
    >>> 
    >>> # Get summary
    >>> model.summary()
    >>> 
    >>> # Get predictions
    >>> predictions = model.predict(model.x.iloc[-1:])
    """
    
    def __init__(self, dat, exo=None, lags=1, ic='BIC', adaptive='none',
                 post=False, mc=False, n_jobs=-1, dfmax=None, horizon=1,
                 trend=False, lambda_values=None):
        """Initialize LassoVAR model."""
        # Validate and store parameters
        self.ic = ic.upper() if ic.upper() in ['BIC', 'AIC'] else 'BIC'
        self.adaptive = adaptive.lower()
        if self.adaptive not in ['none', 'ols', 'lasso', 'ridge', 'group']:
            raise ValueError(f"Invalid adaptive type: {adaptive}")
        
        self.post_ols = post
        self.mc = mc
        self.n_jobs = n_jobs
        self.dfmax = dfmax
        self.horizon = horizon
        self.trend = trend
        self.lambda_values = lambda_values
        
        # Validate lags
        self.lags = validate_lags(lags)
        
        # Coerce data to DataFrames
        self.dat = coerce_to_dataframe(dat)
        self.exo = coerce_to_dataframe(exo) if exo is not None else None
        
        # Check dfmax
        if self.dfmax is None:
            self.dfmax = self.dat.shape[1] * self.lags
        else:
            self.dfmax = int(self.dfmax)
        
        # Prepare VAR data
        y_var = make_var_data(
            self.dat, 
            lags=self.lags,
            horizon=self.horizon,
            exo=self.exo,
            trend=self.trend
        )
        
        self.x = y_var['x']
        self.y = y_var['y']
        self.var_names = list(self.y.columns)
        
        # Model results (to be filled by fit())
        self.coefficients = None
        self.RSS = None
        self.lambda_ = None
        self.spectest = None
        self.estimator = None
        self.nbreq = self.dat.shape[1]
        self.ada_w = None
        self.post = None
        self.fitted_ = False
        
        # Store call information
        self.call = {
            'lags': lags,
            'ic': self.ic,
            'adaptive': self.adaptive,
            'post': self.post_ols,
            'mc': self.mc,
            'dfmax': self.dfmax,
            'horizon': self.horizon,
            'trend': self.trend
        }
    
    def fit(self):
        """
        Fit the Lasso VAR model.
        
        Returns
        -------
        self
            Fitted model
        """
        # Compute adaptive weights if needed
        if self.adaptive != 'none':
            self.ada_w = get_adaptive_weights(
                self.y, self.x, self.adaptive,
                ic=self.ic, mc=self.mc, n_jobs=self.n_jobs,
                dfmax=self.dfmax, trend=self.trend
            )
        
        # Estimate the model
        las_mod = lassovar_equation(
            self.y, self.x, self.ada_w,
            ic=self.ic, mc=self.mc, n_jobs=self.n_jobs,
            alpha=1.0, dfmax=self.dfmax, trend=self.trend,
            lambda_values=self.lambda_values
        )
        
        # Store results
        self.coefficients = las_mod['coefficients']
        self.RSS = las_mod['RSS']
        self.lambda_ = las_mod['lambda']
        self.spectest = las_mod['spectest']
        self.estimator = las_mod['estimator']
        
        # Post-Lasso OLS if requested
        if self.post_ols:
            selected = self.coefficients != 0
            self.post = post_ols_estimation(
                self.y, self.x, selected,
                mc=self.mc, n_jobs=self.n_jobs
            )
        
        self.fitted_ = True
        return self
    
    def predict(self, newdata):
        """
        Make predictions with new data.
        
        Parameters
        ----------
        newdata : pd.DataFrame, np.ndarray, or array-like
            New data for prediction. Should have same number of features
            as training data (excluding intercept).
            
        Returns
        -------
        np.ndarray or pd.DataFrame
            Predictions for each equation
        """
        if not self.fitted_:
            raise RuntimeError("Model must be fitted before prediction. Call fit() first.")
        
        # Convert to array
        if isinstance(newdata, pd.DataFrame):
            newdata_array = newdata.values
        else:
            newdata_array = np.asarray(newdata)
        
        # Handle different dimensions
        if newdata_array.ndim == 1:
            newdata_array = newdata_array.reshape(1, -1)
        
        # Add intercept
        newdata_with_intercept = np.column_stack([
            np.ones(newdata_array.shape[0]),
            newdata_array
        ])
        
        # Make predictions
        predictions = newdata_with_intercept @ self.coefficients
        
        # Return as DataFrame with proper column names
        if isinstance(newdata, pd.DataFrame) or len(predictions) > 1:
            return pd.DataFrame(predictions, columns=self.var_names)
        else:
            return predictions
    
    def residuals(self):
        """
        Extract model residuals.
        
        Returns
        -------
        pd.DataFrame
            Residuals for each equation
        """
        if not self.fitted_:
            raise RuntimeError("Model must be fitted before extracting residuals. Call fit() first.")
        
        # Compute fitted values
        x_with_intercept = np.column_stack([np.ones(len(self.x)), self.x.values])
        fitted = x_with_intercept @ self.coefficients
        
        # Compute residuals
        residuals = self.y.values - fitted
        
        return pd.DataFrame(residuals, columns=self.var_names, index=self.y.index)
    
    def coef(self):
        """
        Get coefficient matrix.
        
        Returns
        -------
        np.ndarray
            Coefficient matrix (P+1 x N), where first row is intercepts
        """
        if not self.fitted_:
            raise RuntimeError("Model must be fitted first. Call fit().")
        return self.coefficients
    
    def summary(self, short=False):
        """
        Print summary statistics of the fitted model.
        
        Parameters
        ----------
        short : bool, optional
            Print shortened summary (default: False)
            
        Returns
        -------
        pd.DataFrame or None
            Summary statistics table (if not short), None otherwise
        """
        if not self.fitted_:
            raise RuntimeError("Model must be fitted first. Call fit().")
        
        print('Call:')
        print(f"  LassoVAR(lags={self.call['lags']}, ic='{self.call['ic']}', " +
              f"adaptive='{self.call['adaptive']}', trend={self.call['trend']})")
        print()
        
        print('Model estimated equation by equation')
        print(f'Selection criterion: {self.ic}')
        print(f'Estimator: {self.estimator}')
        
        if self.ada_w is not None:
            print(f"Adaptive weights Estimator: {self.ada_w['ada']}")
        
        deterministics = 'intercept and trend' if self.trend else 'intercept'
        print(f'Deterministics: {deterministics}.')
        
        T = len(self.y)
        N = len(self.var_names)
        print(f'Dimensions: T = {T}  N = {N}')
        
        # Count non-zero coefficients
        if self.trend:
            # Exclude intercept and trend
            nzeq = np.sum(self.coefficients[1:-1, :] != 0, axis=0)
        else:
            # Exclude only intercept
            nzeq = np.sum(self.coefficients[1:, :] != 0, axis=0)
        
        nbr_nz = np.sum(nzeq)
        total_candidates = self.coefficients[1:, :].size
        
        print(f'\nTotal number of variables selected: {nbr_nz} ' +
              f'({100*nbr_nz/total_candidates:.1f}% of candidates)')
        
        # Create summary table
        eq_sum = pd.DataFrame({
            'Lambda': self.lambda_,
            'non-zero': nzeq,
            'resid var': self.RSS / T
        }, index=self.var_names)
        
        if not short:
            print('\nModel summary statistics:')
            print(eq_sum)
            print('\nDeterministics not included in the non-zero count.')
        else:
            print('Summary statistics')
            print('Average for all equations:')
            print(eq_sum.mean())
        
        if not short:
            return eq_sum
    
    def __repr__(self):
        """String representation."""
        if self.fitted_:
            return (f"LassoVAR(fitted=True, estimator='{self.estimator}', " +
                   f"n_equations={self.nbreq}, ic='{self.ic}')")
        else:
            return (f"LassoVAR(fitted=False, lags={self.lags}, " +
                   f"n_series={self.dat.shape[1]})")
