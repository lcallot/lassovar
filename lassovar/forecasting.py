"""
Forecasting functionality for lassovar.

This module provides functions for pseudo out-of-sample forecasting
experiments with Lasso VAR models.
"""

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import time

from .lassovar import LassoVAR
from .utils import coerce_to_dataframe, make_var_data


def forecast_loop_iteration(fc, dat, fc_train, fc_window, lags, horizon, 
                            ic, exo, adaptive, post, trend, silent, mc_eq, n_jobs):
    """
    Single iteration of the forecasting loop.
    
    Parameters
    ----------
    fc : int
        Forecast iteration number
    dat : pd.DataFrame
        Full dataset
    fc_train : int
        Number of training observations
    fc_window : str
        'fix' or 'expanding'
    lags : int
        Number of lags
    horizon : int
        Forecast horizon
    ic : str
        Information criterion
    exo : pd.DataFrame or None
        Exogenous variables
    adaptive : str
        Adaptive type
    post : bool
        Post-Lasso OLS
    trend : bool
        Include trend
    silent : bool
        Suppress output
    mc_eq : bool
        Parallelize equations
    n_jobs : int
        Number of jobs
        
    Returns
    -------
    dict
        Forecast results for this iteration
    """
    # Determine training sample range
    start_train = 1 + (fc if fc_window == 'fix' else 0)
    end_train = start_train + fc_train - 1 + (fc if fc_window == 'expanding' else 0)
    
    # Extract training data
    train_dat = dat.iloc[start_train:end_train + 1]
    if exo is not None:
        train_exo = exo.iloc[start_train:end_train + 1]
    else:
        train_exo = None
    
    # True forecast value
    fc_dat = dat.iloc[end_train + 1]
    
    # Estimation
    fc_time = time.time()
    
    # Prepare VAR data
    y_var = make_var_data(train_dat, lags, horizon, train_exo, trend)
    
    # Fit model
    model = LassoVAR(
        train_dat,
        exo=train_exo,
        lags=lags,
        ic=ic,
        adaptive=adaptive,
        post=post,
        mc=mc_eq,
        n_jobs=n_jobs if mc_eq else 1,
        horizon=horizon,
        trend=trend
    )
    model.fit()
    
    # Make prediction
    last_x = y_var['x'].iloc[-1:]
    pred = model.predict(last_x)
    
    if isinstance(pred, pd.DataFrame):
        pred = pred.values[0]
    
    # Forecast error
    fcerr = fc_dat.values - pred
    
    # Store results
    result = {
        'pred': pred,
        'fcerr': fcerr,
        'coefficients': model.coefficients,
        'lambda': model.lambda_,
        'spectest': model.spectest.values
    }
    
    if post:
        result['post'] = model.post
    
    elapsed = time.time() - fc_time
    if not silent:
        print(f'fc {fc} completed in {elapsed/60:.2f} minutes.')
    
    return result


def forecast_lassovar(dat, exo=None, fc_train=None, horizon=1, lags=1,
                     fc_window='fix', fc_type='recursive', ic='BIC',
                     adaptive='none', mc=False, n_jobs=-1, silent=False,
                     trend=False, post=False):
    """
    Multiple forecasts of a VAR using Lasso or adaptive Lasso.
    
    This function performs pseudo out-of-sample forecasting experiments.
    For true out-of-sample forecasts, use the predict method of LassoVAR.
    
    Parameters
    ----------
    dat : pd.DataFrame or array-like
        Time series data (T x N)
    exo : pd.DataFrame or array-like, optional
        Exogenous variables (T x K) (default: None)
    fc_train : int
        Number of training observations
    horizon : int, optional
        Forecast horizon (default: 1)
    lags : int, optional
        Number of lags (default: 1)
    fc_window : str, optional
        'fix' or 'expanding' window (default: 'fix')
    fc_type : str, optional
        'recursive' or 'direct' forecasting (default: 'recursive')
        Note: Only recursive currently implemented
    ic : str, optional
        Information criterion ('BIC' or 'AIC', default: 'BIC')
    adaptive : str, optional
        Adaptive Lasso type: 'none', 'ols', 'lasso', 'ridge', 'group'
        (default: 'none')
    mc : bool, optional
        Parallelize across forecasts (default: False)
    n_jobs : int, optional
        Number of parallel jobs (default: -1, all cores)
    silent : bool, optional
        Suppress output (default: False)
    trend : bool, optional
        Include linear trend (default: False)
    post : bool, optional
        Post-Lasso OLS (default: False)
        
    Returns
    -------
    dict
        Dictionary with forecast results:
        - 'err': Forecast errors (T_fc x N)
        - 'pred': Predictions (T_fc x N)
        - 'coefficients': List of coefficient matrices
        - 'lambda': Lambda values (T_fc x N)
        - 'spectest': Specification tests (n_tests x N x T_fc)
        - 'call': Call parameters
    
    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> from lassovar import forecast_lassovar
    >>> 
    >>> # Generate sample data
    >>> data = pd.DataFrame(np.random.randn(100, 3))
    >>> 
    >>> # Perform forecasting experiment
    >>> fc_results = forecast_lassovar(
    ...     data, fc_train=80, horizon=1, lags=1,
    ...     fc_window='expanding', silent=True
    ... )
    >>> 
    >>> # Access forecast errors
    >>> print(fc_results['err'])
    """
    # Validate inputs
    dat = coerce_to_dataframe(dat)
    if exo is not None:
        exo = coerce_to_dataframe(exo)
    
    if fc_train is None:
        raise ValueError("fc_train must be specified")
    
    fc_window = fc_window.lower()
    if fc_window not in ['fix', 'expanding']:
        raise ValueError("fc_window must be 'fix' or 'expanding'")
    
    fc_type = fc_type.lower()
    if fc_type not in ['recursive', 'direct']:
        raise ValueError("fc_type must be 'recursive' or 'direct'")
    
    ic = ic.upper()
    if ic not in ['BIC', 'AIC']:
        ic = 'BIC'
    
    adaptive = adaptive.lower()
    if adaptive not in ['none', 'ols', 'lasso', 'ridge', 'group']:
        adaptive = 'none'
    
    # Number of forecasts
    nbr_fc = len(dat) - fc_train
    
    # Print information
    if not silent:
        print('\n\t-----------------------------\t')
        print('Lassovar forecast')
        estimator_name = 'Adaptive Lasso' if adaptive != 'none' else 'Lasso'
        print(f'Estimator: {estimator_name}')
        if adaptive != 'none':
            print(f'Initial Estimator: {adaptive}')
        print(f'Number of equations: {dat.shape[1]}')
        if mc:
            print(f'Number of cores used: {n_jobs if n_jobs > 0 else "all"}')
        print(f'{fc_window} window forecasts')
        print(f'{horizon}-steps ahead {fc_type} forecasts. ' +
              f'Initial training sample: {fc_train} observations.')
        print(f'Number of forecasts: {nbr_fc}')
        if mc:
            print(f'Forecast level multicore enabled, #cores: {n_jobs if n_jobs > 0 else "all"}')
        print('\n\t-----------------------------\t')
    
    # For direct forecasting, adjust horizon in model estimation
    if fc_type == 'direct':
        model_horizon = horizon
    else:
        model_horizon = 1
    
    # Determine if equation-level parallelization should be used
    mc_eq = not mc  # If parallelizing forecasts, don't parallelize equations
    
    # Run forecasting loop
    n_forecasts = nbr_fc - horizon + 1
    
    if mc and n_jobs != 1:
        fc_results = Parallel(n_jobs=n_jobs)(
            delayed(forecast_loop_iteration)(
                fc, dat, fc_train, fc_window, lags, model_horizon,
                ic, exo, adaptive, post, trend, silent, mc_eq,
                n_jobs if mc_eq else 1
            ) for fc in range(n_forecasts)
        )
    else:
        fc_results = [
            forecast_loop_iteration(
                fc, dat, fc_train, fc_window, lags, model_horizon,
                ic, exo, adaptive, post, trend, silent, False, 1
            ) for fc in range(n_forecasts)
        ]
    
    # Aggregate results
    err = np.array([r['fcerr'] for r in fc_results])
    pred = np.array([r['pred'] for r in fc_results])
    lambdas = np.array([r['lambda'] for r in fc_results])
    coefficients = [r['coefficients'] for r in fc_results]
    
    # Specification tests
    n_tests = fc_results[0]['spectest'].shape[0]
    n_eq = dat.shape[1]
    spectest = np.zeros((n_tests, n_eq, len(fc_results)))
    for i, r in enumerate(fc_results):
        spectest[:, :, i] = r['spectest']
    
    # Post-Lasso OLS if requested
    post_results = None
    if post:
        post_results = [r.get('post') for r in fc_results]
    
    # Return results
    result = {
        'call': {
            'fc_train': fc_train,
            'horizon': horizon,
            'lags': lags,
            'fc_window': fc_window,
            'fc_type': fc_type,
            'ic': ic,
            'adaptive': adaptive,
            'trend': trend,
            'post': post
        },
        'err': pd.DataFrame(err, columns=dat.columns),
        'pred': pd.DataFrame(pred, columns=dat.columns),
        'coefficients': coefficients,
        'lambda': lambdas,
        'spectest': spectest
    }
    
    if post_results is not None:
        result['post'] = post_results
    
    return result
