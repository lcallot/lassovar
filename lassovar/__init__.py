"""
lassovar: Estimation and forecasting of VAR models with the Lasso

This package provides tools for estimating Vector Autoregression models
using Lasso and adaptive Lasso estimators, with support for forecasting
and parallel computation.

Main Components
---------------
LassoVAR : class
    Main class for fitting VAR models with Lasso/adaptive Lasso
    
forecast_lassovar : function
    Perform pseudo out-of-sample forecasting experiments

Examples
--------
Basic usage:

>>> import numpy as np
>>> import pandas as pd
>>> from lassovar import LassoVAR, forecast_lassovar
>>>
>>> # Create sample data
>>> data = pd.DataFrame(np.random.randn(100, 3), columns=['V1', 'V2', 'V3'])
>>>
>>> # Fit a Lasso VAR model
>>> model = LassoVAR(data, lags=1)
>>> model.fit()
>>> model.summary()
>>>
>>> # Perform forecasting
>>> fc_results = forecast_lassovar(data, fc_train=80, horizon=1, lags=1)
"""

__version__ = '0.9.0'
__author__ = 'Laurent Callot (Original R package), Python port'
__license__ = 'MIT'

from .lassovar import LassoVAR
from .forecasting import forecast_lassovar
from .utils import make_var_data, specification_tests
from .adaptive import get_adaptive_weights

__all__ = [
    'LassoVAR',
    'forecast_lassovar',
    'make_var_data',
    'specification_tests',
    'get_adaptive_weights',
]
