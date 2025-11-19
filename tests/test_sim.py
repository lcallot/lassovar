"""
Test basic LassoVAR functionality.

These tests correspond to test_sim.R in the original R package.
"""

import numpy as np
import pandas as pd
import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lassovar import LassoVAR


def generate_var1_data(nobs=20):
    """
    Generate simulated VAR(1) data.
    
    A simple VAR(1) parameter matrix with:
    - Diagonal elements = 0.4
    - A[2,1] = A[1,3] = -0.4
    """
    A = np.zeros((3, 3))
    np.fill_diagonal(A, 0.4)
    A[1, 0] = -0.4  # A[2,1] in 1-indexed
    A[0, 2] = -0.4  # A[1,3] in 1-indexed
    
    # Simulate data
    np.random.seed(42)
    simdata = np.zeros((nobs + 1, 3))
    simdata[0, :] = np.random.randn(3)
    
    for t in range(nobs):
        simdata[t + 1, :] = simdata[t, :] @ A + np.random.randn(3)
    
    simdata = simdata[1:]  # Remove first row (initialization)
    
    # Create DataFrame with column names
    df = pd.DataFrame(simdata, columns=['Cyrus', 'Cambyses', 'Darius'])
    return df


def test_basic_lassovar():
    """Test basic Lasso VAR estimation."""
    simdata = generate_var1_data(nobs=20)
    
    model = LassoVAR(simdata, lags=1)
    model.fit()
    
    assert isinstance(model, LassoVAR)
    assert model.fitted_
    assert model.coefficients is not None
    assert model.coefficients.shape == (2, 3)  # 1 lag * 3 vars + intercept


def test_lassovar_with_trend():
    """Test Lasso VAR with trend."""
    simdata = generate_var1_data(nobs=20)
    
    model = LassoVAR(simdata, lags=1, trend=True)
    model.fit()
    
    assert isinstance(model, LassoVAR)
    assert model.fitted_
    assert model.trend


def test_lassovar_with_aic():
    """Test Lasso VAR with AIC criterion."""
    simdata = generate_var1_data(nobs=20)
    
    model = LassoVAR(simdata, lags=1, ic='AIC')
    model.fit()
    
    assert isinstance(model, LassoVAR)
    assert model.fitted_
    assert model.ic == 'AIC'


def test_lassovar_multiple_lags():
    """Test Lasso VAR with multiple lags."""
    simdata = generate_var1_data(nobs=20)
    
    model = LassoVAR(simdata, lags=2)
    model.fit()
    
    assert isinstance(model, LassoVAR)
    assert model.fitted_
    assert model.lags == 2
    assert model.coefficients.shape[0] == 7  # 2 lags * 3 vars + intercept


def test_lassovar_with_exo():
    """Test Lasso VAR with exogenous variables."""
    nobs = 20
    simdata = generate_var1_data(nobs=nobs)
    exo = pd.DataFrame(np.arange(nobs), columns=['exo_var'])
    
    model = LassoVAR(simdata, lags=1, exo=exo)
    model.fit()
    
    assert isinstance(model, LassoVAR)
    assert model.fitted_
    # intercept + 3 lags + 1 exo
    assert model.coefficients.shape[0] == 5


def test_adaptive_lassovar_ols():
    """Test adaptive Lasso VAR with OLS initial estimator."""
    simdata = generate_var1_data(nobs=20)
    
    model = LassoVAR(simdata, lags=1, adaptive='ols')
    model.fit()
    
    assert isinstance(model, LassoVAR)
    assert model.fitted_
    assert model.ada_w is not None
    assert model.ada_w['ada'] == 'ols'


def test_lassovar_methods():
    """Test that model has required methods."""
    simdata = generate_var1_data(nobs=20)
    
    model = LassoVAR(simdata, lags=1)
    model.fit()
    
    # Test methods exist
    assert hasattr(model, 'predict')
    assert hasattr(model, 'residuals')
    assert hasattr(model, 'summary')
    assert hasattr(model, 'coef')
    
    # Test they can be called
    res = model.residuals()
    assert isinstance(res, pd.DataFrame)
    
    coef = model.coef()
    assert isinstance(coef, np.ndarray)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
