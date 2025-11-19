"""
Test forecasting functionality.

These tests correspond to test_forecast.R in the original R package.
"""

import numpy as np
import pandas as pd
import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lassovar import forecast_lassovar


def generate_var1_data(nobs=50):
    """Generate simulated VAR(1) data."""
    A = np.zeros((3, 3))
    np.fill_diagonal(A, 0.4)
    A[1, 0] = -0.4
    A[0, 2] = -0.4
    
    np.random.seed(42)
    simdata = np.zeros((nobs + 1, 3))
    simdata[0, :] = np.random.randn(3)
    
    for t in range(nobs):
        simdata[t + 1, :] = simdata[t, :] @ A + np.random.randn(3)
    
    simdata = simdata[1:]
    df = pd.DataFrame(simdata, columns=['Cyrus', 'Cambyses', 'Darius'])
    return df


def test_basic_forecast():
    """Test basic forecasting."""
    simdata = generate_var1_data(nobs=50)
    
    result = forecast_lassovar(
        dat=simdata,
        fc_train=48,
        horizon=1,
        fc_window='fix',
        fc_type='recursive',
        silent=True
    )
    
    assert isinstance(result, dict)
    assert 'err' in result
    assert 'pred' in result
    assert 'coefficients' in result
    assert 'lambda' in result


def test_forecast_with_exo():
    """Test forecasting with exogenous variables."""
    nobs = 50
    simdata = generate_var1_data(nobs=nobs)
    exo = pd.DataFrame(np.random.randn(nobs, 1), columns=['exo'])
    
    result = forecast_lassovar(
        dat=simdata,
        exo=exo,
        fc_train=48,
        horizon=1,
        fc_window='fix',
        fc_type='recursive',
        silent=True
    )
    
    assert isinstance(result, dict)
    assert 'err' in result


def test_forecast_expanding_window():
    """Test forecasting with expanding window."""
    simdata = generate_var1_data(nobs=50)
    
    result = forecast_lassovar(
        dat=simdata,
        fc_train=48,
        horizon=1,
        fc_window='expanding',
        fc_type='recursive',
        silent=True
    )
    
    assert isinstance(result, dict)
    assert 'err' in result
    assert result['call']['fc_window'] == 'expanding'


def test_forecast_horizon_2():
    """Test 2-step ahead forecasting."""
    simdata = generate_var1_data(nobs=50)
    
    result = forecast_lassovar(
        dat=simdata,
        fc_train=47,
        horizon=2,
        fc_window='fix',
        fc_type='recursive',
        silent=True
    )
    
    assert isinstance(result, dict)
    assert 'err' in result
    assert result['call']['horizon'] == 2


def test_forecast_direct():
    """Test direct forecasting."""
    simdata = generate_var1_data(nobs=50)
    
    result = forecast_lassovar(
        dat=simdata,
        fc_train=47,
        horizon=2,
        fc_window='fix',
        fc_type='direct',
        silent=True
    )
    
    assert isinstance(result, dict)
    assert 'err' in result
    assert result['call']['fc_type'] == 'direct'


def test_forecast_output_shapes():
    """Test that forecast outputs have correct shapes."""
    simdata = generate_var1_data(nobs=50)
    
    fc_train = 48
    horizon = 1
    
    result = forecast_lassovar(
        dat=simdata,
        fc_train=fc_train,
        horizon=horizon,
        fc_window='fix',
        fc_type='recursive',
        silent=True
    )
    
    n_forecasts = len(simdata) - fc_train - horizon + 1
    
    # Check shapes
    assert result['err'].shape[0] == n_forecasts
    assert result['err'].shape[1] == 3  # 3 equations
    assert result['pred'].shape[0] == n_forecasts
    assert result['pred'].shape[1] == 3


def test_forecast_missing_fc_train():
    """Test that error is raised when fc_train is missing."""
    simdata = generate_var1_data(nobs=50)
    
    with pytest.raises(ValueError):
        forecast_lassovar(dat=simdata, silent=True)


def test_forecast_with_adaptive():
    """Test forecasting with adaptive Lasso."""
    simdata = generate_var1_data(nobs=50)
    
    result = forecast_lassovar(
        dat=simdata,
        fc_train=48,
        horizon=1,
        adaptive='ols',
        silent=True
    )
    
    assert isinstance(result, dict)
    assert result['call']['adaptive'] == 'ols'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
