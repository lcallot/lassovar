"""
Test LassoVAR methods.

These tests correspond to test_methods.R in the original R package.
"""

import numpy as np
import pandas as pd
import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lassovar import LassoVAR


def generate_var1_data(nobs=500):
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


def test_summary_method():
    """Test summary method."""
    simdata = generate_var1_data(nobs=500)
    exo = pd.DataFrame({'exovar': np.arange(500)})
    
    model = LassoVAR(simdata, lags=1, exo=exo)
    model.fit()
    
    # Should return a DataFrame (or None)
    result = model.summary()
    assert result is None or isinstance(result, pd.DataFrame)


def test_predict_method():
    """Test predict method."""
    simdata = generate_var1_data(nobs=500)
    exo = pd.DataFrame({'exovar': np.arange(500)})
    
    model = LassoVAR(simdata, lags=1, exo=exo)
    model.fit()
    
    # Create new data for prediction (last observation's features)
    newdata = model.x.iloc[-2:]
    
    # Get predictions
    predictions = model.predict(newdata)
    
    # Check output
    assert isinstance(predictions, (np.ndarray, pd.DataFrame))
    if isinstance(predictions, pd.DataFrame):
        assert predictions.shape[1] == 3  # 3 equations
        assert len(predictions) == 2  # 2 observations
    else:
        assert predictions.shape[1] == 3


def test_residuals_method():
    """Test residuals method."""
    simdata = generate_var1_data(nobs=500)
    exo = pd.DataFrame({'exovar': np.arange(500)})
    
    model = LassoVAR(simdata, lags=1, exo=exo)
    model.fit()
    
    # Get residuals
    residuals = model.residuals()
    
    # Check output
    assert isinstance(residuals, pd.DataFrame)
    assert residuals.shape[1] == 3  # 3 equations
    assert len(residuals) == len(model.y)  # Same length as y


def test_coef_method():
    """Test coefficient extraction."""
    simdata = generate_var1_data(nobs=500)
    
    model = LassoVAR(simdata, lags=1)
    model.fit()
    
    coef = model.coef()
    
    assert isinstance(coef, np.ndarray)
    assert coef.shape == (4, 3)  # (intercept + 3 lagged vars) x 3 equations


def test_methods_before_fit():
    """Test that methods raise error before fitting."""
    simdata = generate_var1_data(nobs=500)
    
    model = LassoVAR(simdata, lags=1)
    
    # Should raise error before fit
    with pytest.raises(RuntimeError):
        model.predict(model.x.iloc[-1:])
    
    with pytest.raises(RuntimeError):
        model.residuals()
    
    with pytest.raises(RuntimeError):
        model.summary()
    
    with pytest.raises(RuntimeError):
        model.coef()


def test_predict_with_array():
    """Test prediction with numpy array input."""
    simdata = generate_var1_data(nobs=500)
    
    model = LassoVAR(simdata, lags=1)
    model.fit()
    
    # Use numpy array for prediction
    newdata = model.x.values[-1:]
    predictions = model.predict(newdata)
    
    assert predictions is not None
    # Check shape
    if isinstance(predictions, (np.ndarray, pd.DataFrame)):
        assert predictions.shape[0] >= 1


def test_predict_single_observation():
    """Test prediction with single observation."""
    simdata = generate_var1_data(nobs=500)
    
    model = LassoVAR(simdata, lags=1)
    model.fit()
    
    # Single observation as 1D array
    newdata = model.x.values[-1]
    predictions = model.predict(newdata)
    
    assert predictions is not None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
