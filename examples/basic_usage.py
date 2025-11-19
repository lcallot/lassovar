"""
Basic usage examples for lassovar package.
"""

import numpy as np
import pandas as pd
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lassovar import LassoVAR, forecast_lassovar


def example1_basic_lasso():
    """Example 1: Basic Lasso VAR estimation."""
    print("\n" + "="*60)
    print("Example 1: Basic Lasso VAR Estimation")
    print("="*60)
    
    # Generate sample data
    np.random.seed(42)
    data = pd.DataFrame(
        np.random.randn(100, 5),
        columns=['V1', 'V2', 'V3', 'V4', 'V5']
    )
    
    # Fit a VAR(1) model with Lasso
    model = LassoVAR(data, lags=1)
    model.fit()
    
    # Print summary
    print("\nModel Summary:")
    model.summary()
    
    # Get predictions
    last_obs = model.x.iloc[-1:]
    predictions = model.predict(last_obs)
    print("\nPredictions for next period:")
    print(predictions)
    
    # Get residuals
    residuals = model.residuals()
    print(f"\nResiduals shape: {residuals.shape}")
    print(f"Mean residuals:\n{residuals.mean()}")


def example2_adaptive_lasso():
    """Example 2: Adaptive Lasso with OLS initial estimator."""
    print("\n" + "="*60)
    print("Example 2: Adaptive Lasso VAR")
    print("="*60)
    
    # Generate sample data
    np.random.seed(42)
    data = pd.DataFrame(
        np.random.randn(100, 3),
        columns=['Series1', 'Series2', 'Series3']
    )
    
    # Fit adaptive Lasso with OLS initial estimator
    model = LassoVAR(data, lags=2, adaptive='ols', ic='BIC')
    model.fit()
    
    print("\nModel Summary:")
    model.summary(short=True)
    
    print(f"\nEstimator: {model.estimator}")
    print(f"Selected lambda values: {model.lambda_}")


def example3_with_exogenous():
    """Example 3: VAR with exogenous variables."""
    print("\n" + "="*60)
    print("Example 3: VAR with Exogenous Variables")
    print("="*60)
    
    # Generate sample data
    np.random.seed(42)
    nobs = 100
    data = pd.DataFrame(
        np.random.randn(nobs, 3),
        columns=['Y1', 'Y2', 'Y3']
    )
    
    # Generate exogenous variables
    exo = pd.DataFrame({
        'Exo1': np.sin(np.linspace(0, 4*np.pi, nobs)),
        'Exo2': np.arange(nobs)
    })
    
    # Fit model with exogenous variables
    model = LassoVAR(data, lags=1, exo=exo, trend=True)
    model.fit()
    
    print("\nModel Summary:")
    model.summary()


def example4_forecasting():
    """Example 4: Pseudo out-of-sample forecasting."""
    print("\n" + "="*60)
    print("Example 4: Pseudo Out-of-Sample Forecasting")
    print("="*60)
    
    # Generate sample data
    np.random.seed(42)
    nobs = 100
    data = pd.DataFrame(
        np.random.randn(nobs, 3),
        columns=['A', 'B', 'C']
    )
    
    # Perform forecasting experiment
    print("\nRunning forecasting experiment...")
    fc_results = forecast_lassovar(
        data,
        fc_train=80,
        horizon=1,
        lags=1,
        fc_window='expanding',
        fc_type='recursive',
        ic='BIC',
        silent=True
    )
    
    # Display results
    print(f"\nNumber of forecasts: {len(fc_results['err'])}")
    print(f"\nForecast errors (first 5):")
    print(fc_results['err'].head())
    
    print(f"\nMean Absolute Error by series:")
    print(fc_results['err'].abs().mean())
    
    print(f"\nRoot Mean Squared Error by series:")
    print(np.sqrt((fc_results['err']**2).mean()))


def example5_var_simulation():
    """Example 5: Simulate and estimate a VAR(1) process."""
    print("\n" + "="*60)
    print("Example 5: Simulate and Estimate VAR(1)")
    print("="*60)
    
    # Define VAR(1) coefficient matrix
    A = np.array([
        [0.5, 0.1, 0.0],
        [0.2, 0.4, 0.1],
        [0.0, 0.2, 0.3]
    ])
    
    print("\nTrue VAR(1) coefficient matrix:")
    print(A)
    
    # Simulate data
    np.random.seed(42)
    nobs = 200
    data = np.zeros((nobs, 3))
    data[0] = np.random.randn(3)
    
    for t in range(1, nobs):
        data[t] = data[t-1] @ A + np.random.randn(3) * 0.5
    
    df = pd.DataFrame(data, columns=['X1', 'X2', 'X3'])
    
    # Estimate with Lasso VAR
    model = LassoVAR(df, lags=1)
    model.fit()
    
    print("\nEstimated coefficients (excluding intercept):")
    print(model.coefficients[1:, :])
    
    print("\nModel Summary:")
    model.summary(short=True)


if __name__ == '__main__':
    print("\nLassoVAR Package - Usage Examples")
    print("="*60)
    
    # Run all examples
    example1_basic_lasso()
    example2_adaptive_lasso()
    example3_with_exogenous()
    example4_forecasting()
    example5_var_simulation()
    
    print("\n" + "="*60)
    print("All examples completed!")
    print("="*60)
