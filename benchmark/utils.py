#!/usr/bin/env python3
import numpy as np
from scipy.optimize import curve_fit


def fit_power_law(x, y):
    """
    Fit a power law y = k * x^alpha using log-log linear regression.
    
    Args:
        x (array-like): Independent variable
        y (array-like): Dependent variable
        
    Returns:
        dict: Dictionary containing:
            - 'k': Multiplicative constant
            - 'alpha': Power law exponent
            - 'r_squared': Coefficient of determination
            - 'residuals': Residuals in linear space
            - 'log_residuals': Residuals in log space
            - 'y_fitted': Fitted y values at original x points
            - 'x_smooth': Smooth x values for plotting
            - 'y_smooth': Smooth y values for plotting
    """
    # Convert to numpy arrays and take logarithms
    log_x = np.log(x)
    log_y = np.log(y)
    
    # Perform linear regression in log space: log(y) = log(k) + alpha * log(x)
    A = np.vstack([log_x, np.ones(len(log_x))]).T
    coeffs, residuals_sum, _, _ = np.linalg.lstsq(A, log_y, rcond=None)
    
    alpha = coeffs[0]
    log_k = coeffs[1]
    k = np.exp(log_k)
    
    # Calculate R-squared
    y_pred_log = log_k + alpha * log_x
    ss_res = np.sum((log_y - y_pred_log) ** 2)
    ss_tot = np.sum((log_y - np.mean(log_y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    
    # Calculate fitted values at original x points
    y_fitted = k * (x ** alpha)
    
    # Calculate residuals in both log and linear space
    residuals = y - y_fitted
    log_residuals = log_y - y_pred_log
    
    # Generate smooth curve for plotting
    x_smooth = np.logspace(np.log10(x.min()), np.log10(x.max()), 100)
    y_smooth = k * (x_smooth ** alpha)
    
    return {
        'k': k,
        'alpha': alpha,
        'r_squared': r_squared,
        'residuals': residuals,
        'log_residuals': log_residuals,
        'y_fitted': y_fitted,
        'x_smooth': x_smooth,
        'y_smooth': y_smooth
    }


def fit_amdahls_law(threads, speedup):
    """
    Fit Amdahl's Law to speedup data.
    
    Amdahl's Law: Speedup(n) = 1 / (s + (1-s)/n)
    where:
    - n is the number of threads
    - s is the serial fraction (0 ≤ s ≤ 1)
    - (1-s) is the parallel fraction
    
    Args:
        threads (array-like): Number of threads
        speedup (array-like): Measured speedup values
        
    Returns:
        dict: Dictionary containing:
            - 's': Serial fraction (0-1)
            - 'parallel_fraction': Parallel fraction (1-s)
            - 'max_speedup': Theoretical maximum speedup (1/s)
            - 'r_squared': Coefficient of determination
            - 'threads_smooth': Smooth thread values for plotting
            - 'speedup_smooth': Smooth speedup values for plotting
            - 'equation': Human-readable equation string
    """
    def amdahls_law(n, s):
        """Amdahl's Law formula"""
        return 1.0 / (s + (1.0 - s) / n)
    
    try:
        # Fit Amdahl's Law
        # Initial guess: assume 20% serial portion
        popt, pcov = curve_fit(amdahls_law, threads, speedup, 
                              p0=[0.2], 
                              bounds=(0.001, 0.999),  # s must be between 0.1% and 99.9%
                              maxfev=10000)
        
        s = popt[0]
        parallel_fraction = 1 - s
        
        # Calculate R-squared
        speedup_pred = amdahls_law(threads, s)
        ss_res = np.sum((speedup - speedup_pred) ** 2)
        ss_tot = np.sum((speedup - np.mean(speedup)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        # Generate smooth curve for plotting
        threads_smooth = np.logspace(np.log10(threads.min()), np.log10(threads.max()), 100)
        speedup_smooth = amdahls_law(threads_smooth, s)
        
        # Create human-readable equation
        equation = f"1/(s + (1-s)/n), s={s:.3f}"
        
        return {
            's': s,
            'parallel_fraction': parallel_fraction,
            'r_squared': r_squared,
            'threads_smooth': threads_smooth,
            'speedup_smooth': speedup_smooth,
            'equation': equation,
            'fitted_speedup': speedup_pred
        }
        
    except Exception as e:
        print(f"Error fitting Amdahl's Law: {e}")
        # Return fallback values
        return {
            's': 0.5,
            'parallel_fraction': 0.5,
            'max_speedup': 2.0,
            'r_squared': 0.0,
            'threads_smooth': threads,
            'speedup_smooth': np.ones_like(threads),
            'equation': "Fit failed",
            'fitted_speedup': np.ones_like(speedup)
        }