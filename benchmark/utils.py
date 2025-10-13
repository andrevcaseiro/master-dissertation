#!/usr/bin/env python3
import numpy as np
from scipy.optimize import curve_fit

# Global parameter values
N_o = 10000
N_values = [1000000, 2500000, 5000000, 10000000, 25000000, 50000000, 100000000, 50000000]
M_values = [100, 150, 250, 500, 750, 1000, 2500, 5000, 10000, 25000, 5000]
N_o_values = [100, 250, 500, 1000, 2500, 5000, 10000, 25000, 10000]
parameter_values = {
    'trap': [10, 100, 1000, 10000],  # Simple list of N values for trapezoidal method
    'mc': [
        (10000, 100000000, 10000),
        (10000,  50000000, 10000),
        (10000,  25000000, 10000),
        (10000,  10000000, 10000),
        #(5000,  100000000, 10000),
        (5000,   50000000, 10000),
        (5000,   25000000, 10000),
        (5000,   10000000, 10000),
        #(2000,  100000000, 10000),
        #(2000,   50000000, 10000),
        (2000,   25000000, 10000),
        (2000,   10000000, 10000),
        (1000,  100000000, 10000),
        (1000,   50000000, 10000),
        (1000,   25000000, 10000),
        (1000,   10000000, 10000),
        (500,   100000000, 10000),
        (500,    50000000, 10000),
        (500,    25000000, 10000),
        (500,    10000000, 10000),
    ]
}

""" N_o = 10000
N_o_values = [1000, 10000]
N_values = [1000000, 10000000, 25000000, 50000000]
M_values = [1000, 10000]
parameter_values = {
    'trap': [10, 100, 1000, 10000],  # Simple list of N values for trapezoidal method
    'mc': [
        (10000,   50000000, 10000),
    ]
} """
""" N_o = 10000
N_values = [10000, 20000, 50000]
M_values = [100, 250, 500, 1000]
# Print step values for analyzing the effect of print frequency on runtime and error
N_o_values = [100, 500, 1000, 5000, 10000, 20000]
parameter_values = {
    'trap': [100, 1000, 10000],  # Simple list of N values for trapezoidal method
    'mc': [
        (10000,  10000, 1000),
        (10000,  25000, 1000),
        (10000,  50000, 1000),
        #(10000, 100000, 1000),
        #(10000,  10000, 5000),
        #(10000,  20000, 5000),
        #(10000,  50000, 5000),
        #(10000, 100000, 5000),
        
    ]
} """


def calculate_errors(df, ref_df):
    """
    Calculate error metrics between simulation results and reference data.
    
    Args:
        df (pandas.DataFrame): Simulation data with time and voltage columns
        ref_df (pandas.DataFrame): Reference data with time and voltage columns
        
    Returns:
        dict: Dictionary with error metrics
    """

    # Determine which dataset is more dense
    if len(ref_df) > len(df):
        # Reference is more dense, interpolate df to match reference points
        sim_voltage = np.interp(ref_df['time'], df['time'], df['voltage'])
        ref_voltage = ref_df['voltage'].values
        
        # Use reference time points for error calculation
        time_points = ref_df['time']
    else:
        # Simulation is more dense, interpolate reference to match simulation points
        sim_voltage = df['voltage'].values
        ref_voltage = np.interp(df['time'], ref_df['time'], ref_df['voltage'])
        
        # Use simulation time points for error calculation
        time_points = df['time']
    
    abs_error = np.abs(sim_voltage - ref_voltage)
    rel_error = abs_error / (np.abs(ref_voltage) + 1e-10)  # Avoid division by zero
    
    # Calculate time intervals for weighted averages
    time_array = time_points.values if hasattr(time_points, 'values') else time_points
    dt = np.diff(time_array)
    
    # For weighted averages, we need weights for each time point
    # Use trapezoidal rule approach: each point gets half the interval on each side
    weights = np.zeros_like(time_array)
    
    if len(time_array) > 1:
        # First point gets half of the first interval
        weights[0] = dt[0] / 2
        
        # Middle points get half of the intervals on both sides
        for i in range(1, len(time_array) - 1):
            weights[i] = (dt[i-1] + dt[i]) / 2
        
        # Last point gets half of the last interval
        weights[-1] = dt[-1] / 2
        
        # Normalize weights to sum to 1 for proper averaging
        weights = weights / np.sum(weights)
    else:
        weights[0] = 1.0  # Single point case
    
    return {
        'max_error': np.max(abs_error),                           # Maximum absolute error
        'max_rel_error': np.max(rel_error),                      # Maximum relative error
        'avg_error': np.average(abs_error, weights=weights),     # Time-weighted average absolute error
        'avg_rel_error': np.average(rel_error, weights=weights), # Time-weighted average relative error
        'rms_error': np.sqrt(np.average(abs_error**2, weights=weights)),      # Time-weighted RMS error (absolute)
        'rms_rel_error': np.sqrt(np.average(rel_error**2, weights=weights)),  # Time-weighted RMS relative error
        'time_points': time_points
    }

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
    
    def get_equation(param_name="x"):
        return f"{format_scientific_notation(k)} {param_name}^{{{alpha:.2f}}}"
    return {
        'k': k,
        'alpha': alpha,
        'r_squared': r_squared,
        'residuals': residuals,
        'log_residuals': log_residuals,
        'y_fitted': y_fitted,
        'x_smooth': x_smooth,
        'y_smooth': y_smooth,
        'equation': get_equation
    }


def fit_power_law_with_offset(x, y):
    """
    Fit a power law with offset y = k * x^alpha + c using iterative non-linear fitting.
    
    Args:
        x (array-like): Independent variable
        y (array-like): Dependent variable
        
    Returns:
        dict: Dictionary containing:
            - 'k': Multiplicative constant
            - 'alpha': Power law exponent
            - 'c': Constant offset
            - 'r_squared': Coefficient of determination
            - 'residuals': Residuals in linear space
            - 'y_fitted': Fitted y values at original x points
            - 'x_smooth': Smooth x values for plotting
            - 'y_smooth': Smooth y values for plotting
            - 'equation': Human-readable equation string
    """
    def power_law_with_offset(x_val, k, alpha, c):
        """Power law with offset: y = k * x^alpha + c"""
        return k * (x_val ** alpha) + c
    
    # Convert to numpy arrays
    x = np.array(x)
    y = np.array(y)
    
    try:
        # Initial parameter guesses
        # For k: use ratio of y range to x range with some power
        y_range = np.max(y) - np.min(y)
        x_range = np.max(x) - np.min(x)
        k_guess = y_range / (x_range ** 0.5)  # Assume moderate power
        
        # For alpha: start with a reasonable power law exponent
        alpha_guess = 1
        
        # For c: start with minimum y value as offset
        c_guess = np.min(y)
        
        # Set reasonable bounds
        # k should be positive
        # alpha can be negative or positive but reasonable range
        # c can be any value but let's constrain it somewhat
        lower_bounds = [1e-10, -3.0, np.min(y) - np.abs(np.min(y))]
        upper_bounds = [np.inf, 3.0, np.max(y)]
        
        # Fit the power law with offset
        popt, pcov = curve_fit(
            power_law_with_offset, 
            x, y, 
            p0=[k_guess, alpha_guess, c_guess],
            bounds=(lower_bounds, upper_bounds),
            maxfev=10000
        )
        
        k, alpha, c = popt
        
        # Calculate fitted values
        y_fitted = power_law_with_offset(x, k, alpha, c)
        
        # Calculate R-squared
        ss_res = np.sum((y - y_fitted) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        # Calculate residuals
        residuals = y - y_fitted
        
        # Generate smooth curve for plotting
        x_smooth = np.logspace(np.log10(x.min()), np.log10(x.max()), 100)
        y_smooth = power_law_with_offset(x_smooth, k, alpha, c)
        
        # Create equation function that receives parameter name
        def get_equation(param_name="x"):
            return f"{format_scientific_notation(k)} {param_name}^{{{alpha:.2f}}} + {format_scientific_notation(c)}"
        
        return {
            'k': k,
            'alpha': alpha,
            'c': c,
            'r_squared': r_squared,
            'residuals': residuals,
            'y_fitted': y_fitted,
            'x_smooth': x_smooth,
            'y_smooth': y_smooth,
            'equation': get_equation
        }
        
    except Exception as e:
        print(f"Error fitting power law with offset: {e}")
        # Fallback to simple linear fit if power law fails
        coeffs = np.polyfit(x, y, 1)
        k_fallback = coeffs[0]
        c_fallback = coeffs[1]
        alpha_fallback = 1.0
        
        y_fitted = k_fallback * x + c_fallback
        residuals = y - y_fitted
        
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        x_smooth = np.linspace(x.min(), x.max(), 100)
        y_smooth = k_fallback * x_smooth + c_fallback
        
        # Create equation function for fallback case
        def get_equation_fallback(param_name="x"):
            return f"Linear fallback: {k_fallback:.2e} {param_name} + {c_fallback:.2e}"
        
        return {
            'k': k_fallback,
            'alpha': alpha_fallback,
            'c': c_fallback,
            'r_squared': r_squared,
            'residuals': residuals,
            'y_fitted': y_fitted,
            'x_smooth': x_smooth,
            'y_smooth': y_smooth,
            'equation': get_equation_fallback
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

def format_scientific_notation(value):
    """
    Format a number in scientific notation for LaTeX display.
    
    Args:
        value (float): Number to format
        
    Returns:
        str: LaTeX formatted scientific notation string
    """
    if value == 0:
        return "0"
    
    # Get the exponent
    exponent = int(np.floor(np.log10(abs(value))))
    
    # Get the mantissa
    mantissa = value / (10 ** exponent)
    
    # Format based on exponent
    if exponent == 0:
        return f"{mantissa:.1f}"
    elif exponent == 1:
        return f"{value:.0f}"
    else:
        return f"{mantissa:.1f} \\times 10^{{{exponent}}}"