#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import multiprocessing
from tran import run_monte_carlo, run_trapezoidal, run_ngspice


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
    
    return {
        'max_error': np.max(rel_error),
        'avg_error': np.mean(rel_error),
        'rms_error': np.sqrt(np.mean(rel_error**2)),
        'time_points': time_points,
        'abs_error': abs_error,
        'rel_error': rel_error
    }

def plot_error_vs_parameter(results_df, parameter_name, output_dir, fixed_param_info=""):
    """
    Plot error metrics vs a parameter (N or M) with power law trendlines.
    
    Args:
        results_df (pandas.DataFrame): Results dataframe with error metrics
        parameter_name (str): Name of the parameter ('N' or 'M')
        output_dir (Path): Directory to save plots
        fixed_param_info (str): Information about fixed parameter for title
    """
    plt.figure(figsize=(5, 4))  # Smaller, more document-friendly size
    
    x = results_df[parameter_name].values
    
    # Plot data points
    plt.plot(x, results_df['max_error'], marker='o', label='Max Error', linewidth=1.5)
    # plt.plot(x, results_df['avg_error'], marker='x', label='Avg Error', linewidth=1.5)
    plt.plot(x, results_df['rms_error'], marker='s', label='RMS Error', linewidth=1.5)
    
    # Fit and plot trendlines for each error type
    colors = ['C0', 'C1']  # Default matplotlib colors
    error_types = ['max_error', 'rms_error'] # 'Avg Error' occulted
    error_labels = ['Max Error', 'RMS Error']
    
    for i, (error_type, error_label, color) in enumerate(zip(error_types, error_labels, colors)):
        y = results_df[error_type].values
        fit_result = fit_power_law(x, y)
        
        # Plot trendline using pre-calculated smooth values
        plt.plot(fit_result['x_smooth'], fit_result['y_smooth'], '--', color=color, alpha=0.7, linewidth=2,
                label=f'{error_label} Trendline: {fit_result["k"]:.2e} × {parameter_name}^{fit_result["alpha"]:.2f} (R² = {fit_result["r_squared"]:.3f})')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(f'{parameter_name} (log₁₀ scale)')
    plt.ylabel('Error (log₁₀ scale)')
    # plt.title(f'Error vs {parameter_name}{fixed_param_info}')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_dir / f'error_vs_{parameter_name}.pdf', bbox_inches='tight')
    plt.close()

def plot_exec_time_vs_parameter(results_df, parameter_name, output_dir, fixed_param_info="", trapezoidal_time=None):
    """
    Plot execution time vs a parameter (N or M) with power law trendline.
    
    Args:
        results_df (pandas.DataFrame): Results dataframe with execution times
        parameter_name (str): Name of the parameter ('N' or 'M')
        output_dir (Path): Directory to save plots
        fixed_param_info (str): Information about fixed parameter for title
        trapezoidal_time (float, optional): Trapezoidal execution time for reference line
    """
    plt.figure(figsize=(4, 3))  # Compact size for execution time plots
    
    x = results_df[parameter_name].values
    y = results_df['exec_time'].values
    
    # Plot data points
    plt.plot(x, y, marker='o', label='Execution Time', linewidth=1.5, markersize=6)
    
    # Fit and plot trendline
    fit_result = fit_power_law(x, y)
    
    # Plot trendline using pre-calculated smooth values
    plt.plot(fit_result['x_smooth'], fit_result['y_smooth'], '--', color='red', alpha=0.8, linewidth=2,
            label=f'Time Trendline: {fit_result["k"]:.2e} × {parameter_name}^{fit_result["alpha"]:.2f} (R² = {fit_result["r_squared"]:.3f})')
    
    # Add horizontal reference line for trapezoidal execution time
    if trapezoidal_time is not None:
        plt.axhline(y=trapezoidal_time, color='green', linestyle='-', linewidth=2, alpha=0.7,
                   label=f'Trapezoidal Reference: {trapezoidal_time:.3f}s')
    
    plt.xscale('log')
    plt.xlabel(f'{parameter_name} (log₁₀ scale)')
    plt.ylabel('Execution Time (s)')
    # plt.title(f'Execution Time vs {parameter_name}{fixed_param_info}')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_dir / f'exec_time_vs_{parameter_name}.pdf')
    plt.close()

def plot_speedup_analysis(speedup_df, output_dir):
    """
    Plot speedup analysis vs number of threads.
    
    Args:
        speedup_df (pandas.DataFrame): Speedup results dataframe
        output_dir (Path): Directory to save plots
    """
    plt.figure()
    plt.plot(speedup_df['threads'], speedup_df['speedup'], marker='o', label='Measured Speedup')
    plt.plot(speedup_df['threads'], speedup_df['threads'], 'k--', label='Ideal Linear Speedup')
    plt.xscale('log', base=2)
    plt.yscale('log', base=2)
    plt.xlabel('Threads (log₂ scale)')
    plt.ylabel('Speedup (log₂ scale)')
    # plt.title('Scalability: Speedup vs Threads')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_dir / 'speedup_vs_threads.pdf', bbox_inches='tight')
    plt.close()

def plot_voltage_comparison(plot_data, output_dir, filename='voltage_comparison.pdf', title='Voltage Comparison'):
    """
    Plot comparison of voltage data from multiple sources.
    
    Args:
        plot_data (list): List of dictionaries, each containing:
            - 'df': DataFrame with time and voltage columns
            - 'label': Label for the plot legend
            - 'style': Dictionary with matplotlib plot style parameters (e.g., {'color': 'r', 'linestyle': '--', 'linewidth': 2})
        output_dir (Path): Directory to save plots
        filename (str): Filename for the saved plot
        title (str): Title for the plot
    """
    plt.figure(figsize=(8, 5))  # Good size for voltage comparison plots
    
    # Plot all data
    for data in plot_data:
        df = data['df']
        label = data['label']
        style = data.get('style', {})  # Default to empty dict if no style specified
        
        # Set default style parameters if not provided
        default_style = {'linewidth': 1.5, 'alpha': 0.8}
        plot_style = {**default_style, **style}
        
        plt.plot(df['time'], df['voltage'], label=label, **plot_style)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    # plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / filename, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Run parameter sweeps for ODE solvers and plot results')
    parser.add_argument('netlist_path', help='Path to the SPICE netlist file')
    parser.add_argument('--output-dir', type=str, help='Directory to save plots (default: benchmark_dir/\{netlist_name\})')
    parser.add_argument('--include-spice', action='store_true', help='Include SPICE (NGSpice) simulation in comparison plots')

    args = parser.parse_args()

    # Set default output directory if not provided
    if args.output_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        netlist_name = os.path.splitext(os.path.basename(args.netlist_path))[0]
        args.output_dir = os.path.join(script_dir, netlist_name)

    # Configure plot style
    plt.rcParams.update({"legend.fontsize": 6, "savefig.bbox": "tight"})

    # Calculate reference df based on highest values of N and M
    N_values = [1000, 5000, 10000, 40000, 80000, 160000, 320000]
    M_values = [100, 1000, 10000]
    min_N, max_N = min(N_values), max(N_values)
    min_M, max_M = min(M_values), max(M_values)

    # Run the simulation with the highest N and M to get the reference data
    # Assuming run_trapezoidal returns a DataFrame with 'time' and 'voltage' columns
    ref_df, trapezoidal_exec_time, _ = run_trapezoidal(args.netlist_path, 0, 1000000)
    print(f"Trapezoidal reference simulation: exec_time={trapezoidal_exec_time:.4f}s")

    # Run SPICE simulation if requested
    spice_df = None
    if args.include_spice:
        print("Running SPICE simulation...")
        spice_df, _, _ = run_ngspice(args.netlist_path)
        print(f"SPICE simulation completed. Got {len(spice_df)} data points.")
        
        # Calculate error between trapezoidal and SPICE (treating SPICE as ground truth)
        if spice_df is not None and len(spice_df) > 0:
            trap_vs_spice_errors = calculate_errors(ref_df, spice_df)
            print(f"Trapezoidal vs SPICE error metrics:")
            print(f"  Max error: {trap_vs_spice_errors['max_error']:.4e}")
            print(f"  Avg error: {trap_vs_spice_errors['avg_error']:.4e}")
            print(f"  RMS error: {trap_vs_spice_errors['rms_error']:.4e}")
            
            # Create output directory and save error metrics to file
            output_dir = Path(args.output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Save trapezoidal vs SPICE error metrics
            error_metrics = pd.DataFrame([{
                'comparison': 'Trapezoidal_vs_SPICE',
                'max_error': trap_vs_spice_errors['max_error'],
                'avg_error': trap_vs_spice_errors['avg_error'],
                'rms_error': trap_vs_spice_errors['rms_error'],
                'trapezoidal_exec_time': trapezoidal_exec_time
            }])
            error_metrics.to_csv(output_dir / 'trapezoidal_vs_spice_errors.csv', index=False)
            print(f"Error metrics saved to: {output_dir / 'trapezoidal_vs_spice_errors.csv'}")

    # Store all Monte Carlo results in a dictionary with (N,M) as key
    mc_results = {}

    results = []
    for N in N_values:
        # Run simulation for current N and max_M
        sim_df, exec_time, _ = run_monte_carlo(args.netlist_path, 0, N, max_M, print_step=1000)
        mc_results[(N, max_M)] = sim_df  # Store in dictionary

        errors = calculate_errors(sim_df, ref_df)
        results.append({
            'N': N,
            'M': max_M,
            'exec_time': exec_time,
            'max_error': errors['max_error'],
            'avg_error': errors['avg_error'],
            'rms_error': errors['rms_error']
        })
        print(f"N={N}, exec_time={exec_time:.4f}s, max_error={errors['max_error']:.4e}")

    # Store results as DataFrame
    results_df = pd.DataFrame(results)
    # Output directory already created earlier (if SPICE was run) or create it now
    if not Path(args.output_dir).exists():
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = Path(args.output_dir)
    results_df.to_csv(output_dir / 'sweep_N_results.csv', index=False)

    # Plot results for N sweep
    plot_error_vs_parameter(results_df, 'N', output_dir, " (M fixed at max)")
    plot_exec_time_vs_parameter(results_df, 'N', output_dir, " (M fixed at max)", trapezoidal_exec_time)

    results = []
    for M in M_values:
        # Run simulation for current M and max_N
        sim_df, exec_time, _ = run_monte_carlo(args.netlist_path, 0, max_N, M, print_step=1000)
        mc_results[(max_N, M)] = sim_df  # Store in dictionary

        errors = calculate_errors(sim_df, ref_df)
        results.append({
            'N': max_N,
            'M': M,
            'exec_time': exec_time,
            'max_error': errors['max_error'],
            'avg_error': errors['avg_error'],
            'rms_error': errors['rms_error']
        })
        print(f"M={M}, exec_time={exec_time:.4f}s, max_error={errors['max_error']:.4e}")

    # Store results as DataFrame
    results_df = pd.DataFrame(results)
    # Output directory should already exist from earlier steps
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)  # Safe to call multiple times
    results_df.to_csv(output_dir / 'sweep_M_results.csv', index=False)

    # Plot results for M sweep
    plot_error_vs_parameter(results_df, 'M', output_dir, " (N fixed at max)")
    plot_exec_time_vs_parameter(results_df, 'M', output_dir, " (N fixed at max)", trapezoidal_exec_time)

    # Prepare plot data
    plot_data = [
        {
            'df': ref_df,
            'label': 'Reference (Trapezoidal)',
            'style': {'color': 'k', 'linestyle': '-', 'linewidth': 2}
        },
        {
            'df': mc_results[(max_N, max_M)],
            'label': f'Best MC (N={max_N}, M={max_M})',
            'style': {'color': 'g', 'linestyle': '-', 'linewidth': 1.5}
        },
        {
            'df': mc_results[(max_N, min_M)],
            'label': f'Low M (N={max_N}, M={min_M})',
            'style': {'color': 'r', 'linestyle': '--', 'linewidth': 1.5}
        },
        {
            'df': mc_results[(min_N, max_M)],
            'label': f'Low N (N={min_N}, M={max_M})',
            'style': {'color': 'b', 'linestyle': ':', 'linewidth': 1.5}
        }
    ]

    # Add SPICE data if available
    if spice_df is not None:
        plot_data.append({
            'df': spice_df,
            'label': 'SPICE (NGSpice)',
            'style': {'color': 'm', 'linestyle': '-', 'linewidth': 2, 'alpha': 0.7}
        })

    # Plot solution comparison showing actual data
    plot_voltage_comparison(
        plot_data=plot_data,
        output_dir=output_dir,
        filename='solution_comparison.pdf',
        title='Solution Comparison: Effect of N and M Parameters'
    )

    # Scalability (speedup) analysis: vary number of threads from 1 to available CPU cores
    max_threads = multiprocessing.cpu_count()
    thread_values = [2 ** i for i in range(int(np.log2(max_threads)) + 1) if 2 ** i <= max_threads]
    if 1 not in thread_values:
        thread_values = [1] + thread_values  # Ensure 1 is included

    # Use max_N and max_M for a fair comparison
    baseline_time = None
    speedup_results = []
    for threads in thread_values:
        # Assume run_monte_carlo accepts a 'threads' argument
        sim_df, exec_time = run_monte_carlo(args.netlist_path, 0, max_N, max_M, print_step=1000, num_threads=threads)
        if not baseline_time:
            baseline_time = exec_time
        speedup = baseline_time / exec_time
        speedup_results.append({'threads': threads, 'exec_time': exec_time, 'speedup': speedup})
        print(f"Threads={threads}, exec_time={exec_time:.4f}s, speedup={speedup:.2f}")

    speedup_df = pd.DataFrame(speedup_results)
    speedup_df.to_csv(output_dir / 'scalability_results.csv', index=False)

    # Plot speedup analysis
    plot_speedup_analysis(speedup_df, output_dir)

if __name__ == "__main__":
    main()
