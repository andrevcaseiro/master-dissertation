#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from tran import run_monte_carlo, run_trapezoidal, run_ngspice
from scalability import run_scalability_analysis
from utils import N_o, N_values, M_values, N_o_values, fit_power_law, fit_power_law_with_offset, parameter_values, calculate_errors, format_scientific_notation
from time_vs_error_analysis import plot_time_vs_error_analysis


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
    plt.plot(x, results_df['max_error'], "C1o", label='Max Error', linewidth=1.5)
    # plt.plot(x, results_df['avg_error'], marker='x', label='Avg Error', linewidth=1.5)
    plt.plot(x, results_df['rms_error'], "C0s", label='RMS Error', linewidth=1.5)
    
    # Fit and plot trendlines for each error type
    colors = ['C1', 'C0']  # Default matplotlib colors
    error_types = ["max_error", 'rms_error'] 
    error_labels = ["Max error", 'RMS Error']
    
    for i, (error_type, error_label, color) in enumerate(zip(error_types, error_labels, colors)):
        y = results_df[error_type].values
        #fit_result = fit_power_law_with_offset(x, y)
        fit_result = fit_power_law(x, y)
        
        # Plot trendline using pre-calculated smooth values
        plt.plot(fit_result['x_smooth'], fit_result['y_smooth'], '--', color=color, alpha=0.7, linewidth=2,
                label=f'{error_label} Trendline: ${fit_result["equation"](parameter_name)}\;(R^2 = {fit_result["r_squared"]:.3f})$')
    
    plt.xscale('log')
    plt.yscale('log')
    # plt.yscale('log')
    # Create a more descriptive xlabel for print_step
    xlabel = f'${parameter_name}$ ($\log_{{10}}$ scale)'
    plt.xlabel(xlabel)
    plt.ylabel('Error (V)')
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
    plt.figure(figsize=(5, 4))  # Compact size for execution time plots
    
    x = results_df[parameter_name].values
    y = results_df['exec_time'].values
    
    # Plot data points
    plt.plot(x, y, "C0o", label=f'Execution Time{fixed_param_info}', markersize=6)
    
    # Fit and plot trendline
    #fit_result = fit_power_law_with_offset(x, y)
    fit_result = fit_power_law(x, y)
    
    # Plot trendline using pre-calculated smooth values
    plt.plot(fit_result['x_smooth'], fit_result['y_smooth'], '--', color='C0', alpha=0.8, linewidth=2,
            label=f'Time Trendline: ${fit_result["equation"](parameter_name)}\;(R^2 = {fit_result["r_squared"]:.3f})$')
    
    # Add horizontal reference line for trapezoidal execution time
    if trapezoidal_time is not None:
        plt.axhline(y=trapezoidal_time, color='green', linestyle='-', linewidth=2, alpha=0.7,
                   label=f'Trapezoidal Reference: {trapezoidal_time:.3f}s')
    
    plt.xscale('log')
    plt.yscale('log')
    # Create a more descriptive xlabel for print_step
    xlabel = f'${parameter_name}$ (log₁₀ scale)'
    plt.xlabel(xlabel)
    plt.ylabel('Execution Time (s)')
    # plt.title(f'Execution Time vs {parameter_name}{fixed_param_info}')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_dir / f'exec_time_vs_{parameter_name}.pdf')
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
    min_N, max_N = min(N_values), N_values[-1]
    min_M, max_M = min(M_values), M_values[-1]
    min_N_o, max_N_o = min(N_o_values), N_o_values[-1]

    # Run the simulation with the highest N and M to get the reference data
    # Assuming run_trapezoidal returns a DataFrame with 'time' and 'voltage' columns
    ref_df, trapezoidal_exec_time, _ = run_trapezoidal(args.netlist_path, 0, N_o)
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
        sim_df, exec_time, _ = run_monte_carlo(args.netlist_path, 0, N, max_M, print_step=max_N_o)
        mc_results[(N, max_M, max_N_o)] = sim_df  # Store in dictionary with consistent 3-tuple format

        errors = calculate_errors(sim_df, ref_df)
        results.append({
            'N': N,
            'M': max_M,
            'exec_time': exec_time,
            'max_error': errors['max_error'],
            'avg_error': errors['avg_error'],
            'rms_error': errors['rms_error']
        })
        print(f"  N={N}, exec_time={exec_time:.4f}s, max_error={errors['max_error']:.4e}")

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
    plot_error_vs_parameter(results_df, 'N', output_dir, f" ($M={max_M}$, $N_o={max_N_o}$)")
    plot_exec_time_vs_parameter(results_df, 'N', output_dir, f" ($M={max_M}$, $N_o={max_N_o}$)")

    results = []
    for M in M_values:
        # Run simulation for current M and max_N
        sim_df, exec_time, _ = run_monte_carlo(args.netlist_path, 0, max_N, M, print_step=max_N_o)
        mc_results[(max_N, M, max_N_o)] = sim_df  # Store in dictionary with consistent 3-tuple format

        errors = calculate_errors(sim_df, ref_df)
        results.append({
            'N': max_N,
            'M': M,
            'exec_time': exec_time,
            'max_error': errors['max_error'],
            'avg_error': errors['avg_error'],
            'rms_error': errors['rms_error']
        })
        print(f"  M={M}, exec_time={exec_time:.4f}s, max_error={errors['max_error']:.4e}")

    # Store results as DataFrame
    results_df = pd.DataFrame(results)
    # Output directory should already exist from earlier steps
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)  # Safe to call multiple times
    results_df.to_csv(output_dir / 'sweep_M_results.csv', index=False)

    # Plot results for M sweep
    plot_error_vs_parameter(results_df, 'M', output_dir, f" ($N={max_N}$, $N_o={max_N_o}$)")
    plot_exec_time_vs_parameter(results_df, 'M', output_dir, f" ($N={max_N}$, $N_o={max_N_o}$)")

    # Print Step Sweep Analysis
    print("\nRunning print step sweep analysis...")
    results = []
    for print_step in N_o_values:
        # Run simulation for current print_step with max_N and max_M
        sim_df, exec_time, _ = run_monte_carlo(args.netlist_path, 0, max_N, max_M, print_step=print_step)
        mc_results[(max_N, max_M, print_step)] = sim_df  # Store in dictionary with print_step as third key

        errors = calculate_errors(sim_df, ref_df)
        results.append({
            'N_o': print_step,
            'N': max_N,
            'M': max_M,
            'exec_time': exec_time,
            'max_error': errors['max_error'],
            'avg_error': errors['avg_error'],
            'rms_error': errors['rms_error']
        })
        print(f"  print_step={print_step}, exec_time={exec_time:.4f}s, max_error={errors['max_error']:.4e}")

    # Store results as DataFrame
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_dir / 'sweep_print_step_results.csv', index=False)

    # Plot results for print step sweep
    plot_error_vs_parameter(results_df, 'N_o', output_dir, f" ($N={max_N}$, $M={max_M}$)")
    plot_exec_time_vs_parameter(results_df, 'N_o', output_dir, f" ($N={max_N}$, $M={max_M}$)")

    # Prepare plot data
    plot_data = [
        {
            'df': ref_df,
            'label': f'Trapezoidal ($N={N_o}$)',
            'style': {'color': 'C1', 'linestyle': '-', 'linewidth': 2}
        },
        {
            'df': mc_results[(max_N, max_M, max_N_o)],
            'label': f'Best MC ($N={max_N}$, $M={max_M}$, $N_o={max_N_o}$)',
            'style': {'color': 'C0', 'linestyle': '-', 'linewidth': 2}
        },
        {
            'df': mc_results[(max_N, min_M, max_N_o)],
            'label': f'Low $M$ ($N={max_N}$, $M={min_M}$, $N_o={max_N_o}$)',
            'style': {'color': 'C2', 'linestyle': '--', 'linewidth': 1}
        },
        {
            'df': mc_results[(min_N, max_M, max_N_o)],
            'label': f'Low $N$ ($N={min_N}$, $M={max_M}$, $N_o={max_N_o}$)',
            'style': {'color': 'C3', 'linestyle': ':', 'linewidth': 1}
        },
        {
            'df': mc_results[(max_N, max_M, min_N_o)],
            'label': f'Low $N_o$ ($N={max_N}$, $M={max_M}$, $N_o={min_N_o}$)',
            'style': {'color': 'C4', 'linestyle': '-.', 'linewidth': 1}
        }
    ]

    # Add SPICE data if available
    if spice_df is not None:
        plot_data.append({
            'df': spice_df,
            'label': 'Ngspice',
            'style': {'color': 'C5', 'linestyle': '-', 'linewidth': 2, 'alpha': 0.7}
        })

    # Plot solution comparison showing actual data
    plot_voltage_comparison(
        plot_data=plot_data,
        output_dir=output_dir,
        filename='solution_comparison.pdf',
        title='Solution Comparison: Effect of N and M Parameters'
    )

    # Run scalability analysis
    # run_scalability_analysis(args.netlist_path, output_dir)
    
    # Run time vs error analysis for both trapezoidal and Monte Carlo methods
    # Use SPICE as reference if available, otherwise use highest N trapezoidal
    reference_for_analysis = spice_df if spice_df is not None and not spice_df.empty else None
    plot_time_vs_error_analysis(args.netlist_path, 0, parameter_values=parameter_values, output_dir=output_dir, ref_df=reference_for_analysis)

if __name__ == "__main__":
    main()
