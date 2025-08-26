#!/usr/bin/env python3
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import multiprocessing
from tran import run_monte_carlo, run_trapezoidal, run_ngspice


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
    Plot error metrics vs a parameter (N or M).
    
    Args:
        results_df (pandas.DataFrame): Results dataframe with error metrics
        parameter_name (str): Name of the parameter ('N' or 'M')
        output_dir (Path): Directory to save plots
        fixed_param_info (str): Information about fixed parameter for title
    """
    plt.figure()
    plt.plot(results_df[parameter_name], results_df['max_error'], marker='o', label='Max Error')
    plt.plot(results_df[parameter_name], results_df['avg_error'], marker='x', label='Avg Error')
    plt.plot(results_df[parameter_name], results_df['rms_error'], marker='s', label='RMS Error')
    plt.xscale('log', base=2)
    plt.yscale('log', base=2)
    plt.xlabel(f'{parameter_name} (log₂ scale)')
    plt.ylabel('Error (log₂ scale)')
    plt.title(f'Error vs {parameter_name}{fixed_param_info}')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.savefig(output_dir / f'error_vs_{parameter_name}.png')
    plt.close()

def plot_exec_time_vs_parameter(results_df, parameter_name, output_dir, fixed_param_info=""):
    """
    Plot execution time vs a parameter (N or M).
    
    Args:
        results_df (pandas.DataFrame): Results dataframe with execution times
        parameter_name (str): Name of the parameter ('N' or 'M')
        output_dir (Path): Directory to save plots
        fixed_param_info (str): Information about fixed parameter for title
    """
    plt.figure()
    plt.plot(results_df[parameter_name], results_df['exec_time'], marker='o')
    plt.xscale('log', base=2)
    plt.xlabel(f'{parameter_name} (log₂ scale)')
    plt.ylabel('Execution Time (s)')
    plt.title(f'Execution Time vs {parameter_name}{fixed_param_info}')
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.savefig(output_dir / f'exec_time_vs_{parameter_name}.png')
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
    plt.title('Scalability: Speedup vs Threads')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.savefig(output_dir / 'speedup_vs_threads.png')
    plt.show()
    plt.close()

def plot_solution_comparison(ref_df, mc_results, output_dir):
    """
    Plot comparison of actual solution data for different parameter combinations.
    
    Args:
        ref_df (pandas.DataFrame): Reference solution data (already computed)
        mc_results (dict): Dictionary with (N,M) as key and DataFrame as value
        output_dir (Path): Directory to save plots
    """
    # Extract N and M values from dictionary keys
    N_values = sorted(set(key[0] for key in mc_results.keys()))
    M_values = sorted(set(key[1] for key in mc_results.keys()))
    
    min_N, max_N = min(N_values), max(N_values)
    min_M, max_M = min(M_values), max(M_values)
    
    plt.figure(figsize=(12, 8))
    
    # Extract the specific solutions we want to plot
    best_mc_df = mc_results[(max_N, max_M)]
    worst_m_df = mc_results[(max_N, min_M)]
    worst_n_df = mc_results[(min_N, max_M)]
    
    # Plot all solutions
    plt.plot(ref_df['time'], ref_df['voltage'], 'k-', linewidth=2, label='Reference (Trapezoidal)', alpha=0.8)
    plt.plot(best_mc_df['time'], best_mc_df['voltage'], 'g-', linewidth=1.5, label=f'Best MC (N={max_N}, M={max_M})', alpha=0.7)
    plt.plot(worst_m_df['time'], worst_m_df['voltage'], 'r--', linewidth=1.5, label=f'Low M (N={max_N}, M={min_M})', alpha=0.7)
    plt.plot(worst_n_df['time'], worst_n_df['voltage'], 'b:', linewidth=1.5, label=f'Low N (N={min_N}, M={max_M})', alpha=0.7)
    
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (V)')
    plt.title('Solution Comparison: Effect of N and M Parameters')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(output_dir / 'solution_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Run parameter sweeps for ODE solvers and plot results')
    parser.add_argument('netlist_path', help='Path to the SPICE netlist file')
    parser.add_argument('--output-dir', type=str, help='Directory to save plots (default: benchmark_dir/\{netlist_name\})')
    
    args = parser.parse_args()
    
    # Set default output directory if not provided
    if args.output_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        netlist_name = os.path.splitext(os.path.basename(args.netlist_path))[0]
        args.output_dir = os.path.join(script_dir, netlist_name)

    # Calculate reference df based on highest values of N and M
    N_values = [1024 * (2 ** i) for i in range(8)]
    M_values = [1024 * (2 ** i) for i in range(5)]
    max_N = max(N_values)
    max_M = max(M_values)

    # Run the simulation with the highest N and M to get the reference data
    # Assuming run_trapezoidal returns a DataFrame with 'time' and 'voltage' columns
    ref_df, _ = run_trapezoidal(args.netlist_path, 0, max_N)

    # Store all Monte Carlo results in a dictionary with (N,M) as key
    mc_results = {}
    
    results = []
    for N in N_values:
        # Run simulation for current N and max_M
        sim_df, exec_time = run_monte_carlo(args.netlist_path, 0, N, max_M, print_step=1024)
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
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_dir / 'sweep_N_results.csv', index=False)

    # Plot results for N sweep
    plot_error_vs_parameter(results_df, 'N', output_dir, " (M fixed at max)")
    plot_exec_time_vs_parameter(results_df, 'N', output_dir, " (M fixed at max)")

    results = []
    for M in M_values:
        # Run simulation for current M and max_N
        sim_df, exec_time = run_monte_carlo(args.netlist_path, 0, max_N, M, print_step=1024)
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
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_dir / 'sweep_M_results.csv', index=False)

    # Plot results for M sweep
    plot_error_vs_parameter(results_df, 'M', output_dir, " (N fixed at max)")
    plot_exec_time_vs_parameter(results_df, 'M', output_dir, " (N fixed at max)")

    # Plot solution comparison showing actual data
    plot_solution_comparison(ref_df, mc_results, output_dir)
    
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
        sim_df, exec_time = run_monte_carlo(args.netlist_path, 0, max_N, max_M, print_step=1024, num_threads=threads)
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
