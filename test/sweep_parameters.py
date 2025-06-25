#!/usr/bin/env python3
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import multiprocessing
from test_tran import run_monte_carlo, run_trapezoidal, run_ngspice


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

def main():
    parser = argparse.ArgumentParser(description='Run parameter sweeps for ODE solvers and plot results')
    parser.add_argument('netlist_path', help='Path to the SPICE netlist file')
    parser.add_argument('--output-dir', type=str, default='./res/sweep', help='Directory to save plots')
    
    args = parser.parse_args()

    N_values = [1024 * (2 ** i) for i in range(5)]
    M_values = [1024 * (2 ** i) for i in range(5)]

    # Calculate reference df based on highest values of N and M
    max_N = max(N_values)
    max_M = max(M_values)

    # Run the simulation with the highest N and M to get the reference data
    # Assuming run_trapezoidal returns a DataFrame with 'time' and 'voltage' columns
    ref_df, _ = run_trapezoidal(args.netlist_path, 0, max_N)

    results = []
    for N in N_values:
        # Run simulation for current N and max_M
        sim_df, exec_time = run_monte_carlo(args.netlist_path, 0, N, max_M, print_step=1024)
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

    # Plot error vs N
    plt.figure()
    plt.plot(results_df['N'], results_df['max_error'], marker='o', label='Max Error')
    plt.plot(results_df['N'], results_df['avg_error'], marker='x', label='Avg Error')
    plt.plot(results_df['N'], results_df['rms_error'], marker='s', label='RMS Error')
    plt.xscale('log', base=2)  # Use logarithmic scale base 2 for x-axis
    plt.yscale('log', base=2)  # Log scale base 2 for errors
    plt.xlabel('N (log₂ scale)')
    plt.ylabel('Error (log₂ scale)')
    plt.title('Error vs N (M fixed at max)')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)  # Grid for both major and minor ticks
    plt.savefig(output_dir / 'error_vs_N.png')
    plt.close()

    # Plot execution time vs N
    plt.figure()
    plt.plot(results_df['N'], results_df['exec_time'], marker='o')
    plt.xscale('log', base=2)  # Use logarithmic scale base 2 for x-axis
    plt.xlabel('N (log₂ scale)')
    plt.ylabel('Execution Time (s)')
    plt.title('Execution Time vs N (M fixed at max)')
    plt.grid(True, which='both', linestyle='--', alpha=0.7)  # Grid for both major and minor ticks
    plt.savefig(output_dir / 'exec_time_vs_N.png')
    plt.close()

    results = []
    for M in M_values:
        # Run simulation for current M and max_N
        sim_df, exec_time = run_monte_carlo(args.netlist_path, 0, max_N, M, print_step=1024)
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

    # Plot error vs M
    plt.figure()
    plt.plot(results_df['M'], results_df['max_error'], marker='o', label='Max Error')
    plt.plot(results_df['M'], results_df['avg_error'], marker='x', label='Avg Error')
    plt.plot(results_df['M'], results_df['rms_error'], marker='s', label='RMS Error')
    plt.xscale('log', base=2)  # Use logarithmic scale base 2 for x-axis
    plt.yscale('log', base=2)  # Log scale base 2 for errors
    plt.xlabel('M (log₂ scale)')
    plt.ylabel('Error (log₂ scale)')
    plt.title('Error vs M (N fixed at max)')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.7)  # Grid for both major and minor ticks
    plt.savefig(output_dir / 'error_vs_M.png')
    plt.close()

    # Plot execution time vs M
    plt.figure()
    plt.plot(results_df['M'], results_df['exec_time'], marker='o')
    plt.xscale('log', base=2)  # Use logarithmic scale base 2 for x-axis
    plt.xlabel('M (log₂ scale)')
    plt.ylabel('Execution Time (s)')
    plt.title('Execution Time vs M (N fixed at max)')
    plt.grid(True, which='both', linestyle='--', alpha=0.7)  # Grid for both major and minor ticks
    plt.savefig(output_dir / 'exec_time_vs_M.png')
    plt.close()


    # Scalability (speedup) analysis: vary number of threads from 1 to available CPU cores
    max_threads = multiprocessing.cpu_count()
    thread_values = [2 ** i for i in range(int(np.log2(max_threads)) + 1) if 2 ** i <= max_threads]
    if 1 not in thread_values:
        thread_values = [1] + thread_values  # Ensure 1 is included

    speedup_results = []
    # Use max_N and max_M for a fair comparison
    baseline_time = None

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

    # Plot speedup
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
    
    

if __name__ == "__main__":
    main()
