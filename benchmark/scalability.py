#!/usr/bin/env python3
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing
from tran import run_monte_carlo, run_trapezoidal
from utils import N_values, M_values, N_o_values, fit_amdahls_law
from scipy.optimize import curve_fit



def plot_speedup_analysis(mc_speedup_df, trap_speedup_df, output_dir):
    """
    Plot speedup analysis vs number of threads for both Monte Carlo and Trapezoidal methods.

    Args:
        mc_speedup_df (pandas.DataFrame): Monte Carlo speedup results dataframe
        trap_speedup_df (pandas.DataFrame): Trapezoidal speedup results dataframe
        output_dir (Path): Directory to save plots
    """
    plt.figure(figsize=(4, 3))
    
    # Plot ideal linear speedup
    threads = mc_speedup_df["threads"]
    plt.plot(threads, threads, "k-", label="Ideal Linear Speedup", lw=1)

    # Plot Monte Carlo speedup
    mc_x = mc_speedup_df["threads"].values
    mc_y = mc_speedup_df["speedup"].values
    plt.plot(mc_x, mc_y, "C0o", label="Monte Carlo Speedup")
    
    # Fit and plot Monte Carlo Amdahl's Law
    mc_fit = fit_amdahls_law(mc_x, mc_y)
    plt.plot(mc_fit['threads_smooth'], mc_fit['speedup_smooth'], "C0--", lw=1.5, alpha=0.7,
             label=f'MC Amdahl\'s Law: s={mc_fit["s"]*100:.1f}% (R²={mc_fit["r_squared"]:.3f})')
    
    # Plot Trapezoidal speedup
    trap_x = trap_speedup_df["threads"].values
    trap_y = trap_speedup_df["speedup"].values
    plt.plot(trap_x, trap_y, "C1s", label="Trapezoidal Speedup")
    
    # Fit and plot Trapezoidal Amdahl's Law
    trap_fit = fit_amdahls_law(trap_x, trap_y)
    plt.plot(trap_fit['threads_smooth'], trap_fit['speedup_smooth'], "C1--", lw=1.5, alpha=0.7,
             label=f'Trap Amdahl\'s Law: s={trap_fit["s"]*100:.1f}% (R²={trap_fit["r_squared"]:.3f})')
    
    plt.xscale("log", base=2)
    plt.yscale("log", base=2)
    for axis in [plt.gca().xaxis, plt.gca().yaxis]:
        axis.set_major_formatter(FuncFormatter(lambda x, p: f'{int(x)}'))
    plt.xlabel("Threads (log₂ scale)")
    plt.ylabel("Speedup (log₂ scale)")
    # plt.title('Scalability: Speedup vs Threads')
    plt.legend()
    plt.grid(True, which="both", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_dir / "speedup_vs_threads.pdf", bbox_inches="tight")
    plt.close()
    
    # Print detailed analysis
    print(f"\n=== Amdahl's Law Analysis ===")
    print(f"Monte Carlo:")
    print(f"  Serial fraction (s): {mc_fit['s']:.3f} ({mc_fit['s']*100:.1f}%)")
    print(f"  Parallel fraction: {mc_fit['parallel_fraction']:.3f} ({mc_fit['parallel_fraction']*100:.1f}%)")
    print(f"  R² goodness of fit: {mc_fit['r_squared']:.3f}")
    
    print(f"\nTrapezoidal:")
    print(f"  Serial fraction (s): {trap_fit['s']:.3f} ({trap_fit['s']*100:.1f}%)")
    print(f"  Parallel fraction: {trap_fit['parallel_fraction']:.3f} ({trap_fit['parallel_fraction']*100:.1f}%)")
    print(f"  R² goodness of fit: {trap_fit['r_squared']:.3f}")
    
    return mc_fit, trap_fit


def fit_execution_time_model(threads, exec_time):
    """
    Fit execution time model: T(p) = T_serial + T_parallel/p
    
    Args:
        threads (array-like): Number of threads
        exec_time (array-like): Execution times
        
    Returns:
        dict: Dictionary containing fitted parameters and analysis
    """
    def exec_time_model(p, t_serial, t_parallel):
        return t_serial + t_parallel / p
    
    try:
        # Fit the model
        popt, pcov = curve_fit(exec_time_model, threads, exec_time, 
                              p0=[0.1 * exec_time[0], 0.9 * exec_time[0]], 
                              bounds=([0, 0], [exec_time[0], 10 * exec_time[0]]))
        
        t_serial, t_parallel = popt
        
        # Calculate R-squared
        exec_time_pred = exec_time_model(threads, t_serial, t_parallel)
        ss_res = np.sum((exec_time - exec_time_pred) ** 2)
        ss_tot = np.sum((exec_time - np.mean(exec_time)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        # Generate smooth curve for plotting
        threads_smooth = np.logspace(np.log10(threads.min()), np.log10(threads.max()), 100)
        exec_time_smooth = exec_time_model(threads_smooth, t_serial, t_parallel)
        
        # Calculate parallel efficiency
        parallel_fraction = t_parallel / (t_serial + t_parallel)
        serial_fraction = t_serial / (t_serial + t_parallel)
        
        return {
            't_serial': t_serial,
            't_parallel': t_parallel,
            'serial_fraction': serial_fraction,
            'parallel_fraction': parallel_fraction,
            'r_squared': r_squared,
            'threads_smooth': threads_smooth,
            'exec_time_smooth': exec_time_smooth,
            'fitted_exec_time': exec_time_pred
        }
        
    except Exception as e:
        print(f"Error fitting execution time model: {e}")
        return {
            't_serial': exec_time[0] * 0.5,
            't_parallel': exec_time[0] * 0.5,
            'serial_fraction': 0.5,
            'parallel_fraction': 0.5,
            'r_squared': 0.0,
            'threads_smooth': threads,
            'exec_time_smooth': exec_time,
            'fitted_exec_time': exec_time
        }


def plot_execution_time_analysis(mc_speedup_df, trap_speedup_df, output_dir):
    """
    Plot execution time analysis vs number of threads for both Monte Carlo and Trapezoidal methods.

    Args:
        mc_speedup_df (pandas.DataFrame): Monte Carlo speedup results dataframe
        trap_speedup_df (pandas.DataFrame): Trapezoidal speedup results dataframe
        output_dir (Path): Directory to save plots
    """
    plt.figure(figsize=(4, 3))
    
    # Plot Monte Carlo execution time
    mc_x = mc_speedup_df["threads"].values
    mc_y = mc_speedup_df["exec_time"].values
    plt.plot(mc_x, mc_y, "C0o", label="Monte Carlo")
    
    # Fit and plot Monte Carlo execution time model: T(p) = T_serial + T_parallel/p
    mc_fit = fit_execution_time_model(mc_x, mc_y)
    plt.plot(mc_fit['threads_smooth'], mc_fit['exec_time_smooth'], "C0--", lw=1.5, alpha=0.7,
             label=f'MC Model: $T_s={mc_fit["t_serial"]:.3f}\\text{{s}}, T_p={mc_fit["t_parallel"]:.3f}\\text{{s}}\quad (R^2={mc_fit["r_squared"]:.3f})$')
    
    # Plot Trapezoidal execution time
    trap_x = trap_speedup_df["threads"].values
    trap_y = trap_speedup_df["exec_time"].values
    plt.plot(trap_x, trap_y, "C1s", label="Trapezoidal")
    
    # Fit and plot Trapezoidal execution time model
    trap_fit = fit_execution_time_model(trap_x, trap_y)
    plt.plot(trap_fit['threads_smooth'], trap_fit['exec_time_smooth'], "C1--", lw=1.5, alpha=0.7,
             label=f'Trap Model: $T_s={trap_fit["t_serial"]:.3f}\\text{{s}}, T_p={trap_fit["t_parallel"]:.3f}\\text{{s}}\quad (R^2={trap_fit["r_squared"]:.3f})$')
    
    plt.xscale("log", base=2)
    plt.yscale("log")
    for axis in [plt.gca().xaxis, plt.gca().yaxis]:
        axis.set_major_formatter(FuncFormatter(lambda x, p: f'{int(x)}' if x >= 1 else f'{x:.1f}'))
    plt.xlabel("Threads (log₂ scale)")
    plt.ylabel("Execution Time (seconds, log scale)")
    plt.legend()
    plt.grid(True, which="both", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_dir / "execution_time_vs_threads.pdf", bbox_inches="tight")
    plt.close()
    
    # Print detailed analysis
    print(f"\n=== Execution Time Model Analysis: T(p) = T_serial + T_parallel/p ===")
    print(f"Monte Carlo:")
    print(f"  Serial time (T_serial): {mc_fit['t_serial']:.3f}s ({mc_fit['serial_fraction']*100:.1f}%)")
    print(f"  Parallel time (T_parallel): {mc_fit['t_parallel']:.3f}s ({mc_fit['parallel_fraction']*100:.1f}%)")
    print(f"  R² goodness of fit: {mc_fit['r_squared']:.3f}")
    print(f"  T(p) = {mc_fit['t_serial']:.3f} + {mc_fit['t_parallel']:.3f}/p")
    
    print(f"\nTrapezoidal:")
    print(f"  Serial time (T_serial): {trap_fit['t_serial']:.3f}s ({trap_fit['serial_fraction']*100:.1f}%)")
    print(f"  Parallel time (T_parallel): {trap_fit['t_parallel']:.3f}s ({trap_fit['parallel_fraction']*100:.1f}%)")
    print(f"  R² goodness of fit: {trap_fit['r_squared']:.3f}")
    print(f"  T(p) = {trap_fit['t_serial']:.3f} + {trap_fit['t_parallel']:.3f}/p")
    
    return mc_fit, trap_fit


def run_scalability_analysis(netlist_path, output_dir):
    """
    Run scalability (speedup) analysis by varying number of threads for both Monte Carlo and Trapezoidal methods.

    Args:
        netlist_path (str): Path to the SPICE netlist file
        output_dir (Path): Directory to save results

    Returns:
        tuple: (mc_speedup_df, trap_speedup_df) DataFrames with scalability results
    """
    max_N = N_values[-1]
    max_M = M_values[-1]
    max_N_o = N_o_values[-1]

    # Scalability (speedup) analysis: vary number of threads from 1 to available CPU cores
    max_threads = multiprocessing.cpu_count() / 2 # assume two virtual threads per phyical core
    thread_values = [2**i for i in range(int(np.log2(max_threads)) + 1) if 2**i <= max_threads]
    if 1 not in thread_values:
        thread_values = [1] + thread_values  # Ensure 1 is included

    # Monte Carlo scalability analysis
    print("Running Monte Carlo scalability analysis...")
    mc_baseline_time = None
    mc_speedup_results = []

    for threads in thread_values:
        sim_df, exec_time, _ = run_monte_carlo(netlist_path, 0, max_N, max_M, print_step=max_N_o, num_threads=threads)
        if not mc_baseline_time:
            mc_baseline_time = exec_time
        speedup = mc_baseline_time / exec_time
        mc_speedup_results.append({"threads": threads, "exec_time": exec_time, "speedup": speedup})
        print(f"  Monte Carlo - Threads={threads}, exec_time={exec_time:.4f}s, speedup={speedup:.2f}")

    # Trapezoidal scalability analysis
    print("Running Trapezoidal scalability analysis...")
    trap_baseline_time = None
    trap_speedup_results = []
    
    # Use a reasonable number of steps for trapezoidal method
    num_steps = 10000

    for threads in thread_values:
        sim_df, exec_time, _ = run_trapezoidal(netlist_path, 0, num_steps, method="pardiso", num_threads=threads)
        if not trap_baseline_time:
            trap_baseline_time = exec_time
        speedup = trap_baseline_time / exec_time
        trap_speedup_results.append({"threads": threads, "exec_time": exec_time, "speedup": speedup})
        print(f"  Trapezoidal - Threads={threads}, exec_time={exec_time:.4f}s, speedup={speedup:.2f}")

    # Create DataFrames
    mc_speedup_df = pd.DataFrame(mc_speedup_results)
    trap_speedup_df = pd.DataFrame(trap_speedup_results)
    
    # Save results
    mc_speedup_df.to_csv(output_dir / "monte_carlo_scalability_results.csv", index=False)
    trap_speedup_df.to_csv(output_dir / "trapezoidal_scalability_results.csv", index=False)

    # Plot combined speedup analysis with Amdahl's Law fitting
    mc_fit, trap_fit = plot_speedup_analysis(mc_speedup_df, trap_speedup_df, output_dir)
    
    # Plot execution time analysis with power law fitting
    mc_time_fit, trap_time_fit = plot_execution_time_analysis(mc_speedup_df, trap_speedup_df, output_dir)
    
    # Save Amdahl's Law analysis
    amdahl_analysis = pd.DataFrame([
        {
            'method': 'Monte Carlo',
            'serial_fraction': mc_fit['s'],
            'parallel_fraction': mc_fit['parallel_fraction'],
            'r_squared': mc_fit['r_squared']
        },
        {
            'method': 'Trapezoidal',
            'serial_fraction': trap_fit['s'],
            'parallel_fraction': trap_fit['parallel_fraction'],
            'r_squared': trap_fit['r_squared']
        }
    ])
    amdahl_analysis.to_csv(output_dir / "amdahl_analysis.csv", index=False)

    return mc_speedup_df, trap_speedup_df, mc_fit, trap_fit, mc_time_fit, trap_time_fit
