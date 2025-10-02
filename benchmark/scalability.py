#!/usr/bin/env python3
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing
from tran import run_monte_carlo, run_trapezoidal
from utils import fit_amdahls_law

# Global parameter values
N_values = [100000, 1000000, 10000000]
M_values = [100, 1000, 10000]


def plot_speedup_analysis(mc_speedup_df, trap_speedup_df, output_dir):
    """
    Plot speedup analysis vs number of threads for both Monte Carlo and Trapezoidal methods.

    Args:
        mc_speedup_df (pandas.DataFrame): Monte Carlo speedup results dataframe
        trap_speedup_df (pandas.DataFrame): Trapezoidal speedup results dataframe
        output_dir (Path): Directory to save plots
    """
    plt.figure(figsize=(8, 5))
    
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


def run_scalability_analysis(netlist_path, output_dir):
    """
    Run scalability (speedup) analysis by varying number of threads for both Monte Carlo and Trapezoidal methods.

    Args:
        netlist_path (str): Path to the SPICE netlist file
        output_dir (Path): Directory to save results

    Returns:
        tuple: (mc_speedup_df, trap_speedup_df) DataFrames with scalability results
    """
    max_N = max(N_values)
    max_M = max(M_values)

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
        sim_df, exec_time, _ = run_monte_carlo(netlist_path, 0, max_N, max_M, print_step=10000, num_threads=threads)
        if not mc_baseline_time:
            mc_baseline_time = exec_time
        speedup = mc_baseline_time / exec_time
        mc_speedup_results.append({"threads": threads, "exec_time": exec_time, "speedup": speedup})
        print(f"Monte Carlo - Threads={threads}, exec_time={exec_time:.4f}s, speedup={speedup:.2f}")

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
        print(f"Trapezoidal - Threads={threads}, exec_time={exec_time:.4f}s, speedup={speedup:.2f}")

    # Create DataFrames
    mc_speedup_df = pd.DataFrame(mc_speedup_results)
    trap_speedup_df = pd.DataFrame(trap_speedup_results)
    
    # Save results
    mc_speedup_df.to_csv(output_dir / "monte_carlo_scalability_results.csv", index=False)
    trap_speedup_df.to_csv(output_dir / "trapezoidal_scalability_results.csv", index=False)

    # Plot combined speedup analysis with Amdahl's Law fitting
    mc_fit, trap_fit = plot_speedup_analysis(mc_speedup_df, trap_speedup_df, output_dir)
    
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

    return mc_speedup_df, trap_speedup_df, mc_fit, trap_fit
