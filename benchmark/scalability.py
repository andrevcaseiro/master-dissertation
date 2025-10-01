#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing
from tran import run_monte_carlo

# Global parameter values
N_values = [100000, 1000000, 10000000]
M_values = [100, 1000, 10000]


def plot_speedup_analysis(speedup_df, output_dir):
    """
    Plot speedup analysis vs number of threads.

    Args:
        speedup_df (pandas.DataFrame): Speedup results dataframe
        output_dir (Path): Directory to save plots
    """
    plt.figure()
    plt.plot(speedup_df["threads"], speedup_df["speedup"], marker="o", label="Measured Speedup")
    plt.plot(speedup_df["threads"], speedup_df["threads"], "k--", label="Ideal Linear Speedup")
    plt.xscale("log", base=2)
    plt.yscale("log", base=2)
    plt.xlabel("Threads (log₂ scale)")
    plt.ylabel("Speedup (log₂ scale)")
    # plt.title('Scalability: Speedup vs Threads')
    plt.legend()
    plt.grid(True, which="both", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(output_dir / "speedup_vs_threads.pdf", bbox_inches="tight")
    plt.close()


def run_scalability_analysis(netlist_path, output_dir):
    """
    Run scalability (speedup) analysis by varying number of threads.

    Args:
        netlist_path (str): Path to the SPICE netlist file
        output_dir (Path): Directory to save results

    Returns:
        pandas.DataFrame: DataFrame with scalability results
    """
    max_N = max(N_values)
    max_M = max(M_values)

    # Scalability (speedup) analysis: vary number of threads from 1 to available CPU cores
    max_threads = multiprocessing.cpu_count()
    thread_values = [2**i for i in range(int(np.log2(max_threads)) + 1) if 2**i <= max_threads]
    if 1 not in thread_values:
        thread_values = [1] + thread_values  # Ensure 1 is included

    # Use max_N and max_M for a fair comparison
    baseline_time = None
    speedup_results = []

    for threads in thread_values:
        # Assume run_monte_carlo accepts a 'threads' argument
        sim_df, exec_time, _ = run_monte_carlo(netlist_path, 0, max_N, max_M, print_step=10000, num_threads=threads)
        if not baseline_time:
            baseline_time = exec_time
        speedup = baseline_time / exec_time
        speedup_results.append({"threads": threads, "exec_time": exec_time, "speedup": speedup})
        print(f"Threads={threads}, exec_time={exec_time:.4f}s, speedup={speedup:.2f}")

    speedup_df = pd.DataFrame(speedup_results)
    speedup_df.to_csv(output_dir / "scalability_results.csv", index=False)

    # Plot speedup analysis
    plot_speedup_analysis(speedup_df, output_dir)

    return speedup_df
