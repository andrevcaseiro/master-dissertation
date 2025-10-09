#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path
from tran import run_trapezoidal, run_monte_carlo
from utils import calculate_errors


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


def run_time_vs_error_analysis(netlist_path, final_time, parameter_values, ref_df=None):
    """
    Run simulations with different parameter values and collect time/error data.
    
    Args:
        netlist_path (str): Path to the SPICE netlist file
        final_time (float): Final simulation time
        parameter_values (dict): Dictionary with keys:
            - 'trap': List of N values for trapezoidal method
            - 'mc': List of tuples (n_o, n, m) for Monte Carlo method
        ref_df (pd.DataFrame, optional): Reference solution to use for error calculation.
                                        If None, uses highest N trapezoidal solution.
        
    Returns:
        tuple: (trap_df, mc_df) - DataFrames for trapezoidal and Monte Carlo results
    """
    trap_results = []
    mc_results = []
    
    # Get reference solution
    if ref_df is None:
        # Use highest N value from trapezoidal method as reference
        max_n = max(parameter_values['trap'])
        print(f"Running reference solution with highest N={max_n}...")
        ref_df, _, _ = run_trapezoidal(
            netlist_path=netlist_path,
            final_time=final_time,
            num_steps=max_n,
            method="pardiso"
        )
        
        if ref_df.empty:
            print("Error: Could not get reference solution with highest N")
            return pd.DataFrame(), pd.DataFrame()
    else:
        max_n = max(parameter_values['trap'])
        print("Using provided reference solution (SPICE)")
    
    # Run trapezoidal simulations
    for n in parameter_values['trap']:
        try:
            sim_df, exec_time, _ = run_trapezoidal(
                netlist_path=netlist_path,
                final_time=final_time,
                num_steps=n,
                method="pardiso"
            )
            
            if not sim_df.empty and exec_time is not None:

                error_metrics = calculate_errors(sim_df, ref_df)
                
                trap_results.append({
                    'N': n,
                    'exec_time': exec_time,
                    'max_abs_error': error_metrics['max_error'],
                    'max_rel_error': error_metrics['max_rel_error'],
                    'avg_abs_error': error_metrics['avg_error'],
                    'avg_rel_error': error_metrics['avg_rel_error'],
                    'rms_abs_error': error_metrics['rms_error'],
                    'rms_rel_error': error_metrics['rms_rel_error']
                })
                
                print(f"  Trapezoidal N={n}: time={exec_time:.4f}s, rms_rel_error={error_metrics['rms_rel_error']:.2e}")
        
        except Exception as e:
            print(f"Error running trapezoidal simulation for N={n}: {e}")
    
    # Run Monte Carlo simulations
    for n_o, n, m in parameter_values['mc']:
        try:
            sim_df, exec_time, _ = run_monte_carlo(
                netlist_path=netlist_path,
                final_time=final_time,
                num_steps=n,
                num_samples=m,
                print_step=n_o,
            )
            
            if not sim_df.empty and exec_time is not None:
                error_metrics = calculate_errors(sim_df, ref_df)
                
                mc_results.append({
                    'n_o': n_o,
                    'N': n,
                    'M': m,
                    'exec_time': exec_time,
                    'max_abs_error': error_metrics['max_error'],
                    'max_rel_error': error_metrics['max_rel_error'],
                    'avg_abs_error': error_metrics['avg_error'],
                    'avg_rel_error': error_metrics['avg_rel_error'],
                    'rms_abs_error': error_metrics['rms_error'],
                    'rms_rel_error': error_metrics['rms_rel_error']
                })
                
                print(f"  Monte Carlo n_o={n_o}, N={n}, M={m}: time={exec_time:.4f}s, rms_rel_error={error_metrics['rms_rel_error']:.2e}")
        
        except Exception as e:
            print(f"Error running Monte Carlo simulation for n_o={n_o}, N={n}, M={m}: {e}")
    
    return pd.DataFrame(trap_results), pd.DataFrame(mc_results)


def plot_time_vs_error(trap_df, mc_df, netlist_name, output_dir, error_column, ylabel):
    """
    Plot execution time vs error for both trapezoidal and Monte Carlo methods.
    
    Args:
        trap_df (pd.DataFrame): Trapezoidal results with all error metric columns
        mc_df (pd.DataFrame): Monte Carlo results with all error metric columns
        netlist_name (str): Name of the netlist for the plot title
        output_dir (Path): Directory to save the plot
        error_column (str): Column name for the error metric to plot
        ylabel (str): Label for the y-axis
    """
    if trap_df.empty and mc_df.empty:
        print("No results to plot")
        return
    
    plt.figure(figsize=(10, 6))
    
    # Plot trapezoidal method results with squares
    if not trap_df.empty:
        trap_df = trap_df.sort_values('exec_time')
        plt.scatter(trap_df['exec_time'], trap_df[error_column], 
                   color='C1', s=60, marker='s', label='Trapezoidal')
        
        # Add N value annotations for trapezoidal
        for _, row in trap_df.iterrows():
            plt.annotate(f'$N={format_scientific_notation(row["N"])}$', 
                        (row['exec_time'], row[error_column]),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=10, alpha=0.8)
    
    # Plot Monte Carlo method results with color scale based on N_o
    if not mc_df.empty:
        mc_df = mc_df.sort_values('exec_time')
        
        # Create color map based on N_o values
        n_o_values = mc_df['n_o'].values
        if len(set(n_o_values)) > 1:  # Only use colormap if there are different N_o values
            # Use logarithmic normalization for color mapping
            from matplotlib.colors import LogNorm
            
            scatter = plt.scatter(mc_df['exec_time'], mc_df[error_column],
                               c=n_o_values, cmap='viridis', s=60, label='Monte Carlo',
                               norm=LogNorm(vmin=n_o_values.min(), vmax=n_o_values.max()))
            
            # Add colorbar for N_o values with logarithmic scale
            cbar = plt.colorbar(scatter)
            cbar.set_label('$N_o$')
        else:
            # Use single color if all N_o values are the same
            plt.scatter(mc_df['exec_time'], mc_df[error_column],
                       color='C0', s=60, label='Monte Carlo')
        
        # Add n_o,N,M value annotations for Monte Carlo
        for _, row in mc_df.iterrows():
            n_o_str = format_scientific_notation(row["n_o"])
            n_str = format_scientific_notation(row["N"])
            m_str = format_scientific_notation(row["M"])
            plt.annotate(f'$N_o={n_o_str},N={n_str},M={m_str}$', 
                        (row['exec_time'], row[error_column]),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8, alpha=0.8)
    
    plt.xlabel('Execution Time (s)')
    plt.ylabel(ylabel)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / f'time_vs_{error_column}.pdf', bbox_inches='tight')
    plt.show()


def plot_time_vs_error_analysis(netlist_path, final_time, parameter_values={}, output_dir=None, ref_df=None):
    """
    Run the complete time vs error analysis.
    
    Args:
        netlist_path (str): Path to the SPICE netlist file
        final_time (float): Final simulation time
        parameter_values (dict): Dictionary with keys:
            - 'trap': List of N values for trapezoidal method
            - 'mc': List of tuples (n_o, n, m) for Monte Carlo method
        output_dir (Path): Directory to save the plot
        ref_df (pd.DataFrame, optional): Reference solution to use for error calculation.
                                        If None, uses highest N trapezoidal solution.
    """
    # Run simulations
    trap_df, mc_df = run_time_vs_error_analysis(netlist_path, final_time, parameter_values, ref_df)
    
    # Plot results for all error metrics
    netlist_name = Path(netlist_path).stem
    if output_dir is None:
        output_dir = Path('.')
    
    print(trap_df)
    print(mc_df)

    # Define all error metrics to plot
    error_metrics = [
        ('max_abs_error', 'Maximum Absolute Error'),
        ('max_rel_error', 'Maximum Relative Error'),
        ('avg_abs_error', 'Average Absolute Error'),
        ('avg_rel_error', 'Average Relative Error'),
        ('rms_abs_error', 'RMS Absolute Error'),
        ('rms_rel_error', 'RMS Relative Error')
    ]
    
    # Generate plots for all error metrics
    for error_column, ylabel in error_metrics:
        plot_time_vs_error(trap_df, mc_df, netlist_name, output_dir, error_column, ylabel)