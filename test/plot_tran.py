import argparse
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path


def run_ngspice(circuit):
    """Run ngspice and return DataFrame with results."""
    result = subprocess.run(['ngspice', '-b', circuit], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError("ngspice failed!")
    
    data = []
    # Skip until we find the data section (starts with Index header)
    found_data = False
    for line in result.stdout.splitlines():
        if line.startswith("Index   time"):
            found_data = True
            continue
        if found_data and line.strip() and not line.startswith("---"):
            try:
                parts = line.strip().split()
                if len(parts) == 3:  # [index, time, voltage]
                    time = float(parts[1])
                    voltage = float(parts[2])
                    data.append((time, voltage))
            except (ValueError, IndexError):
                continue
    
    return pd.DataFrame(data, columns=["time", "voltage"])

def run_monte_carlo(circuit, M, N, threads=1):
    """Run Monte Carlo transient analysis and return DataFrame with results and execution time."""
    cmd = f"./main tran {circuit} -M {M} -N {N} -p 500".split()
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError("Monte Carlo simulation failed!")
    
    data = []
    exec_time = None
    for line in result.stdout.splitlines():
        if "Monte Carlo simulation time:" in line:
            exec_time = float(line.split(":")[1].split()[0])
            continue
        try:
            parts = line.strip().split()
            if len(parts) == 2:
                time = float(parts[0])
                voltage = float(parts[1])
                data.append((time, voltage))
        except (ValueError, IndexError):
            continue
    
    return pd.DataFrame(data, columns=["time", "voltage"]), exec_time

def run_trapezoidal(circuit, N, threads=1):
    """Run trapezoidal method and return DataFrame with results and execution time."""
    cmd = f"./main tran {circuit} --solver trapezoidal -N {N}".split()
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError("Trapezoidal simulation failed!")
    
    data = []
    exec_time = None
    for line in result.stdout.splitlines():
        if "Trapezoidal simulation time:" in line:
            exec_time = float(line.split(":")[1].split()[0])
            continue
        try:
            parts = line.strip().split()
            if len(parts) == 2:
                time = float(parts[0])
                voltage = float(parts[1])
                data.append((time, voltage))
        except (ValueError, IndexError):
            continue
    
    return pd.DataFrame(data, columns=["time", "voltage"]), exec_time

def calculate_errors(mc_df, ref_df):
    """Calculate error metrics between Monte Carlo and reference results."""
    # Interpolate reference data to match Monte Carlo time points
    ref_voltage = np.interp(mc_df['time'], ref_df['time'], ref_df['voltage'])
    
    abs_error = np.abs(mc_df['voltage'] - ref_voltage)
    rel_error = abs_error / (np.abs(ref_voltage) + 1e-10)  # Avoid division by zero
    
    return {
        'max_error': np.max(rel_error),
        'avg_error': np.mean(rel_error),
        'rms_error': np.sqrt(np.mean(rel_error**2)),
        'abs_error': abs_error,
        'ref_voltage': ref_voltage
    }

def plot_results(mc_df, trap_df, ref_df, errors, circuit_name):
    """Create comparison plot with error annotation."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Determine if trapezoidal is being used as reference
    use_trap_as_ref = ref_df is trap_df
    
    # Plot all methods
    if use_trap_as_ref:
        ax.plot(trap_df['time'], trap_df['voltage'], 'g-', label='Trapezoidal (Reference)', linewidth=1)
    else:
        ax.plot(trap_df['time'], trap_df['voltage'], 'r-', label='Trapezoidal', linewidth=1)
        ax.plot(ref_df['time'], ref_df['voltage'], 'g-', label='NgSpice (Reference)', linewidth=1)
    
    ax.plot(mc_df['time'], mc_df['voltage'], 'b-', label='Monte Carlo', linewidth=1)
    
    # Find and annotate maximum error
    max_error_idx = np.argmax(errors['abs_error'])
    max_time = mc_df['time'].iloc[max_error_idx]
    max_val = mc_df['voltage'].iloc[max_error_idx]
    
    # Determine annotation direction
    dir = 1 if mc_df['voltage'].iloc[max_error_idx] > errors['ref_voltage'][max_error_idx] else -1
    
    # Add error annotation
    ax.annotate(f'Max rel. error: {errors["max_error"]:.2e}',
                xy=(max_time, max_val),
                xytext=(max_time, max_val + dir * 0.07 * (mc_df['voltage'].max() - mc_df['voltage'].min())),
                ha='center',
                va='bottom' if dir == 1 else 'top',
                arrowprops=dict(arrowstyle="->", color='red'))
    
    ax.set_title(f'Transient Analysis Comparison - {circuit_name}')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Voltage (V)')
    ax.grid(True)
    ax.legend()
    
    plt.tight_layout()

    # Save plot
    output_dir = Path("test/res")
    output_dir.mkdir(parents=True, exist_ok=True)
    use_trap_as_ref = ref_df is trap_df
    ref_type = "trap_ref" if use_trap_as_ref else "ngspice_ref"
    plt.savefig(output_dir / f"{circuit_name}_{ref_type}_comparison.pdf", format="pdf", bbox_inches="tight")

    # Show plot
    plt.show()

    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Run and plot transient analysis results')
    parser.add_argument('circuit', help='Path to the circuit file')
    parser.add_argument('--M', type=int, default=1000, help='Number of Monte Carlo iterations')
    parser.add_argument('--N', type=int, default=100000, help='Number of time steps')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--use-trap-ref', action='store_true', 
                        help='Use trapezoidal method as reference instead of NgSpice')
    args = parser.parse_args()
    circuit_name = Path(args.circuit).stem

    print("Running simulations...")
    
    # Run trapezoidal simulation first
    print("1. Trapezoidal method...")
    trap_df, trap_time = run_trapezoidal(args.circuit, args.N, args.threads)
    print(trap_df.head())
    
    print("2. Monte Carlo method...")
    mc_df, mc_time = run_monte_carlo(args.circuit, args.M, args.N, args.threads)
    print(mc_df.head())
    
    # Run NgSpice only if not using trapezoidal as reference
    ref_df = None
    if not args.use_trap_ref:
        print("3. NgSpice reference...")
        ref_df = run_ngspice(args.circuit)
        print(ref_df.head())
    else:
        print("Using Trapezoidal method as reference (skipping NgSpice)")
        ref_df = trap_df  # Use trapezoidal as reference
    
    print("\nExecution Times:")
    print(f"Monte Carlo:  {mc_time:.4f} s")
    print(f"Trapezoidal: {trap_time:.4f} s")
    if mc_time and trap_time:
        speedup = trap_time / mc_time
        print(f"Speedup (Trap/MC): {speedup:.2f}x")
    
    print("\nCalculating errors...")
    errors = calculate_errors(mc_df, ref_df)
    
    print(f"\nError Metrics:")
    print(f"Max Relative Error:  {errors['max_error']:.2e}")
    print(f"Avg Relative Error:  {errors['avg_error']:.2e}")
    print(f"RMS Relative Error:  {errors['rms_error']:.2e}")
    
    print("\nGenerating plot...")
    plot_results(mc_df, trap_df, ref_df, errors, circuit_name)
    use_trap_as_ref = ref_df is trap_df
    ref_type = "trap_ref" if use_trap_as_ref else "ngspice_ref"
    print(f"Plot saved as test/res/{circuit_name}_{ref_type}_comparison.pdf")

if __name__ == '__main__':
    main()
