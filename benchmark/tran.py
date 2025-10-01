#!/usr/bin/env python3
import os
import subprocess
import argparse
import pandas as pd
from datetime import datetime


def run_and_save(command, filepath, env_vars=None):
    """
    Run a command and save its output to a file.

    Args:
        command (str): The command to run
        filepath (str): Full path to save the output
        env_vars (dict, optional): Environment variables to set. Defaults to None.

    Returns:
        str: Output of the command
    """
    print("> " + command)
    try:
        with open(filepath, "r") as f:
            # Get file creation/modification time
            file_stat = os.stat(filepath)
            cached_time = datetime.fromtimestamp(file_stat.st_mtime)
            time_diff = datetime.now() - cached_time
            
            # Format time difference in a human-readable way
            time_units = [
                ("day", time_diff.days),
                ("hour", time_diff.seconds // 3600),
                ("minute", time_diff.seconds // 60),
                ("second", time_diff.seconds)
            ]
            
            for unit, value in time_units:
                if value > 0:
                    time_ago = f"{value} {unit}{'s' if value != 1 else ''} ago"
                    break
            
            print(f"> Using cached result from {cached_time.strftime('%Y-%m-%d %H:%M:%S')} ({time_ago})")
            return f.read()
    except FileNotFoundError:
        # Setup environment variables if provided
        env = os.environ.copy()
        if env_vars:
            env.update(env_vars)

        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True, env=env)
        output = result.stdout
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, "w") as f:
            f.write(output)
        return output


def run_monte_carlo(netlist_path, final_time, num_steps, num_samples=1000, seed=-1, print_step=1000, num_threads=None, verbose=False):
    """
    Run transient analysis with Monte Carlo ODE solver.
    
    Args:
        netlist_path (str): Path to the SPICE netlist file
        final_time (float): Final simulation time
        num_steps (int): Number of time steps
        num_samples (int, optional): Number of Monte Carlo samples. Defaults to 1000.
        seed (int, optional): Random seed (-1 for random). Defaults to -1.
        print_step (int, optional): Print every N-th step. Defaults to 1000.
        num_threads (int, optional): Number of OpenMP threads to use. Defaults to None (use system default).
        verbose (bool, optional): Print detailed information. Defaults to False.
    
    Returns:
        tuple: (DataFrame with time and voltage columns, execution time in seconds, DC analysis time in seconds)
    """
    # Set up environment variables for OpenMP
    env_vars = {}
    threads_str = ""
    if num_threads is not None:
        env_vars["OMP_NUM_THREADS"] = str(num_threads)
        threads_str = f"_threads{num_threads}"
    
    cmd = f"./main tran {netlist_path} --solver monte-carlo --method pardiso -t {final_time} -N {num_steps} -M {num_samples} -s {seed} -p {print_step}"
    
    if verbose:
        cmd += " --verbose"
    
    # Remove file extension from basename for output path
    netlist_basename = os.path.splitext(os.path.basename(netlist_path))[0]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, f"{netlist_basename}/monte_carlo/t{final_time}_N{num_steps}_M{num_samples}_s{seed}_p{print_step}{threads_str}.out")
    
    output = run_and_save(cmd, output_path, env_vars)
    
    # Parse output to get data and execution time
    data = []
    exec_time = None
    dc_time = None

    for line in output.splitlines():
        if "ODE solver time:" in line:
            exec_time = float(line.split(":")[1].split()[0])
            continue
        if "DC analysis time:" in line:
            dc_time = float(line.split(":")[1].split()[0])
            continue
        try:
            parts = line.strip().split()
            if len(parts) == 2:
                time = float(parts[0])
                voltage = float(parts[1])
                data.append((time, voltage))
        except (ValueError, IndexError):
            continue
    
    return pd.DataFrame(data, columns=["time", "voltage"]), exec_time, dc_time


def run_trapezoidal(netlist_path, final_time, num_steps, method="pardiso", num_threads=None, verbose=False):
    """
    Run transient analysis with Trapezoidal ODE solver.
    
    Args:
        netlist_path (str): Path to the SPICE netlist file
        final_time (float): Final simulation time
        num_steps (int): Number of time steps
        method (str, optional): Linear solver method (lu, cg, slu). Defaults to "slu".
        num_threads (int, optional): Number of OpenMP threads to use. Defaults to None (use system default).
        verbose (bool, optional): Print detailed information. Defaults to False.
    
    Returns:
        tuple: (DataFrame with time and voltage columns, execution time in seconds, DC analysis time in seconds)
    """
    # Set up environment variables for OpenMP
    env_vars = {}
    threads_str = ""
    if num_threads is not None:
        env_vars["OMP_NUM_THREADS"] = str(num_threads)
        threads_str = f"_threads{num_threads}"
    
    cmd = f"./main tran {netlist_path} --solver trapezoidal --method {method} -t {final_time} -N {num_steps}"
    
    if verbose:
        cmd += " --verbose"
    
    # Remove file extension from basename for output path
    netlist_basename = os.path.splitext(os.path.basename(netlist_path))[0]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, f"{netlist_basename}/trapezoidal/t{final_time}_N{num_steps}_{method}{threads_str}.out")

    output = run_and_save(cmd, output_path, env_vars)
    
    # Parse output to get data and execution time
    data = []
    exec_time = None
    dc_time = None
    
    for line in output.splitlines():
        if "ODE solver time:" in line:
            exec_time = float(line.split(":")[1].split()[0])
            continue
        if "DC analysis time:" in line:
            dc_time = float(line.split(":")[1].split()[0])
            continue
        try:
            parts = line.strip().split()
            if len(parts) == 2:
                time = float(parts[0])
                voltage = float(parts[1])
                data.append((time, voltage))
        except (ValueError, IndexError):
            continue
    
    return pd.DataFrame(data, columns=["time", "voltage"]), exec_time, dc_time


def run_ngspice(netlist_path, verbose=False):
    """
    Run SPICE simulation using ngspice.
    
    Args:
        netlist_path (str): Path to the SPICE netlist file
        verbose (bool, optional): Print detailed information. Defaults to False.
    
    Returns:
        tuple: (DataFrame with time and voltage columns, execution time in seconds)
    """
    # Remove file extension from basename for output path
    netlist_basename = os.path.splitext(os.path.basename(netlist_path))[0]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, f"{netlist_basename}/ngspice/res.out")
    
    cmd = f"ngspice -b {netlist_path}"
    if verbose:
        cmd += " -v"
    
    # Run ngspice
    output = run_and_save(cmd, output_path)
    
    # Parse output to get data
    data = []
    # Skip until we find the data section (starts with Index header)
    found_data = False
    for line in output.splitlines():
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
    
    # ngspice doesn't provide execution time directly, so we return None
    return pd.DataFrame(data, columns=["time", "voltage"]), None, None


def main():
    """
    Main function to run transient analysis with different solvers when script is executed directly.
    """
    parser = argparse.ArgumentParser(description="Run transient analysis on SPICE netlists with different ODE solvers")
    parser.add_argument("netlist_path", help="Path to the SPICE netlist file")
    parser.add_argument("-s", "--solver", choices=["monte-carlo", "trapezoidal", "ngspice"], default="monte-carlo",
                        help="Solver method to use (default: monte-carlo)")
    parser.add_argument("-t", "--time", type=float, required=True,
                        help="Final simulation time")
    parser.add_argument("-N", "--steps", type=int, required=True,
                        help="Number of time steps")
    parser.add_argument("-M", "--samples", type=int, default=1000,
                        help="Number of Monte Carlo samples (only for monte-carlo solver)")
    parser.add_argument("--seed", type=int, default=-1,
                        help="Random seed (-1 for random, only for monte-carlo solver)")
    parser.add_argument("-p", "--print-step", type=int, default=1000,
                        help="Print every N-th step (only for monte-carlo solver)")
    parser.add_argument("--method", choices=["lu", "cg", "slu"], default="slu",
                        help="Linear solver method for trapezoidal integration (default: slu)")
    parser.add_argument("--threads", type=int, 
                        help="Number of OpenMP threads to use")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print detailed information")
    
    args = parser.parse_args()
    
    if args.solver == "monte-carlo":
        output = run_monte_carlo(
            netlist_path=args.netlist_path,
            final_time=args.time,
            num_steps=args.steps,
            num_samples=args.samples,
            seed=args.seed,
            print_step=args.print_step,
            num_threads=args.threads,
            verbose=args.verbose
        )
    elif args.solver == "trapezoidal":
        output = run_trapezoidal(
            netlist_path=args.netlist_path,
            final_time=args.time,
            num_steps=args.steps,
            method=args.method,
            num_threads=args.threads,
            verbose=args.verbose
        )
    else:  # ngspice
        output = run_ngspice(
            netlist_path=args.netlist_path,
            verbose=args.verbose
        )
    
    print(output)


if __name__ == "__main__":
    main()
