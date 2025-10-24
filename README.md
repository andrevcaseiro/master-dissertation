# Master dissertation

Power grid simulation using a parallel Monte Carlo algorithm.

## CLI tool

The CLI tool supports the following commands:

- `solve` Solve $\mathbf{x}'=A\mathbf{x}+b$ with $\mathbf{x}(0)=\mathbf{x}_0$ at time $T$.
- `solve-sequence` Solve $\mathbf{x}'=A\mathbf{x}+b$ with $\mathbf{x}(0)=\mathbf{x}_0$ at times $0$ to $T$.
- `exp` Calculate $e^A$.
- `tran` Perform transient analysis on SPICE netlists using various ODE solvers.

 

## Prerequisites

This project links against external native libraries that must be available at build and at runtime:

- HighFM (libHighFM)
- SuperLU_MT (OpenMP variant)
- Intel oneAPI MKL and OpenMP runtime

Environment setup:

- Before building or running, source Intel oneAPI environment so MKL and the OpenMP runtime are on your include/library paths:
	- Bash: `source ~/intel/oneapi/setvars.sh` (or your install path, e.g. `/opt/intel/oneapi/setvars.sh`)
	- You need to do this in every new shell session prior to `make` or executing `./main`.
- Ensure the dynamic linker can find the HighFM and SuperLU_MT shared libraries if they are installed in non-system locations:
	- Example: `export LD_LIBRARY_PATH=/path/to/HighFM/lib:/path/to/superlu_mt/lib:$LD_LIBRARY_PATH`

Notes:

- The provided Makefile assumes MKL is installed under an Intel oneAPI directory and uses OpenMP threading.
- Depending on your distribution, you may also need `fmt` and `tbb` development packages installed from your package manager.

## Build and Run

After setting up the environment (see Prerequisites above), you can build and run the CLI:

```bash
make
```

Then execute the tool:

```
./main
```

## Project Structure

```
├── src/                    # Source code
│   ├── cli.cpp            # Command-line interface implementation
│   ├── *_ode_solver.cpp   # ODE solver implementations (trapezoidal, Monte Carlo, HighFM)
│   ├── spice/             # SPICE netlist parsing and analysis
│   │   ├── netlist.cpp    # Netlist parser with voltage source elimination
│   │   ├── mna.cpp        # Modified Nodal Analysis implementation
│   │   ├── dc_solver.cpp  # DC solver
│   │   └── ode.cpp        # ODE system generation from circuits
│   ├── matrix/            # Sparse matrix implementations
│   │   ├── coo_matrix.cpp # Coordinate format sparse matrix
│   │   └── csrd_matrix.cpp# Compressed sparse row sparse matrix with constant time diagonal acess
│   ├── highfm/            # HighFM algorithm components
│   ├── super_lu/          # SuperLU solver interface
│   ├── subcommands/       # CLI subcommands (solve, solve-sequence, exp, tran)
│   └── utils/             # Utility functions
│       ├── union_find.cpp # Union-Find data structure for graph algorithms
│       └── progress_bar.cpp# Progress visualization
├── external/              # External dependencies
│   ├── CLI11.hpp          # Command-line parsing library
│   ├── Eigen/             # Linear algebra library
│   └── indicators/        # Progress bar library
├── benchmark/             # Python benchmarks (scalability, time-vs-error, transient)
├── test/                  # Test suites
├── matlab/               # MATLAB reference implementations
└── res/                  # Results and output data
```

## Key Components

### SPICE Circuit Analysis

- **Netlist Parser**: Supports voltage sources, current sources, resistors, and capacitors
- **Voltage Source Elimination**: Advanced algorithm using Union-Find for efficient node merging
- **Modified Nodal Analysis (MNA)**: Converts circuits to ODE systems
- **DC Operating Point**: Finds steady-state solutions

### ODE Solvers

- **Monte Carlo Method**: Probabilistic approach for large-scale systems
- **Trapezoidal Method**: Implicit solver for stiff systems using different implementations

### Matrix Operations

- **Sparse Matrix Support**: COO and CSRD formats for memory efficiency
- **Linear Solvers**: Integration with SuperLU and Pardiso
- **Matrix Exponential**: Various algorithms for $e^{At}$ computation

## MATLAB

The `matlab` folder contains MATLAB scripts to generate some of the matrices and vectors in the `test` folder, as well as some scripts to generate expected results.

## Tests

The `test` folder contains multiple sets of input matrices/vectors and expected output calculated using MATLAB.
Furthermore, this folder includes bash and python scripts to generate data/plots analysing execution time/error with parameter variation.

## Benchmarking

The `benchmark/` folder contains Python scripts used to evaluate performance, scalability, and accuracy against the compiled CLI (`./main`):

- `scalability.py`: sweep threads/problem sizes and record speedup.
- `time_vs_error_analysis.py`: study the trade-off between runtime and error.
- `tran.py`: convenience runner for transient SPICE experiments.
- `main.py` and `utils.py`: shared entry points/utilities.

These scripts expect the `./main` binary to be built and available. Check `--help` on each script for available options and output locations.
