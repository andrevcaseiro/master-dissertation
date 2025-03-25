# Master dissertation
Power grid simulation using a parallel Monte Carlo algorithm.

## CLI tool

The CLI tool supports the following commands:
- `solve` Solve $\mathbf{x}'=A\mathbf{x}+b$ with $\mathbf{x}(0)=\mathbf{x}_0$ at time $T$.
- `solve-sequence` Solve $\mathbf{x}'=A\mathbf{x}+b$ with $\mathbf{x}(0)=\mathbf{x}_0$ at times $0$ to $T$.
- `exp` Calculate $e^A$.

To install dependencies and compile the tool, run:
```bash
make
```

To execute the tool, run:
```
./main
```

## Matlab

The `matlab` folder contains Matlab scripts to generate some of the matrices and vectors in the `test` folder, as well as some scripts to generate expected results.

## Tests

The `test` folder contains multiple sets of input matrices/vectors and expected output calculated using Matlab.
Furthermore, this folder includes bash and python scripts to generate data/plots analysing execution time/error with parameter variation.
