import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd

def pulse(t, v1, v2, td, tr, pw, tf, per):
    if t < td:
      return v1
    cycle_time = (t - td) % per

    if cycle_time < tr:
       return v1 + (v2 - v1) * (cycle_time / tr)

    if cycle_time < tr + pw:
       return v2

    if cycle_time < tr + pw + tf:
       return v2 - (v2 - v1) * ((cycle_time - tr - pw) / tf)

    return v1

# Matrix A
A = np.array([[-0.250, 0.250, 0.000],
              [ 0.250,-0.500, 0.250],
              [ 0.000, 0.250,-0.500]])

# Time-dependent b(t)
def b(t):
    return np.array([pulse(t, 0.5, 2.5, 0.1, 0.2, 0.2, 0.2, 5), 0, 0])

# System of equations: dx/dt = A * x + b(t)
def system(t, x):
    return A @ x + b(t)

# Initial condition
x0 = [6.000, 4.000, 2.000]

# Time span
t_span = (0, 9.5)
t_eval = np.linspace(*t_span, 200)

# Solve the system
solution = solve_ivp(system, t_span, x0, t_eval=t_eval, max_step=0.1)

# Extract solutions
t = solution.t

t = solution.t              # Time points
y = solution.y              # Solution values (shape: [n_vars, n_times])

# Transpose y to align each row with a time point
df = pd.DataFrame(y.T, columns=[f'v(v{i+1})' for i in range(y.shape[0])])
df.insert(0, 'time', t)  # Insert time as the first column


if __name__ == "__main__":
    print(solution.y[0][-1])

    # Plot
    plt.figure(figsize=(10, 5))

    plt.subplot(1, 2, 1)
    for i in range(len(x0)):
        plt.plot(t, solution.y[i], label=f'x{i}(t)')
    plt.xlabel('Time t')
    plt.ylabel('Values')
    plt.title('x(t) and y(t)')
    plt.legend()
    plt.grid(True)

    t_test = np.linspace(0, 10, 500)
    pulse_vals = [pulse(ti, 0, 5, 0.1, 0.2, 0.2, 0.2, 1) for ti in t_test]


    plt.subplot(1,2,2)
    plt.plot(t_test, pulse_vals)
    plt.title("Pulse Function")
    plt.xlabel("Time (t)")
    plt.ylabel("Pulse Value")
    plt.grid(True)

    plt.tight_layout()
    plt.show()
