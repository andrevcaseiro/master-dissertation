import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the CSV file
main_df = pd.read_csv("res/time-expm-m-n.csv")

# Select fixed N
df = main_df[main_df["N"] == 8]

# Extract relevant columns
M = df["M"].values
TIME = df["TIME"].values

# Fit a quadratic curve: TIME = aM + b
a, b = np.polyfit(M, TIME, 1)

# Generate fitted values
M_fit = np.linspace(min(M), max(M), 100)  # Smooth x values for the curve
TIME_fit = a * M_fit + b  # Compute the quadratic function

# Calculate R2
TIME_pred = a * M + b  # Predicted values
SS_res = np.sum((TIME - TIME_pred) ** 2)  # Residual sum of squares
SS_tot = np.sum((TIME - np.mean(TIME)) ** 2)  # Total sum of squares
R2 = 1 - (SS_res / SS_tot)

# Create the plot
plt.figure(figsize=(8, 6))
plt.scatter(M, TIME, color="blue", label="Measured Data")  # Original data points
plt.plot(M_fit, TIME_fit, color="red", label="Linear Fit")  # Fitted curve

plt.xscale("log", base=2)

# Labels and title
plt.xlabel("M values (log scale)")
plt.ylabel("Execution Time (s)")
plt.title("Execution Time vs. M with Linear Fit")
plt.legend(title=f"t={a}M+{b}\nR^2: {R2}", loc="best")

# Save as a PDF for LaTeX
plt.savefig("res/execution_time_fit_m.pdf", format="pdf", bbox_inches="tight")


# Select fixed M
df = main_df[main_df["M"] == 8224]

# Extract relevant columns
N = df["N"].values
TIME = df["TIME"].values

# Fit a quadratic curve: TIME = aM + b
a, b = np.polyfit(N, TIME, 1)

# Generate fitted values
N_fit = np.linspace(min(N), max(N), 100)  # Smooth x values for the curve
TIME_fit = a * N_fit + b  # Compute the quadratic function

# Calculate R2
TIME_pred = a * N + b  # Predicted values
SS_res = np.sum((TIME - TIME_pred) ** 2)  # Residual sum of squares
SS_tot = np.sum((TIME - np.mean(TIME)) ** 2)  # Total sum of squares
R2 = 1 - (SS_res / SS_tot)

# Create the plot
plt.figure(figsize=(8, 6))
plt.scatter(N, TIME, color="blue", label="Measured Data")  # Original data points
plt.plot(N_fit, TIME_fit, color="red", label="Linear Fit")  # Fitted curve

plt.xscale("log", base=2)

# Labels and title
plt.xlabel("N values (log scale)")
plt.ylabel("Execution Time (s)")
plt.title("Execution Time vs. N with Linear Fit")
plt.legend(title=f"t={a}M+{b}\nR2: {R2}", loc="best")

# Save as a PDF for LaTeX
plt.savefig("res/execution_time_fit_n.pdf", format="pdf", bbox_inches="tight")
