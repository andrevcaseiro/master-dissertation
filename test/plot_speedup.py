import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the CSV file
df = pd.read_csv("res/time-expm-threads.csv")

# Extract relevant columns
THREADS = df["THREADS"].values
TIME = df["TIME"].values
SPEEDUP = TIME[0] / TIME

# Fit a quadratic curve: TIME = aM + b
a, b = np.polyfit(THREADS, SPEEDUP, 1)

# Generate fitted values
THREADS_fit = np.linspace(min(THREADS), max(THREADS), 100)  # Smooth x values for the curve
SPEEDUP_fit = a * THREADS_fit + b  # Compute the quadratic function

# Calculate R2
SPEEDUP_pred = a * THREADS + b  # Predicted values
SS_res = np.sum((SPEEDUP - SPEEDUP_pred) ** 2)  # Residual sum of squares
SS_tot = np.sum((SPEEDUP - np.mean(SPEEDUP)) ** 2)  # Total sum of squares
R2 = 1 - (SS_res / SS_tot)

# Create the plot
plt.figure(figsize=(8, 6))
plt.scatter(THREADS, SPEEDUP, color="blue", label="Measured Data")  # Original data points
plt.plot(THREADS_fit, SPEEDUP_fit, color="red", label="Linear Fit")  # Fitted curve

plt.xscale("log", base=2)

# Labels and title
plt.xlabel("No. of threads")
plt.ylabel("Execution Time (s)")
plt.title("Execution Time vs. No. of threads with Linear Fit")
plt.legend(title=f"SPEEDUP={a}THREADS+{b}\nR^2: {R2}", loc="best")

# Save as a PDF for LaTeX
plt.savefig("res/execution_time_fit_threads.pdf", format="pdf", bbox_inches="tight")
