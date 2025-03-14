import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load the CSV file
df_full = pd.read_csv("res/test-expm-m-n.csv")

df = df_full[df_full["N"] == 128]

# Compute the average error for each (M, N) pair
df_grouped = df.groupby(["M", "N"]).mean().reset_index()

M = df["M"]
ERROR = df["ERROR"]
grouped_M = df_grouped["M"]
avg_ERROR = df_grouped["ERROR"]

# Create the plot
plt.figure(figsize=(8, 6))
plt.scatter(M, ERROR, color="blue", label="Measured Data")  # Original data points
plt.plot(grouped_M, avg_ERROR, color="red", label="Average Error")

plt.xscale("log", base=2)

# Labels and title
plt.xlabel("M values (log scale)")
plt.ylabel("Relative error (%)")
plt.title("Relative error vs. M")
#plt.legend(title=f"t={a}M+{b}\nR2: {R2}", loc="best")

# Save as a PDF for LaTeX
plt.savefig("res/error_fit_m.pdf", format="pdf", bbox_inches="tight")

def f(x, a, b):
  return a*x**(-2) + b*x**(-1)

df = df_full[df_full["M"] == 263168]

# Compute the average error for each (M, N) pair
df_grouped = df.groupby(["M", "N"]).mean().reset_index()

N = df["N"]
ERROR = df["ERROR"]
grouped_N = df_grouped["N"].values.astype(float)
avg_ERROR = df_grouped["ERROR"]

popt, _ = curve_fit(f, grouped_N, avg_ERROR)
ERROR_pred = f(grouped_N, *popt)

# Create the plot
plt.figure(figsize=(8, 6))
plt.scatter(N, ERROR, color="blue", label="Measured Data")  # Original data points
plt.plot(grouped_N, avg_ERROR, color="lightblue", label="Average Error")
plt.plot(grouped_N, ERROR_pred, color="red", label="Fit")

plt.xscale("log", base=2)

# Labels and title
plt.xlabel("N values (log scale)")
plt.ylabel("Relative error (%)")
plt.title("Relative error vs. N")
plt.legend(loc="best")

# Save as a PDF for LaTeX
plt.savefig("res/error_fit_n.pdf", format="pdf", bbox_inches="tight")
