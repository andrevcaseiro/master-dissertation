import subprocess
import matplotlib.pyplot as plt
import re
import pandas as pd
import numpy as np

filename = "test/spice/line5.spice"

# Run and parse NgSpice
result = subprocess.run(['ngspice', '-b', filename], capture_output=True, text=True)
if result.returncode != 0:
    print("ngspice failed!")
    exit(1)

data_line_pattern = re.compile(r'^\s*(\d+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)')

spice_data = []
for line in result.stdout.splitlines():
    match = data_line_pattern.match(line)
    if match:
        index = int(match.group(1))
        time = float(match.group(2))
        value = float(match.group(3))
        spice_data.append((index, time, value))

spice_df = pd.DataFrame(spice_data, columns=["index", "time", "voltage"])

# Run and parse Trap
cmd = ['./main', 'solve-trap', filename, '0', '-', '0', '--dc-solver', 'SparseLU']
result = subprocess.run(cmd, capture_output=True, text=True)

data = []
for line in result.stdout.splitlines():
    parts = line.split()
    if len(parts) == 2:
        time = float(parts[0])
        voltage = float(parts[1])
        data.append((time, voltage))

trap_df = pd.DataFrame(data, columns=["time", "voltage"])

# Run and parse Monte Carlo tool
cmd = ['./main', 'solve-spice', filename, '10000', '10000', '-', '0', '-s', '--output-freq', '100', '--dc-solver', 'SparseLU']
result = subprocess.run(cmd, capture_output=True, text=True)

data = []
for line in result.stdout.splitlines():
    parts = line.split()
    if len(parts) == 2:
        time = float(parts[0])
        voltage = float(parts[1])
        data.append((time, voltage))

df = pd.DataFrame(data, columns=["time", "voltage"])

# Calculate error
df["expected"] = np.interp(df['time'], spice_df['time'], spice_df['voltage'])

# Calculate error (e.g., absolute error, squared error, etc.)
df['abs_error'] = np.abs(df['voltage'] - df['expected'])
df['squared_error'] = (df['voltage'] - df['expected']) ** 2

# Optional: summary metrics
mae = df['abs_error'].mean()
mse = df['squared_error'].mean()
rmse = np.sqrt(mse)

print(f"MAE:  {mae:.4f}\nMSE:  {mse:.4f}\nRMSE: {rmse:.4f}")

# Plot both data frames

fig, ax = plt.subplots()

ax.plot(trap_df['time'], trap_df['voltage'], label='Trap', color="red", linewidth=2)
ax.plot(df['time'], df['expected'], label='NgSpice', color="orange", linewidth=2)
ax.plot(df['time'], df['voltage'], label='Monte Carlo', color="steelblue", linewidth=2)
#from plot_expected import df as pythondf
#ax.plot(pythondf['time'], pythondf['v(v5)'], label='Python', color="steelblue", linewidth=2)

#ax.fill_between(df['time'],df['expected'], df['voltage'], color='gray', alpha=0.3, label='Error band')
max_error_idx = df['abs_error'].idxmax()
max_time = df.loc[max_error_idx, 'time']
max_val = df.loc[max_error_idx, 'voltage']
max_error = df.loc[max_error_idx, 'abs_error']

# Scatter and annotation
dir = 1 if df.loc[max_error_idx, 'voltage'] > df.loc[max_error_idx, 'expected'] else -1
ax.annotate(f'Max error: {max_error:.2f}',
            xy=(max_time, max_val),
            xytext=(max_time, max_val + dir * 0.07 * (df['voltage'].max() - df['voltage'].min())),
            ha='center',
            va='bottom' if dir == 1 else 'top',
            arrowprops=dict(arrowstyle="->", color='red'),
            )


ax.set_title('Comparison of Two Series')
ax.set_xlabel("Time (s)")
ax.set_ylabel("Voltage (V)")
ax.grid(True)
ax.legend()

plt.tight_layout()
plt.savefig("test/res/spice.pdf", format="pdf", bbox_inches="tight")
plt.show()
