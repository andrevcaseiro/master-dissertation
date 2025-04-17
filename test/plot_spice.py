import subprocess
import matplotlib.pyplot as plt
import re
import pandas as pd
from plot_expected import df as python_df

""" result = subprocess.run(['make', 'clean'], capture_output=True, text=True)
print(result.stdout)

result = subprocess.run(['make', '-j'], capture_output=True, text=True)
print(result.stdout)

if result.returncode != 0:
    print("Make failed!")
    exit(1) """

# Run and parse NgSpice
result = subprocess.run(['ngspice', '-b', 'test/spice/simple.spice'], capture_output=True, text=True)
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

spice_df = pd.DataFrame(spice_data, columns=["index", "time", "v(v1)"])


# Run and parse this tool
result = subprocess.run(['./main', 'solve-spice', './test/spice/simple.spice', '100000', '1000', 'v1', '100', '-s'], capture_output=True, text=True)

data = []
for line in result.stdout.splitlines():
    parts = line.split()
    if len(parts) == 2:
        time = float(parts[0])
        voltage = float(parts[1])
        data.append((time, voltage))

df = pd.DataFrame(data, columns=["time", "v(v1)"])


# Plot both data frames

fig, ax = plt.subplots()

spice_df.plot(x='time', y='v(v1)', label='NgSpice', linewidth=3, ax=ax)
df.plot(x='time', y='v(v1)', label='Monte Carlo', linewidth=2, ax=ax)
python_df.plot(x='time', y='v(v1)', label='Scipy', linewidth=1, ax=ax)

ax.set_title('Comparison of Two Series')
ax.set_xlabel("Time (s)")
ax.set_ylabel("Voltage v(v1)")
ax.grid(True)
ax.legend()

plt.tight_layout()
plt.savefig("test/res/spice.pdf", format="pdf", bbox_inches="tight")
plt.show()
