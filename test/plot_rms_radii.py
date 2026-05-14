#!/usr/bin/env python3
"""
Plot RMS radii comparison from gemstore bottomonium_1S0 output.
Dynamically extracts two sets of RMS radii from CLI output.
"""

import subprocess
import re
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 2:
    script_name = Path(sys.argv[0]).name
    print(f"Usage: python3 {script_name} <input.json>", file=sys.stderr)
    sys.exit(1)

input_json = sys.argv[1]


def parse_rms_table(table_text):
    values = []
    for line in table_text.splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        if not parts[0].isdigit():
            continue
        try:
            values.append(float(parts[2]))
        except ValueError:
            continue
    return values

# Run gemstore and capture output
gemstore_bin = './gemstore' if Path('./gemstore').is_file() else 'gemstore'
result = subprocess.run([gemstore_bin, '--compute', input_json], 
                       capture_output=True, text=True)

if result.returncode != 0:
    print(result.stdout, end="")
    print(result.stderr, end="", file=sys.stderr)
    sys.exit(result.returncode)

output = result.stdout + result.stderr

# Extract both MESON RESULTS SUMMARY tables
tables = re.findall(r'MESON RESULTS SUMMARY.*?(?=GLOBAL STATISTICS|Copyright|\Z)', output, re.DOTALL)

rms_radii_set1 = []
rms_radii_set2 = []

# Parse first table
if len(tables) > 0:
    rms_radii_set1 = parse_rms_table(tables[0])

# Parse second table
if len(tables) > 1:
    rms_radii_set2 = parse_rms_table(tables[1])

# If only one table found, set2 = set1
if len(tables) == 1:
    rms_radii_set2 = rms_radii_set1.copy()

if not rms_radii_set1 or not rms_radii_set2:
    print('Failed to parse RMS radii from gemstore output.', file=sys.stderr)
    sys.exit(1)

if len(rms_radii_set1) != len(rms_radii_set2):
    print(
        f'Parsed mismatched RMS lengths: set1={len(rms_radii_set1)}, set2={len(rms_radii_set2)}',
        file=sys.stderr,
    )
    sys.exit(1)

rms_radii_set1 = np.array(rms_radii_set1)
rms_radii_set2 = np.array(rms_radii_set2)
state_numbers = np.arange(1, len(rms_radii_set1) + 1)

# Create figure with subplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Both sets on the same plot
ax = axes[0, 0]
ax.plot(state_numbers, rms_radii_set1, 'o-', label='Set 1', linewidth=2, markersize=6)
ax.plot(state_numbers, rms_radii_set2, 's--', label='Set 2', linewidth=2, markersize=6)
ax.set_xlabel('State Number')
ax.set_ylabel('RMS Radius (fm)')
ax.set_title('RMS Radii Comparison: Two Sets Overlaid')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: Difference between sets
ax = axes[0, 1]
difference = rms_radii_set2 - rms_radii_set1
ax.bar(state_numbers, difference, color='steelblue', alpha=0.7)
ax.axhline(y=0, color='red', linestyle='--', linewidth=1)
ax.set_xlabel('State Number')
ax.set_ylabel('Difference (Set 2 - Set 1) (fm)')
ax.set_title('RMS Radii Difference Between Sets')
ax.grid(True, alpha=0.3, axis='y')

# Plot 3: Set 1 only
ax = axes[1, 0]
ax.plot(state_numbers, rms_radii_set1, 'o-', color='tab:blue', linewidth=2, markersize=8)
ax.fill_between(state_numbers, rms_radii_set1, alpha=0.3, color='tab:blue')
ax.set_xlabel('State Number')
ax.set_ylabel('RMS Radius (fm)')
ax.set_title('Set 1: RMS Radii')
ax.grid(True, alpha=0.3)

# Plot 4: Set 2 only
ax = axes[1, 1]
ax.plot(state_numbers, rms_radii_set2, 's-', color='tab:orange', linewidth=2, markersize=8)
ax.fill_between(state_numbers, rms_radii_set2, alpha=0.3, color='tab:orange')
ax.set_xlabel('State Number')
ax.set_ylabel('RMS Radius (fm)')
ax.set_title('Set 2: RMS Radii')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/serial/Code/Ongoing/gemstore/rms_radii_comparison.png', dpi=300, bbox_inches='tight')
print("Plot saved to rms_radii_comparison.png")

# Print statistics
print("\n" + "="*60)
print("RMS RADII STATISTICS")
print("="*60)
print(f"\nSet 1 Statistics:")
print(f"  Mean:       {np.mean(rms_radii_set1):.6f} fm")
print(f"  Std Dev:    {np.std(rms_radii_set1):.6f} fm")
print(f"  Min:        {np.min(rms_radii_set1):.6f} fm (State {np.argmin(rms_radii_set1) + 1})")
print(f"  Max:        {np.max(rms_radii_set1):.6f} fm (State {np.argmax(rms_radii_set1) + 1})")

print(f"\nSet 2 Statistics:")
print(f"  Mean:       {np.mean(rms_radii_set2):.6f} fm")
print(f"  Std Dev:    {np.std(rms_radii_set2):.6f} fm")
print(f"  Min:        {np.min(rms_radii_set2):.6f} fm (State {np.argmin(rms_radii_set2) + 1})")
print(f"  Max:        {np.max(rms_radii_set2):.6f} fm (State {np.argmax(rms_radii_set2) + 1})")

print(f"\nDifference Statistics:")
print(f"  Mean diff:  {np.mean(difference):.6f} fm")
print(f"  Max diff:   {np.max(np.abs(difference)):.6f} fm (State {np.argmax(np.abs(difference)) + 1})")
print(f"  States with differences > 0.1 fm:")

for i, (s1, s2, d) in enumerate(zip(rms_radii_set1, rms_radii_set2, difference), 1):
    if np.abs(d) > 0.1:
        print(f"    State {i:2d}: {s1:8.6f} → {s2:8.6f} (Δ = {d:+8.6f} fm)")
