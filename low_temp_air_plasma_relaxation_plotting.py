"""
File which takes in the relaxation data set from the low_two_temperature_reacting_air_kinetics.d relaxation unit test.
This was made before I figured out how to use os, so you just need to specify the filename below instead of the directory
"""

import pandas as pd
import matplotlib.pyplot as plt

# Read the data, delimiter=',' tells pandas to split by commas
# skipinitialspace=True is needed because D code does format("%e, %e") causes spaces to appear after commas
filename = "//wsl.localhost/Ubuntu-20.04/home/tom_cond/gdtk/src/kinetics/two_temperature_reacting_air_kinetics_test_results.data"

try:
    df = pd.read_csv(filename, delimiter=',', skipinitialspace=True)
    
    # Strip quotes from column names if they were included in the header
    df.columns = df.columns.str.replace('"', '').str.strip()
    
    print("Successfully read data. Columns found:")
    print(df.columns)

except FileNotFoundError:
    print(f"Error: Could not find file '{filename}'.")
    exit()

# We use subplots because Temperature and Mass Fraction have vastly different scales
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Subplot 1: Temperatures
# Plot T (Transrotational) and Tve (Vibroelectronic)
ax1.plot(df['time (s)'], df['T (K)'], label='T (Trans-Rot)', color='red', linewidth=2)
ax1.plot(df['time (s)'], df['Tve (K)'], label='Tve (Vib-Elec)', color='blue', linestyle='--', linewidth=2)

ax1.set_ylabel('Temperature (K)')
ax1.set_title('Relaxation of Two-Temperature Air')
ax1.grid(True, which='both', linestyle='--', alpha=0.7)
ax1.legend()

# Subplot 2: Mass Fractions
# Plot Species (N, O, NO, N2, O2)
ax2.plot(df['time (s)'], df['N mass fraction'], label='N', color='green')
ax2.plot(df['time (s)'], df['O mass fraction'], label='O', color='pink')
ax2.plot(df['time (s)'], df['NO mass fraction'], label='NO', color='black')
ax2.plot(df['time (s)'], df['N2 mass fraction'], label='N2', color='blue')
ax2.plot(df['time (s)'], df['O2 mass fraction'], label='O2', color='red')


ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Mass Fraction')
ax2.set_yscale('log') # Log scale is better for species traces
ax2.grid(True, which='both', linestyle='--', alpha=0.7)
ax2.legend()

# Final Layout Adjustments
plt.tight_layout()

# Save the plot
plt.savefig('relaxation_results.png')
print("Plot saved as 'relaxation_results.png'")

# Show the plot
plt.show()