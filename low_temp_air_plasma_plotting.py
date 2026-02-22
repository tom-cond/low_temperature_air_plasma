"""
File which considers the 3 different test cases and creates plots for them
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

case = 'energy investigation'

if case == 'energy investigation':
    """
    Case which takes in the 4 data sets from the low_two_temperature_reacting_air_kinetics.d energy investigation unit test.
    Change DATA_DIR to whatever directory will contain said .data files (for windows systems using WSL, it will most
    likely be easier if you save (in the unit test main()) them to the kinetics file in the gdtk)
    """

    # Config
    # Update this to the directory containing your .data files
    DATA_DIR = "//wsl.localhost/Ubuntu-20.04/home/tom_cond/gdtk/src/kinetics"

    CASE_LABELS = [
        "Case 0: E=0 V/m",
        "Case 1: E=15000 V/m",
        "Case 2: E=17500 V/m",
        "Case 3: E=20000 V/m",
    ]

    # Fixed params common to all cases
    COMMON_PARAMS = "v_elec = 2000 m/s, MB = 1 T, Qb = 30 MW/m^3, T=Tve = 500K"

    # Load data
    dfs = []
    for v in range(4):
        fname = os.path.join(DATA_DIR, f"two_temperature_reacting_air_kinetics_test_results_v{v}.data")
        df = pd.read_csv(fname, delimiter=',', skipinitialspace=True)
        df.columns = df.columns.str.replace('"', '').str.strip()
        # Drop the duplicate t=0 row that appears at the start of each case
        df = df.drop_duplicates(subset='time (s)', keep='last').reset_index(drop=True)
        dfs.append(df)
        print(f"Loaded {fname}: {len(df)} rows")

    # Temp plot
    fig_T, axes_T = plt.subplots(2, 2, figsize=(14, 7.5), sharex=False) # We want a 2x2 set of plots
    axes_T = axes_T.flatten()

    for i, (ax, df, label) in enumerate(zip(axes_T, dfs, CASE_LABELS)):
        ax.plot(df['time (s)'] * 1e6, df['T (K)'],
                label='T (Trans-Rot)', color='tomato', linewidth=1.8)
        ax.plot(df['time (s)'] * 1e6, df['Tve (K)'],
                label='Tve (Vib-Elec)', color='royalblue', linestyle='--', linewidth=1.8)
        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.set_xlabel('Time (µs)')
        ax.set_ylabel('Temperature (K)')
        ax.grid(True, which='both', linestyle='--', alpha=0.5)
        ax.legend(fontsize=9)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    fig_T.suptitle(
        f'Two-Temperature Air Kinetics, Temperature Evolution\n({COMMON_PARAMS})',
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout(rect=[0, 0, 1, 0.975])
    fig_T.savefig('kinetics_temperatures.png', dpi=150, bbox_inches='tight')
    print("Saved kinetics_temperatures.png")

    # Mass fraction plot
    NEUTRAL_SPECIES = {
        'N2 mass fraction':  ('N2', 'steelblue' , '-' ),
        'O2 mass fraction':  ('O2', 'firebrick' , '-' ),
        'N mass fraction':   ('N' , 'darkorange', '--'),
        'O mass fraction':   ('O' , 'orchid'    , '--'),
        'NO mass fraction':  ('NO', 'dimgrey'   , ':' ),
    }

    ION_SPECIES = {
        'e- mass fraction':  ('e-' , 'black'     , '-' ),
        'N2+ mass fraction': ('N2+', 'steelblue' , '--'),
        'O2+ mass fraction': ('O2+', 'firebrick' , '--'),
        'N+ mass fraction':  ('N+' , 'darkorange', ':' ),
        'O+ mass fraction':  ('O+' , 'orchid'    , ':' ),
        'O2- mass fraction': ('O2-', 'seagreen'  , '-.'),
        'NO+ mass fraction': ('NO+', 'dimgrey'   , '-.'),
    }

    def plot_species_panel(ax, df, species_dict, title, log_scale=True):
        for col, (name, color, ls) in species_dict.items():
            vals = df[col]
            if vals.max() > 0: # Skip any species which do not have a mass fraction (useful for baseline cases)
                            # Might change to machine epsilon as there may be floating point mass fractions due to stiffness
                ax.plot(df['time (s)'] * 1e6, vals,
                        label=name, color=color, linestyle=ls, linewidth=1.4)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlabel('Time (µs)')
        ax.set_ylabel('Mass Fraction')
        if log_scale:
            ax.set_yscale('log')
        ax.grid(True, which='both', linestyle='--', alpha=0.5)
        ax.legend(fontsize=8, ncol=2)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator()) # Shoutout matplotlib

    # Neutral species subplots
    fig_n, axes_n = plt.subplots(2, 2, figsize=(14, 7.5))
    axes_n = axes_n.flatten()

    for i, (ax, df, label) in enumerate(zip(axes_n, dfs, CASE_LABELS)):
        plot_species_panel(ax, df, NEUTRAL_SPECIES, label, log_scale=True)

    fig_n.suptitle(
        f'Two-Temperature Air Kinetics, Neutral Species Mass Fractions\n({COMMON_PARAMS})',
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout(rect=[0, 0, 1, 0.975])
    fig_n.savefig('kinetics_neutral_species.png', dpi=150, bbox_inches='tight')
    print("Saved kinetics_neutral_species.png")

    # Charged species subplots
    fig_i, axes_i = plt.subplots(2, 2, figsize=(14, 7.5))
    axes_i = axes_i.flatten()

    for i, (ax, df, label) in enumerate(zip(axes_i, dfs, CASE_LABELS)):
        plot_species_panel(ax, df, ION_SPECIES, label, log_scale=True)

    fig_i.suptitle(
        f'Two-Temperature Air Kinetics, Charged Species Mass Fractions\n({COMMON_PARAMS})',
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout(rect=[0, 0, 1, 0.975])
    fig_i.savefig('kinetics_charged_species.png', dpi=150, bbox_inches='tight')
    print("Saved kinetics_charged_species.png")

    plt.show()

elif case == 'relaxation':
    """
    Case which takes in the relaxation data set from the low_two_temperature_reacting_air_kinetics.d relaxation unit test.
    This was made before I figured out how to use os, so you just need to specify the filename below instead of the directory
    """

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


elif case == 'sensitivity investigation':
    """
    File which takes in the 4 data sets from the low_two_temperature_reacting_air_kinetics.d sensitivity investigation unit test.
    Each file corresponds to a fixed initial gas temperature (T/Tve ratio), with Tve fixed at 300 K.
    Each subplot shows the electron mass fraction at t = 25 microseconds as a function of electric field strength.
    Change DATA_DIR to whatever directory will contain said .data files.
    """

    # Config
    DATA_DIR = "//wsl.localhost/Ubuntu-20.04/home/tom_cond/gdtk/src/kinetics"

    # Tve = 300 K fixed, T varies with T/Tve ratios of 1, 10, 20, 30
    CASE_LABELS = [
        "Case 0: T = 300 K,  Tve = 300 K  (T/Tve = 1)",
        "Case 1: T = 3000 K, Tve = 300 K  (T/Tve = 10)",
        "Case 2: T = 6000 K, Tve = 300 K  (T/Tve = 20)",
        "Case 3: T = 9000 K, Tve = 300 K  (T/Tve = 30)",
    ]

    COMMON_PARAMS = "v_elec = 2000 m/s, B = 1 T, Qb = 30 MW/m^3, t = 25 µs"

    # Load data
    dfs = []
    for v in range(4):
        fname = os.path.join(DATA_DIR, f"two_temperature_reacting_air_kinetics_test_results_v{v}.data")
        df = pd.read_csv(fname, delimiter=',', skipinitialspace=True)
        df.columns = df.columns.str.replace('"', '').str.strip()
        dfs.append(df)
        print(f"Loaded {fname}: {len(df)} rows")

    # Electron mass fraction vs electric field, one subplot per T/Tve ratio
    fig_s, axes_s = plt.subplots(2, 2, figsize=(14, 7.5), sharex=False)
    axes_s = axes_s.flatten()

    for i, (ax, df, label) in enumerate(zip(axes_s, dfs, CASE_LABELS)):
        ax.plot(df['E_field (V/m)'], df['e_mass_frac'],
                label='e-', color='#51247A', linewidth=1.5) #, marker='o', markersize=2)
        ax.plot(df['E_field (V/m)'], df['O2-_mass_fraction'],
            label='O2-', color='#2EA836', linewidth=1.5, linestyle='--') #, marker='o', markersize=0.5)
        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.set_xlabel('Electric Field (V/m)')
        ax.set_ylabel('Electron Mass Fraction')
        ax.set_yscale('log')
        ax.grid(True, which='both', linestyle='--', alpha=0.5)
        ax.legend(fontsize=9)
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    fig_s.suptitle(
        f'Negatively Charged Species Mass Fraction vs Electric Field Strength per T/Tve\n({COMMON_PARAMS})',
        fontsize=13, fontweight='bold'
    )
    plt.tight_layout(rect=[0, 0, 1, 0.975])
    fig_s.savefig('kinetics_sensitivity.png', dpi=150, bbox_inches='tight')
    print("Saved kinetics_sensitivity.png")

    plt.show()

else:
    raise Exception(f'No plotting script implemented for a {case} case')