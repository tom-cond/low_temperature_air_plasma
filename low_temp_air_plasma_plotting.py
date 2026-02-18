import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

# ── Configuration ──────────────────────────────────────────────────────────────
# Update this to the directory containing your .data files
DATA_DIR = "//wsl.localhost/Ubuntu-20.04/home/tom_cond/gdtk/src/kinetics"

CASE_LABELS = [
    "Case 0: E=12500 V/m",
    "Case 1: E=15000 V/m",
    "Case 2: E=17500 V/m",
    "Case 3: E=20000 V/m",
]

# Fixed parameters common to all cases
COMMON_PARAMS = "v_elec=2000 m/s, B=1 T, Qb=30 MW/m³"

# ── Load data ──────────────────────────────────────────────────────────────────
dfs = []
for v in range(4):
    fname = os.path.join(DATA_DIR, f"two_temperature_reacting_air_kinetics_test_results_v{v}.data")
    df = pd.read_csv(fname, delimiter=',', skipinitialspace=True)
    df.columns = df.columns.str.replace('"', '').str.strip()
    # Drop the duplicate t=0 row that appears at the start of each case
    df = df.drop_duplicates(subset='time (s)', keep='last').reset_index(drop=True)
    dfs.append(df)
    print(f"Loaded {fname}: {len(df)} rows")

# ── Temperature figure ─────────────────────────────────────────────────────────
fig_T, axes_T = plt.subplots(2, 2, figsize=(14, 9), sharex=False)
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
    f'Two-Temperature Air Kinetics — Temperature Evolution\n({COMMON_PARAMS})',
    fontsize=13, fontweight='bold', y=1.01
)
plt.tight_layout()
fig_T.savefig('kinetics_temperatures.png', dpi=150, bbox_inches='tight')
print("Saved kinetics_temperatures.png")

# ── Mass fraction figure ───────────────────────────────────────────────────────
NEUTRAL_SPECIES = {
    'N2 mass fraction':  ('N₂',  'steelblue',    '-'),
    'O2 mass fraction':  ('O₂',  'firebrick',     '-'),
    'N mass fraction':   ('N',   'darkorange',    '--'),
    'O mass fraction':   ('O',   'orchid',        '--'),
    'NO mass fraction':  ('NO',  'dimgrey',       ':'),
}

ION_SPECIES = {
    'e- mass fraction':  ('e⁻',   'black',        '-'),
    'N2+ mass fraction': ('N₂⁺',  'steelblue',    '--'),
    'O2+ mass fraction': ('O₂⁺',  'firebrick',    '--'),
    'N+ mass fraction':  ('N⁺',   'darkorange',   ':'),
    'O+ mass fraction':  ('O⁺',   'orchid',       ':'),
    'O2- mass fraction': ('O₂⁻',  'seagreen',     '-.'),
    'NO+ mass fraction': ('NO⁺',  'dimgrey',      '-.'),
}

def plot_species_panel(ax, df, species_dict, title, log_scale=True):
    for col, (name, color, ls) in species_dict.items():
        vals = df[col]
        if vals.max() > 0:
            ax.plot(df['time (s)'] * 1e6, vals,
                    label=name, color=color, linestyle=ls, linewidth=1.4)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.set_xlabel('Time (µs)')
    ax.set_ylabel('Mass Fraction')
    if log_scale:
        ax.set_yscale('log')
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    ax.legend(fontsize=8, ncol=2)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

# ── Neutral species subplots ───────────────────────────────────────────────────
fig_n, axes_n = plt.subplots(2, 2, figsize=(14, 9))
axes_n = axes_n.flatten()

for i, (ax, df, label) in enumerate(zip(axes_n, dfs, CASE_LABELS)):
    plot_species_panel(ax, df, NEUTRAL_SPECIES, label, log_scale=True)

fig_n.suptitle(
    f'Two-Temperature Air Kinetics — Neutral Species Mass Fractions\n({COMMON_PARAMS})',
    fontsize=13, fontweight='bold', y=1.01
)
plt.tight_layout()
fig_n.savefig('kinetics_neutral_species.png', dpi=150, bbox_inches='tight')
print("Saved kinetics_neutral_species.png")

# ── Charged species subplots ───────────────────────────────────────────────────
fig_i, axes_i = plt.subplots(2, 2, figsize=(14, 9))
axes_i = axes_i.flatten()

for i, (ax, df, label) in enumerate(zip(axes_i, dfs, CASE_LABELS)):
    plot_species_panel(ax, df, ION_SPECIES, label, log_scale=True)

fig_i.suptitle(
    f'Two-Temperature Air Kinetics — Charged Species Mass Fractions\n({COMMON_PARAMS})',
    fontsize=13, fontweight='bold', y=1.01
)
plt.tight_layout()
fig_i.savefig('kinetics_charged_species.png', dpi=150, bbox_inches='tight')
print("Saved kinetics_charged_species.png")

plt.show()