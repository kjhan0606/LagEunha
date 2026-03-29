import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# # # #  # # 
try:
    df = pd.read_csv('gce_data.csv')
except FileNotFoundError:
    print("Error: 'gce_data.csv' not found. Please run the C simulation first.")
    exit()

# # # #  # # 
plt.figure(figsize=(10, 8))
plt.style.use('default')  # # #  'seaborn-v0_8-whitegrid'

# Scatter plot with color mapping for Time
sc = plt.scatter(df['Fe_H'], df['O_Fe'], c=df['Time_Myr'], 
                 cmap='viridis', s=50, edgecolor='k', alpha=0.8, label='Simulation Track')

# # # #  # # # 
plt.colorbar(sc, label='Time (Myr)')
plt.xlabel('[Fe/H]', fontsize=14)
plt.ylabel('[O/Fe]', fontsize=14)
plt.title('Galactic Chemical Evolution: [O/Fe] vs [Fe/H]', fontsize=16)

# SN Ia "Knee" # #  # #  (# # )
# # #  100 Myr ~ 1 Gyr # # # #  # # #  # # # #  # 
plt.axhline(y=0.0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
plt.axvline(x=0.0, color='gray', linestyle='--', alpha=0.5, linewidth=1)

# # #  # # 
plt.text(-2.5, 0.6, 'SN II Dominated\n(Alpha Enhanced)', fontsize=12, color='blue', ha='left')
plt.text(-0.5, 0.1, 'SN Ia Onset\n(Iron Production)', fontsize=12, color='red', ha='left')

# Grid #  # #  # # 
plt.grid(True, linestyle=':', alpha=0.6)
plt.xlim(-4.0, 0.5)  # # # # #  # #  # # 
plt.ylim(-0.2, 1.0)  # # # # #  # #  # # 

# # # # #  # #  # #  # # 
indices = [int(len(df)*0.1), int(len(df)*0.5), int(len(df)*0.9)]
for i in indices:
    if i+5 < len(df):
        plt.annotate('', xy=(df['Fe_H'].iloc[i+5], df['O_Fe'].iloc[i+5]), 
                     xytext=(df['Fe_H'].iloc[i], df['O_Fe'].iloc[i]),
                     arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

plt.tight_layout()
plt.savefig('ofe_vs_feh_plot.png', dpi=300)
plt.show()

print("Plot saved as 'ofe_vs_feh_plot.png'")
