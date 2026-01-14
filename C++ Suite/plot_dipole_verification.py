import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_radial():
    df = pd.read_csv("dipole_radial_verification.csv")
    
    # Filter for D values
    D_vals = df['D'].unique()
    
    fig, axes = plt.subplots(len(D_vals), 1, figsize=(8, 12), sharex=True)
    if len(D_vals) == 1: axes = [axes]
    
    for i, D in enumerate(D_vals):
        ax = axes[i]
        subset = df[df['D'] == D]
        
        # Plot N=0
        n0 = subset[subset['N'] == 0]
        ax.plot(n0['r'], n0['Re_Radial'], label=f"N=0 (l_eff={n0['L_eff'].iloc[0]:.2f})", linestyle='-')
        ax.plot(n0['r'], n0['Re_Ref_Bessel'], label="Ref l=0", linestyle='--')
        
        # Plot N=1
        n1 = subset[subset['N'] == 1]
        ax.plot(n1['r'], n1['Re_Radial'], label=f"N=1 (l_eff={n1['L_eff'].iloc[0]:.2f})", linestyle='-')
        ax.plot(n1['r'], n1['Re_Ref_Bessel'], label="Ref l=1", linestyle='--')
        
        ax.set_title(f"Radial Functions (D={D})")
        ax.legend()
        ax.set_ylabel("Amplitude")
        
    axes[-1].set_xlabel("r (au)")
    plt.tight_layout()
    plt.savefig("dipole_radial_plot.png")
    print("Saved dipole_radial_plot.png")

def plot_angular():
    df = pd.read_csv("dipole_angular_verification.csv")
    D_vals = df['D'].unique()
    
    fig, axes = plt.subplots(len(D_vals), 1, figsize=(8, 12), sharex=True)
    if len(D_vals) == 1: axes = [axes]
    
    for i, D in enumerate(D_vals):
        ax = axes[i]
        subset = df[df['D'] == D]
        
        # N=0
        n0 = subset[subset['N'] == 0]
        ax.plot(n0['theta'], n0['Re_Omega'], label="N=0 Omega", linestyle='-')
        ax.plot(n0['theta'], n0['Re_Ref_Ylm'], label="Ref Y_00", linestyle='--')
        
        # N=1
        n1 = subset[subset['N'] == 1]
        ax.plot(n1['theta'], n1['Re_Omega'], label="N=1 Omega", linestyle='-')
        ax.plot(n1['theta'], n1['Re_Ref_Ylm'], label="Ref Y_10", linestyle='--')
        
        ax.set_title(f"Angular Functions (D={D})")
        ax.legend()
        ax.set_ylabel("Amplitude")

    axes[-1].set_xlabel("Theta (rad)")
    plt.tight_layout()
    plt.savefig("dipole_angular_plot.png")
    print("Saved dipole_angular_plot.png")

if __name__ == "__main__":
    plot_radial()
    plot_angular()
