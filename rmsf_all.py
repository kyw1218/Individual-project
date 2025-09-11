import matplotlib
matplotlib.use('Agg')  # For headless plotting
import numpy as np
import matplotlib.pyplot as plt
import os

# Dictionary of systems with their replica files
rmsf_simulations = [
    ["rmsf-weca-rep1.xvg","rmsf-weca-rep2.xvg", "rmsf-weca-rep3.xvg"],
    ["rmsf-weca-udp-rep1.xvg", "rmsf-weca-udp-rep2.xvg", "rmsf-weca-udp-rep3.xvg"],
    ["rmsf_weca_c55_rep1.xvg", "rmsf_weca_c55_rep2.xvg", "rmsf_weca_c55_rep3.xvg"],
    ["rmsf-weca-glc-c55-rep1.xvg", "rmsf-weca-glc-c55-rep2.xvg", "rmsf-weca-glc-c55-rep3.xvg"]
]

simulation_labels = ["WECA", "WECA-GLC-MG", "WECA-C55", "WECA-GLC-C55"]
colors = ["blue", "green", "orange", "red"]

plt.figure(figsize=(10, 6))
x_values = None

for sim_idx, (replica_files, label, color) in enumerate(zip(rmsf_simulations, simulation_labels, colors)):
    sim_rmsf_replicas = []
    for fname in replica_files:
        try:
            x, y = np.loadtxt(fname, comments=("#", "@"), unpack=True)
            sim_rmsf_replicas.append(y)
            if x_values is None:
                x_values = x
        except Exception as e:
            print(f"Error reading {fname}: {e}")
    
    if sim_rmsf_replicas:
        sim_rmsf_replicas = np.array(sim_rmsf_replicas)
        sim_mean_rmsf = np.mean(sim_rmsf_replicas, axis=0)
        sim_std_rmsf = np.std(sim_rmsf_replicas, axis=0)
        
        # Optional: use standard error
        # sim_std_rmsf /= np.sqrt(len(sim_rmsf_replicas))

        color = colors[sim_idx]
        plt.plot(x_values, sim_mean_rmsf, label=label, color=color)
        plt.fill_between(x_values, sim_mean_rmsf - sim_std_rmsf, sim_mean_rmsf + sim_std_rmsf,
                         color=color, alpha=0.3)

        # Annotate max RMSF
        max_rmsf = np.max(sim_mean_rmsf)
        max_res = x_values[np.argmax(sim_mean_rmsf)]
        print(f"{label}: Max RMSF = {max_rmsf:.4f} nm at residue {int(max_res)}")

# Customize plot
plt.xlim(1, 350)
plt.ylim(0, 0.8)
plt.xlabel("Residue Index")
plt.ylabel("RMSF (nm)")
plt.title("RMSF for WECA Residues 1-350")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("rmsf_WECA_1-350.png", dpi=300)

