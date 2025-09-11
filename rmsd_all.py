import matplotlib
matplotlib.use('Agg')  # No GUI
import numpy as np
import matplotlib.pyplot as plt

# Define RMSD replica files grouped by simulation
# Each sublist corresponds to one simulation's 3 replicas
rmsd_simulations = [
    ["rmsd-WECA-rep1.xvg","rmsd-WECA-rep2.xvg", "rmsd-WECA-rep3.xvg"],
    ["rmsd-weca-udp-rep1.xvg", "rmsd-weca-udp-rep2.xvg", "rmsd-weca-udp-rep3.xvg"],
    ["rmsd_weca_c55_rep1.xvg", "rmsd_weca_c55_rep2.xvg", "rmsd_weca_c55_rep3.xvg"],
    ["rmsd-weca-glc-c55-rep1.xvg", "rmsd-weca-glc-c55-rep2.xvg", "rmsd-weca-glc-c55-rep3.xvg"]
]

simulation_labels = ["WECA", "WECA-GLC-MG", "WECA-C55", "WECA-GLC-C55"]
colors = ["blue", "green", "orange", "red"]

# Start plot
plt.figure(figsize=(9, 6))

# Loop over each simulation group
for sim_idx, (replica_files, label, color) in enumerate(zip(rmsd_simulations, simulation_labels, colors)):
    rmsd_data = []
    x_values = None

    for fname in replica_files:
        try:
            x, y = np.loadtxt(fname, comments=("#", "@"), unpack=True)
            if x_values is None:
                x_values = x
            rmsd_data.append(y)
        except Exception as e:
            print(f"Error reading {fname}: {e}")

    # Compute mean and std across replicas
    rmsd_array = np.array(rmsd_data)
    mean_rmsd = np.mean(rmsd_array, axis=0)
    std_rmsd = np.std(rmsd_array, axis=0)

    # Plot with error band
    plt.plot(x_values, mean_rmsd, label=label, color=color)
    plt.fill_between(x_values, mean_rmsd - std_rmsd, mean_rmsd + std_rmsd,
                     color=color, alpha=0.3)

    # Optionally print max RMSD for this simulation
    max_rmsd = np.max(mean_rmsd)
    max_rmsd_time = x_values[np.argmax(mean_rmsd)]
    print(f"{label} - Max RMSD: {max_rmsd:.4f} at {max_rmsd_time:.1f} ps")

# Final formatting
plt.xlabel("Time (ps)")
plt.ylabel("RMSD (nm)")
plt.title("RMSD Comparison for different WECA Simulations")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("rmsd_WECA.png", dpi=300)
