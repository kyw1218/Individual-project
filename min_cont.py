import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

def read_xvg(filename):
    """Reads a GROMACS .xvg file and returns time and data arrays."""
    time = []
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                time.append(float(parts[0]))
                data.append(float(parts[1]))
    return np.array(time), np.array(data)

mindist_simulations = [
    ["mindist_weca_rep1.xvg",     "mindist_weca_rep2.xvg",     "mindist_weca_rep3.xvg"],
    ["mindist_weca_glc_rep1.xvg", "mindist_weca_glc_rep2.xvg", "mindist_weca_glc_rep3.xvg"],
]

numcont_simulations = [
    ["numcont_weca_rep1.xvg",     "numcont_weca_rep2.xvg",     "numcont_weca_rep3.xvg"],
    ["numcont_weca_glc_rep1.xvg", "numcont_weca_glc_rep2.xvg", "numcont_weca_glc_rep3.xvg"],
]

labels = ["WECA", "WECA-GLC-MG"]

for mindist_group, numcont_group, label in zip(mindist_simulations, numcont_simulations, labels):
    # -- Mindist data
    mindist_data = [read_xvg(f)[1] for f in mindist_group]
    time = read_xvg(mindist_group[0])[0]  # Assume same time in all repeats
    mindist_array = np.array(mindist_data)
    mindist_mean = np.mean(mindist_array, axis=0)
    mindist_sem = np.std(mindist_array, axis=0) / np.sqrt(len(mindist_array))

    # -- Numcont data
    numcont_data = [read_xvg(f)[1] for f in numcont_group]
    numcont_array = np.array(numcont_data)
    numcont_mean = np.mean(numcont_array, axis=0)
    numcont_sem = np.std(numcont_array, axis=0) / np.sqrt(len(numcont_array))

    # ---- Plotting dual Y-axis ----
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Left Y-axis: mindist
    ax1.set_xlabel("Time (ps)")
    ax1.set_ylabel("Average Protein-Membrane\nDistance (nm)", color='green')
    ax1.plot(time, mindist_mean, color='green')
    ax1.fill_between(time, mindist_mean - mindist_sem, mindist_mean + mindist_sem, color='green', alpha=0.3)
    ax1.tick_params(axis='y', labelcolor='green')

    # Right Y-axis: numcont
    ax2 = ax1.twinx()
    ax2.set_ylabel("Average PC Contacts", color='purple')
    ax2.plot(time, numcont_mean, color='purple')
    ax2.fill_between(time, numcont_mean - numcont_sem, numcont_mean + numcont_sem, color='purple', alpha=0.3)
    ax2.tick_params(axis='y', labelcolor='purple')

    # Title and layout
    plt.title(f"{label}: Distance vs. Contacts")
    plt.tight_layout()

    # Save figure
    plt.savefig(f"{label}_mindist_numcont_dual_axis.png", dpi=300)

