import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

BASE = "/Users/bhatti.39/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Current_Research/Research Caspase-4/New graphs/XVGs_extracted/XVGs"

DRUGS = ["DB00519", "DB01068", "DB05316", "DB06202", "DB08882", "Donepezil"]
REPLICATES = ["R1", "R2", "R3"]

# Colors matching the reference figure style
REP_COLORS = {"R1": "#E8512A", "R2": "#2878B8", "R3": "#2CA02C"}  # red, blue, green
AVG_COLOR = "#555555"


def read_xvg(path):
    """Parse a GROMACS .xvg file, returning (time_array, value_array)."""
    times, vals = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(("#", "@")) or not line:
                continue
            parts = line.split()
            times.append(float(parts[0]))
            vals.append(float(parts[1]))
    return np.array(times), np.array(vals)


fig, axes = plt.subplots(2, 3, figsize=(15, 8))
axes = axes.flatten()

for idx, drug in enumerate(DRUGS):
    ax = axes[idx]
    all_vals = []

    for rep in REPLICATES:
        fpath = os.path.join(BASE, drug, rep, "rmsd_backbone.xvg")
        if not os.path.exists(fpath):
            print(f"Missing: {fpath}")
            continue
        t, v = read_xvg(fpath)
        avg_rep = np.mean(v)
        all_vals.append(v)
        label = f"{rep} (avg: {avg_rep:.2f})"
        ax.plot(t, v, color=REP_COLORS[rep], linewidth=0.7, alpha=0.85, label=label)

    # Overall average dashed line
    if all_vals:
        # Align lengths to shortest replicate
        min_len = min(len(v) for v in all_vals)
        stacked = np.array([v[:min_len] for v in all_vals])
        overall_avg = np.mean(stacked)
        ax.axhline(overall_avg, color=AVG_COLOR, linestyle="--", linewidth=1.2,
                   label=f"Overall Avg: {overall_avg:.2f}")

    ax.set_title(drug, fontsize=11, fontweight="bold")
    ax.set_xlabel("Time (ns)", fontsize=9)
    ax.set_ylabel("RMSD (nm)", fontsize=9)
    ax.set_ylim(0, 0.6)
    ax.set_xlim(left=0)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.7)
    ax.grid(False)

fig.suptitle("Backbone RMSD — MD Simulations (Triplicates)", fontsize=13, fontweight="bold", y=1.01)
plt.tight_layout()

out_path = "/Users/bhatti.39/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Current_Research/Research Caspase-4/New graphs/RMSD_triplicates.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved: {out_path}")
