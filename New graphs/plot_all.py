import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

BASE = "/Users/bhatti.39/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Current_Research/Research Caspase-4/New graphs/XVGs_extracted/XVGs"
OUT  = "/Users/bhatti.39/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Current_Research/Research Caspase-4/New graphs"

DRUGS         = ["DB00519", "DB01068", "DB05316", "DB06202", "DB08882"]
REPLICATES    = ["R1", "R2", "R3"]
PANEL_LETTERS = ["A", "B", "C", "D", "E"]

# Distinct, well-separated colors (purple/red/blue/green don't visually clash)
DONEPEZIL_COLOR = "#9467BD"   # purple — reference compound
REP_COLORS      = {"R1": "#5B9BD5", "R2": "#ED7D31", "R3": "#70AD47"}  # blue, orange, green
AVG_COLOR       = "#555555"


def open_box(ax):
    """Drop the top and right spines for an open-box look."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_xvg(path, col=1, ps_to_ns=False):
    times, vals = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(("#", "@")) or not line:
                continue
            parts = line.split()
            try:
                t = float(parts[0])
                v = float(parts[col])
            except (IndexError, ValueError):
                continue
            if ps_to_ns:
                t /= 1000.0
            times.append(t)
            vals.append(v)
    return np.array(times), np.array(vals)


def rolling_mean(arr, window=200):
    kernel = np.ones(window) / window
    return np.convolve(arr, kernel, mode="same")


def load_donepezil_ref(file_key, col=1, ps_to_ns=False, smooth=False,
                        smooth_window=200, is_rmsf=False):
    """
    Return (x, mean_y, avg_scalar) for Donepezil across R1/R2/R3.
    For time-series: x = time array (shortest replicate), mean_y = mean across reps.
    For RMSF: x = residue array (shortest), mean_y = mean across reps.
    """
    all_x, all_v = [], []
    for rep in REPLICATES:
        fpath = os.path.join(BASE, "Donepezil", rep, file_key)
        if not os.path.exists(fpath):
            continue
        x, v = read_xvg(fpath, col=col, ps_to_ns=ps_to_ns)
        all_x.append(x)
        all_v.append(v)

    if not all_v:
        return None, None, None

    min_len   = min(len(v) for v in all_v)
    x_ref     = all_x[0][:min_len]
    stacked   = np.array([v[:min_len] for v in all_v])
    mean_y    = np.mean(stacked, axis=0)
    avg_scalar = np.mean(mean_y)

    if smooth:
        mean_y = rolling_mean(mean_y, smooth_window)

    return x_ref, mean_y, avg_scalar


def make_figure(title, ylabel, filename,
                file_key, col=1, ps_to_ns=False,
                smooth=False, smooth_window=200,
                is_rmsf=False):

    # Pre-load Donepezil reference
    don_x, don_y, don_avg = load_donepezil_ref(
        file_key, col=col, ps_to_ns=ps_to_ns,
        smooth=smooth, smooth_window=smooth_window, is_rmsf=is_rmsf
    )

    # ---- Pass 1: load every series once, so the y-axis can be sized to the
    #      actual data (shared across panels, with light padding) instead of
    #      a fixed range that leaves a lot of empty white space ----
    drug_data, combined = {}, []
    for drug in DRUGS:
        reps = {}
        for rep in REPLICATES:
            fpath = os.path.join(BASE, drug, rep, file_key)
            if not os.path.exists(fpath):
                print(f"  Missing: {fpath}")
                continue
            x, v = read_xvg(fpath, col=col, ps_to_ns=ps_to_ns)
            v_plot = rolling_mean(v, smooth_window) if smooth else v
            reps[rep] = (x, v, v_plot)
            combined.append(v_plot)
        drug_data[drug] = reps
    if don_y is not None:
        combined.append(don_y)

    combined    = np.concatenate(combined)
    y_lo        = np.percentile(combined, 1)
    y_hi        = combined.max()
    pad         = (y_hi - y_lo) * 0.08
    shared_ylim = (y_lo - pad, y_hi + pad)

    # ---- Layout: 3 panels across the top row, 2 centered on the bottom row ----
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 6, left=0.055, right=0.97, top=0.90, bottom=0.10,
                          hspace=0.40, wspace=0.55)
    axes_flat = [
        fig.add_subplot(gs[0, 0:2]),
        fig.add_subplot(gs[0, 2:4]),
        fig.add_subplot(gs[0, 4:6]),
        fig.add_subplot(gs[1, 1:3]),
        fig.add_subplot(gs[1, 3:5]),
    ]

    for idx, drug in enumerate(DRUGS):
        ax = axes_flat[idx]
        all_vals = []

        # --- Donepezil reference line (drawn first so it sits behind) ---
        if don_y is not None:
            if is_rmsf:
                don_label = f"Donepezil (avg: {don_avg:.3f})"
            else:
                don_label = f"Donepezil (avg: {don_avg:.2f})"
            ax.plot(don_x, don_y,
                    color=DONEPEZIL_COLOR, linewidth=1.2, alpha=0.9,
                    label=don_label, zorder=1)

        # --- Drug replicates ---
        for rep in REPLICATES:
            if rep not in drug_data[drug]:
                continue
            x, v, v_plot = drug_data[drug][rep]
            avg_rep = np.mean(v)
            all_vals.append(v)
            if is_rmsf:
                label = f"{rep} (avg: {avg_rep:.3f})"
            else:
                label = f"{rep} (avg: {avg_rep:.2f})"
            ax.plot(x, v_plot, color=REP_COLORS[rep],
                    linewidth=1.0, alpha=0.95, label=label, zorder=2)

        # --- Overall drug average dashed line ---
        if all_vals and not is_rmsf:
            min_len     = min(len(v) for v in all_vals)
            stacked     = np.array([v[:min_len] for v in all_vals])
            overall_avg = np.mean(stacked)
            ax.axhline(overall_avg, color=AVG_COLOR, linestyle="--",
                       linewidth=1.2, label=f"Overall Avg: {overall_avg:.2f}",
                       zorder=3)

        ax.set_ylim(shared_ylim)
        ax.set_title(f"({PANEL_LETTERS[idx]}) {drug}", fontsize=12, fontweight="bold", loc="left")
        ax.set_xlabel("Residue" if is_rmsf else "Time (ns)", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_xlim(left=100 if is_rmsf else 0)
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7.5, loc="upper left", framealpha=0.7)
        ax.grid(False)
        open_box(ax)

    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.965)

    out_path = os.path.join(OUT, filename)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Generate all plots
# ---------------------------------------------------------------------------

print("=== RMSD ===")
make_figure(
    title    = "Backbone RMSD — MD Simulations (Triplicates)",
    ylabel   = "RMSD (nm)",
    filename = "RMSD_triplicates.png",
    file_key = "rmsd_backbone.xvg",
    col=1, ps_to_ns=False,
)

print("=== SASA ===")
make_figure(
    title    = "Solvent Accessible Surface Area — MD Simulations (Triplicates)",
    ylabel   = "SASA (nm²)",
    filename = "SASA_triplicates.png",
    file_key = "sasa.xvg",
    col=1, ps_to_ns=True,
)

print("=== Radius of Gyration ===")
make_figure(
    title    = "Radius of Gyration — MD Simulations (Triplicates)",
    ylabel   = "Rg (nm)",
    filename = "Gyration_triplicates.png",
    file_key = "gyrate.xvg",
    col=1, ps_to_ns=True,
)

print("=== RMSF per Residue ===")
make_figure(
    title    = "RMSF per Residue — MD Simulations (Triplicates)",
    ylabel   = "RMSF (nm)",
    filename = "RMSF_triplicates.png",
    file_key = "rmsf_residues.xvg",
    col=1, ps_to_ns=False, is_rmsf=True,
)

print("=== Protein–Ligand Distance ===")
make_figure(
    title    = "Protein–Ligand COM Distance — MD Simulations (Triplicates)",
    ylabel   = "Distance (nm)",
    filename = "ProtLigDist_triplicates.png",
    file_key = "prot_lig_dist.xvg",
    col=1, ps_to_ns=True,
)

print("=== Total Energy ===")
make_figure(
    title          = "Total Energy — MD Simulations (Triplicates)",
    ylabel         = "Total Energy (kJ/mol)",
    filename       = "TotalEnergy_triplicates.png",
    file_key       = "total_energy.xvg",
    col=3, ps_to_ns=True,
    smooth=True, smooth_window=500,
)

print("\nAll plots done.")
