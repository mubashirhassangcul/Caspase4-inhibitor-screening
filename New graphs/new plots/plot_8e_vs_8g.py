import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = HERE

COMPOUNDS = ["8e", "8g"]
COLORS    = {"8e": "#5B9BD5", "8g": "#ED7D31"}  # blue, orange


def open_box(ax):
    """Drop the top and right spines for an open-box look."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


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


def find_legend_col(path, label):
    """Return the data column index (1-based, time is col 0) whose
    '@ sN legend "label"' line matches `label`. Needed because the total
    energy files carry different sets/orders of energy terms per compound."""
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("@") and "legend" in line and line.split()[1].startswith("s"):
                idx = int(line.split()[1][1:])
                if f'"{label}"' in line:
                    return idx + 1
    raise ValueError(f"Legend '{label}' not found in {path}")


def rolling_mean(arr, window=200):
    kernel = np.ones(window) / window
    return np.convolve(arr, kernel, mode="same")


def make_figure(title, ylabel, filename, file_key,
                 col=1, ps_to_ns=False, legend_label=None,
                 smooth=False, smooth_window=200, is_rmsf=False):

    data = {}
    for compound in COMPOUNDS:
        fpath = os.path.join(HERE, compound, file_key)
        if not os.path.exists(fpath):
            print(f"  Missing: {fpath}")
            continue
        use_col = col
        if legend_label is not None:
            use_col = find_legend_col(fpath, legend_label)
        x, v = read_xvg(fpath, col=use_col, ps_to_ns=ps_to_ns)
        v_plot = rolling_mean(v, smooth_window) if smooth else v
        data[compound] = (x, v, v_plot)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    for compound in COMPOUNDS:
        if compound not in data:
            continue
        x, v, v_plot = data[compound]
        avg = np.mean(v)
        fmt = f"{avg:.3f}" if is_rmsf else f"{avg:.2f}"
        ax.plot(x, v_plot, color=COLORS[compound], linewidth=1.2, alpha=0.95,
                label=f"{compound} (avg: {fmt})", zorder=2)

    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel("Residue" if is_rmsf else "Time (ns)", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(labelsize=9)
    ax.legend(fontsize=9, loc="best", framealpha=0.7)
    ax.grid(False)
    open_box(ax)

    fig.tight_layout()
    out_path = os.path.join(OUT, filename)
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Generate all plots
# ---------------------------------------------------------------------------

print("=== RMSD ===")
make_figure(
    title    = "Backbone RMSD — 8e vs 8g",
    ylabel   = "RMSD (nm)",
    filename = "RMSD_8e_vs_8g.png",
    file_key = "rmsd_backbone.xvg",
    col=1, ps_to_ns=False,
)

print("=== SASA ===")
make_figure(
    title    = "Solvent Accessible Surface Area — 8e vs 8g",
    ylabel   = "SASA (nm²)",
    filename = "SASA_8e_vs_8g.png",
    file_key = "sasa.xvg",
    col=1, ps_to_ns=True,
)

print("=== Radius of Gyration ===")
make_figure(
    title    = "Radius of Gyration — 8e vs 8g",
    ylabel   = "Rg (nm)",
    filename = "Gyration_8e_vs_8g.png",
    file_key = "gyrate.xvg",
    col=1, ps_to_ns=True,
)

print("=== RMSF per Residue ===")
make_figure(
    title    = "RMSF per Residue — 8e vs 8g",
    ylabel   = "RMSF (nm)",
    filename = "RMSF_8e_vs_8g.png",
    file_key = "rmsf_residues.xvg",
    col=1, ps_to_ns=False, is_rmsf=True,
)

print("=== Total Energy ===")
make_figure(
    title        = "Total Energy — 8e vs 8g",
    ylabel       = "Total Energy (kJ/mol)",
    filename     = "TotalEnergy_8e_vs_8g.png",
    file_key     = "total_energy.xvg",
    ps_to_ns=True, legend_label="Total Energy",
    smooth=True, smooth_window=500,
)

print("\nAll plots done.")
