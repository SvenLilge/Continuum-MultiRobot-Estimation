#!/usr/bin/env python3
"""Run the Cosserat-prior evaluation and make all of its figures.

    python scripts/evaluate.py run     # run the C++ driver for every R_v and force (~1 min)
    python scripts/evaluate.py plot    # figures + summary table from the saved data (seconds)
    python scripts/evaluate.py all     # both

Run from anywhere; paths are relative to the repository root. Build the driver
first (cmake -DUSE_LOCAL_TDCR=ON ...). Needs numpy, pandas and matplotlib only.

Outputs (doc/figures/):
    data/trials_rv<R_v>.csv       one row per trial and method, incl. per-node errors
    data/shape_rv<R_v>_fz<F>.csv  estimated shapes of one trial (for the overlays)
    summary.csv                   mean, std, median, IQR, convergence per method/force
    fig1_problem                  the estimation problem: model, real shape, force, sensor
    fig2_pipeline                 how ground truth, measurement and estimates are produced
    fig3_ranking                  which method is best: average error vs. ground truth
    fig4_shapes                   estimated shapes drawn against the ground truth
    fig5_curvature                why: the bending curvature each method produces
    fig6_error_vs_force           how each method's error grows with the unknown force
    fig7_error_along_backbone     where along the rod the error occurs
    fig8_error_distribution       spread over the 50 trials (box plots)
    fig9_runtime                  time per estimate and convergence rate
Every figure is a 300 dpi PNG of the same size (IEEE page width x 3.4 in)
with the same 8 pt text, so they look alike side by side.

Plot conventions follow common practice for comparing estimators: per-trial
distributions on the same noise seeds (median + IQR, box plots with 5-95th
percentile whiskers) rather than mean +/- std alone, a colorblind-safe
Okabe-Ito palette with marker/line-style redundancy.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

REPO = Path(__file__).resolve().parents[1]
DRIVER = REPO / "examples" / "evaluation_cosserat_priors"
BASE_CONFIG = REPO / "config" / "7_evaluation_s_shape.yaml"
CONFIG_DIR = REPO / "build" / "eval_configs"
OUT_DIR = REPO / "doc" / "figures"
DATA_DIR = OUT_DIR / "data"

DEFAULT_RV = ["0.01", "0.1", "0.5"]
DEFAULT_FORCES = ["0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"]

# One size and one text size for every figure, so they match when shown together.
PAGE = 7.16                  # IEEE two-column page width [inch]
FIG_SIZE = (PAGE, 3.4)
FS = 8                       # all text [pt]

# Okabe-Ito colors (validated colorblind-safe for all pairs). Method 2 uses a
# light-to-dark ramp of its hue for increasing R_v (looser trust in the model).
INK, INK_2, GRID = "#1a1a1a", "#555555", "#dddddd"
MODEL_GRAY = "#8c8c8c"
ERROR_LABEL = "Error vs. real shape [mm]"
M1_COLOR, M3_COLOR = "#D55E00", "#0072B2"
M2_RAMP = ["#79c9ad", "#009E73", "#00553e"]


# --------------------------------------------------------------------------- run

def write_config(rv):
    """Copy the base config with R_v = R_u = rv; return its path."""
    text = BASE_CONFIG.read_text()
    text = re.sub(r"(\bR_v:\s*)[0-9.eE+-]+", rf"\g<1>{rv}", text)
    text = re.sub(r"(\bR_u:\s*)[0-9.eE+-]+", rf"\g<1>{rv}", text)
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    path = CONFIG_DIR / f"rv_{rv}.yaml"
    path.write_text(text)
    return path


def run_driver(config, fz, args, extra):
    cmd = [str(DRIVER), str(config), "--fz", fz, "--sigma", str(args.sigma)] + extra
    res = subprocess.run(cmd, capture_output=True, text=True)
    # Exit code 2 only means the last trial hit the iteration limit.
    if res.returncode not in (0, 2):
        sys.exit(f"Driver failed ({res.returncode}): {' '.join(cmd)}\n{res.stderr[-2000:]}")


def cmd_run(args):
    if not DRIVER.exists():
        sys.exit(f"Driver not found: {DRIVER}\nBuild it with: cmake -G Ninja -DCMAKE_BUILD_TYPE=Release "
                 "-DUSE_LOCAL_TDCR=ON -S . -B build && cmake --build build")
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for rv in args.rv:
        config = write_config(rv)
        trials_csv = DATA_DIR / f"trials_rv{rv}.csv"
        trials_csv.unlink(missing_ok=True)
        for fz in args.forces:
            print(f"  R_v = {rv:<5} F_z = {fz} N  ({args.trials} trials)", flush=True)
            run_driver(config, fz, args, ["--n-trials", str(args.trials), "--trials-csv", str(trials_csv)])
        for fz in args.shape_forces:
            run_driver(config, fz, args, ["--n-trials", "1", "--csv", str(DATA_DIR / f"shape_rv{rv}_fz{fz}.csv")])
    print(f"Data written to {DATA_DIR.relative_to(REPO)}")


# --------------------------------------------------------------------------- data

def load_trials():
    files = sorted(DATA_DIR.glob("trials_rv*.csv"), key=lambda p: float(p.stem[len("trials_rv"):]))
    if not files:
        sys.exit(f"No data in {DATA_DIR}. Run: python scripts/evaluate.py run")
    frames = []
    for path in files:
        df = pd.read_csv(path)
        df["rv"] = path.stem[len("trials_rv"):]
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    err_cols = sorted([c for c in df.columns if re.fullmatch(r"e\d+", c)], key=lambda c: int(c[1:]))
    df[["rmse", "max_err"] + err_cols] *= 1000.0          # m -> mm
    df["fz"] = df["fz"].round(4)
    return df, err_cols


def build_series(df):
    """List of plotted series: (key, label, short label, style, rows).

    Names are the ones a paper reader sees; they must match the text.
    Method 2's noise setting sigma (R_v = R_u in the config) is how much the
    estimator trusts the model: smaller sigma = more trust.
    """
    rvs = sorted(df["rv"].unique(), key=float)
    first = df[df["rv"] == rvs[0]]          # Baseline, force input and the model do not use R_v
    ramp = M2_RAMP if len(rvs) == len(M2_RAMP) else [M2_RAMP[1]] * len(rvs)
    trust = dict(zip(rvs, ["high", "medium", "low"])) if len(rvs) == 3 else {}
    series = [
        ("model", "Model only (no estimator)", "Model\nonly",
         dict(color=MODEL_GRAY, linestyle="--", marker=None), first[first["method"] == "0_model"]),
        ("m1", "Baseline: tip measurement only", "Baseline",
         dict(color=M1_COLOR, linestyle="--", marker="o"), first[first["method"] == "1_old"]),
    ]
    for rv, color in zip(rvs, ramp):
        rows = df[(df["rv"] == rv) & (df["method"] == "2_strain")]
        level = f"{trust[rv]} trust " if rv in trust else ""
        series.append((f"m2_{rv}", rf"Strain as measurement, {level}($\sigma$ = {rv})",
                       f"Strain\n{level.replace(' ', chr(10)).strip()}" if level else f"Strain\nσ={rv}",
                       dict(color=color, linestyle="-", marker="s"), rows))
    series.append(("m3", "Force as input", "Force\ninput",
                   dict(color=M3_COLOR, linestyle="-", marker="^"), first[first["method"] == "3_force"]))
    return series


def write_summary(series):
    rows = []
    for key, _, _, _, d in series:
        for fz, g in d.groupby("fz"):
            r = g["rmse"]
            rows.append(dict(series=key, fz=fz, n=len(g), rmse_mean=r.mean(), rmse_std=r.std(ddof=1) if len(g) > 1 else 0.0,
                             rmse_median=r.median(), rmse_q25=r.quantile(0.25), rmse_q75=r.quantile(0.75),
                             time_median_ms=g["time_ms"].median(), converged_frac=g["converged"].mean()))
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / "summary.csv", index=False, float_format="%.4f")

    table = summary.pivot(index="fz", columns="series", values="rmse_mean")
    order = [s[0] for s in series]
    table = table[[c for c in order if c in table.columns]]
    print("\nMean RMSE [mm] per force (rows) and series (columns):\n")
    print("| F (N) | " + " | ".join(table.columns) + " |")
    print("|---:|" + "---:|" * len(table.columns))
    for fz, row in table.iterrows():
        print(f"| {fz:.2f} | " + " | ".join(f"{v:.2f}" for v in row.values) + " |")


# --------------------------------------------------------------------------- style

def set_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": FS, "axes.labelsize": FS, "axes.titlesize": FS,
        "legend.fontsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS,
        "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK_2,
        "xtick.color": INK_2, "ytick.color": INK_2,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.direction": "in", "ytick.direction": "in",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "lines.linewidth": 1.3, "lines.markersize": 4,
        "legend.frameon": False,
        "figure.constrained_layout.use": True,
        "savefig.dpi": 300,
    })


def log_axis(ax, axis="y"):
    fmt = FuncFormatter(lambda v, _: f"{v:g}")
    loc = LogLocator(base=10, subs=(1.0, 2.0, 5.0))
    target = ax.yaxis if axis == "y" else ax.xaxis
    (ax.set_yscale if axis == "y" else ax.set_xscale)("log")
    target.set_major_locator(loc)
    target.set_major_formatter(fmt)
    target.set_minor_formatter(NullFormatter())


def save(fig, name):
    fig.savefig(OUT_DIR / f"{name}.png")
    plt.close(fig)
    print(f"  {name}.png")


def legend_below(fig, axes, ncol):
    handles, labels = axes.get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=ncol, handlelength=2.6, columnspacing=1.2)


# --------------------------------------------------------------------------- figures

SHOW_HEADLINES = True        # figure titles stating the question; off with --paper


def headline(fig, title, subtitle=None):
    """Title (bold) + optional subtitle above the figure, left-aligned.

    Makes each PNG understandable on its own. Use --paper to leave them out
    when the figure goes into a paper, where the caption carries this text.
    """
    if not SHOW_HEADLINES:
        return
    bold = r"$\mathbf{" + title.replace(" ", r"\ ") + "}$"
    text = bold + ("\n" + subtitle if subtitle else "")
    fig.suptitle(text, x=0.0, ha="left", fontsize=FS, color=INK, linespacing=1.4)


def mark_ground_truth(ax):
    """Label the zero line: zero distance from the real shape = the ground truth.

    The label sits just below the line, where errors (>= 0) never go; call
    room_below_zero() after all data is plotted to make space for it.
    """
    ax.axhline(0, color=INK, linewidth=1.4, zorder=1)
    ax.annotate("ground truth (zero error)", (1, 0), xycoords=("axes fraction", "data"),
                xytext=(-2, -3), textcoords="offset points", ha="right", va="top", fontsize=FS, color=INK)


def room_below_zero(ax):
    top = ax.get_ylim()[1]
    ax.set_ylim(-0.12 * top, top)


def panel(ax, i, text):
    """Paper-style panel label: (a), (b), ... on the left above the axes."""
    ax.set_title(f"({chr(97 + i)}) {text}", loc="left")


def load_shape(fz, rv):
    """Shape of one trial: node rows and the noisy tip measurement (node = -1)."""
    path = DATA_DIR / f"shape_rv{rv}_fz{fz}.csv"
    if not path.exists():
        return None, None
    raw = pd.read_csv(path)
    return raw[raw["node"] >= 0], raw[raw["node"] == -1]


def bending_plane(shapes):
    """The two axes with the largest spread of the true and predicted shapes."""
    allpts = pd.concat([s[[f"{w}_{c}" for w in ("prior", "truth") for c in "xyz"]] for s in shapes])
    spread = {c: np.ptp(allpts[[f"prior_{c}", f"truth_{c}"]].values) for c in "xyz"}
    return tuple(sorted(sorted(spread, key=spread.get)[-2:]))


def fig_problem(series, fz, rv):
    """Concept figure: what the estimator must recover, and from what."""
    s, tip = load_shape(fz, rv)
    if s is None:
        print("  (skipping fig1_problem, no shape data)")
        return
    a, b = bending_plane([s])
    fig, ax = plt.subplots(figsize=FIG_SIZE, layout="constrained")
    headline(fig, "The estimation problem",
             "Recover the real shape (black) from one noisy tip position (star)")
    P = lambda col: (s[f"{col}_{a}"].values * 1000.0, s[f"{col}_{b}"].values * 1000.0)
    xm, zm = P("prior")
    xt, zt = P("truth")
    ax.plot(xm, zm, color=MODEL_GRAY, linestyle="--", linewidth=1.5)
    ax.plot(xt, zt, color=INK, linewidth=2.2, marker="o", markersize=3)
    tx, tz = tip[f"m1_old_{a}"].iloc[0] * 1000.0, tip[f"m1_old_{b}"].iloc[0] * 1000.0
    ax.plot(tx, tz, marker="*", markersize=12, markerfacecolor="white", markeredgecolor=INK,
            markeredgewidth=0.9, linestyle="none", zorder=5)
    f = force_in_plot_frame(series, fz)
    d = np.array([f[a], f[b]]) / np.hypot(f[a], f[b])
    end = np.array([xt[-1], zt[-1]]) + 22 * d
    ax.plot(*end, alpha=0)
    ax.annotate("", xy=end, xytext=(xt[-1], zt[-1]),
                arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=1.2, mutation_scale=9))
    ax.plot(0, 0, marker="s", color=INK, markersize=6, linestyle="none")
    line = dict(arrowstyle="-", color=INK_2, linewidth=0.6, shrinkA=1, shrinkB=3)
    ax.annotate("Model prediction\n(knows the tendons,\nnot the force)", xy=(xm[5], zm[5]), xytext=(4, 66),
                fontsize=FS, va="center", arrowprops=line)
    ax.annotate("Tip sensor, ±1 mm:\nthe only measurement", xy=(tx, tz), xytext=(236, 92),
                fontsize=FS, ha="right", va="center", arrowprops=line)
    ax.annotate(rf"Unknown force $F$ = {float(fz):.2f} N", xy=end, xytext=(236, 16),
                fontsize=FS, ha="right", va="center", arrowprops=line)
    ax.annotate("Real shape: what the estimator\nmust recover from the star", xy=(xt[3], zt[3]), xytext=(40, -16),
                fontsize=FS, va="center", arrowprops=line)
    ax.annotate("fixed\nbase", (0, 0), xytext=(-2, 5), textcoords="offset points", fontsize=FS,
                color=INK_2, ha="right", va="bottom")
    ax.set_xlim(-30, 273)                       # wide enough to fill the canvas at equal aspect
    ax.set_ylim(-28, 102)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"{a} [mm]")
    ax.set_ylabel(f"{b} [mm]")
    save(fig, "fig1_problem")


def fig_pipeline():
    """Concept figure: how ground truth, measurement and estimates are produced."""
    from matplotlib.patches import FancyBboxPatch
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    headline(fig, "How the simulated experiment works")
    ax.set_xlim(0, 100)
    ax.set_ylim(-5, 33)
    ax.axis("off")

    def box(x, y, w, h, text, edge=INK_2, face="#f4f4f2", lw=0.8, weight="normal"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                    facecolor=face, edgecolor=edge, linewidth=lw))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=FS, color=INK, weight=weight)

    def arrow(x0, y0, x1, y1, text=None, dy=1.2):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=INK_2, linewidth=0.8, mutation_scale=7,
                                    shrinkA=0, shrinkB=0))
        if text:
            ax.text((x0 + x1) / 2, (y0 + y1) / 2 + dy, text, ha="center", va="bottom", fontsize=FS, color=INK_2)

    box(0.5, 13, 11, 7, "Known tendon\ntensions $q$")
    box(17, 23, 17, 7.5, "Cosserat model\nwithout force")
    box(17, 2.5, 17, 7.5, "Cosserat model\nwith unknown force $F$", edge=INK, lw=1.1)
    box(40, 2.5, 13, 7.5, "Tip position\n+ noise (1 mm)")
    box(60, 9.5, 17, 14, "Estimator\n\n" + "Baseline\nStrain as measurement\nForce as input", lw=1.0)
    box(84, 9.5, 15.5, 14, "Compare\nwith ground truth\n\n→ shape error [mm]")
    arrow(11.9, 18.5, 16.6, 26.5)
    arrow(11.9, 14.5, 16.6, 6.5)
    arrow(34.4, 26.75, 59.6, 21.5, "model prediction")
    arrow(34.4, 6.25, 39.6, 6.25)
    arrow(53.4, 6.25, 59.6, 12.0)
    ax.text(56.0, 4.2, "noisy tip\nmeasurement", ha="left", va="top", fontsize=FS, color=INK_2)
    arrow(77.4, 16.5, 83.6, 16.5, "estimate")
    ax.plot([25.5, 25.5, 91.75], [2.0, -1.8, -1.8], color=INK, linewidth=1.1, solid_joinstyle="round")
    ax.annotate("", xy=(91.75, 9.0), xytext=(91.75, -1.8),
                arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=1.1, mutation_scale=7, shrinkA=0, shrinkB=0))
    ax.text(58, -2.4, "ground-truth shape (hidden from the estimator)", ha="center", va="top",
            fontsize=FS, color=INK)
    save(fig, "fig2_pipeline")


def curvature_from_positions(s, col, a, b):
    """Signed bending curvature [1/m] at interior nodes, from the turning of the chords."""
    P = s[[f"{col}_{a}", f"{col}_{b}"]].values
    d = np.diff(P, axis=0)
    theta = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
    ds = np.linalg.norm(d, axis=1)
    return np.diff(theta) / (0.5 * (ds[:-1] + ds[1:]))


def fig_curvature(series, forces, shape_rv, junction):
    """Why the methods differ: the bending curvature each one produces along the rod."""
    loaded = {fz: load_shape(fz, shape_rv)[0] for fz in forces}
    if any(v is None for v in loaded.values()):
        print("  (skipping fig5_curvature, missing shape data)")
        return
    a, b = bending_plane(list(loaded.values()))
    by_key = {s[0]: s for s in series}
    lines = [
        ("truth", "Ground truth", dict(color=INK, linestyle="-", marker="o", linewidth=2.0, markersize=3)),
        ("prior", by_key["model"][1], dict(color=MODEL_GRAY, linestyle="--", marker=None)),
        ("m1_old", by_key["m1"][1], by_key["m1"][3]),
        ("m2_strain", by_key[f"m2_{shape_rv}"][1], by_key[f"m2_{shape_rv}"][3]),
        ("m3_force", by_key["m3"][1], by_key["m3"][3]),
    ]
    fig, axes = plt.subplots(1, len(forces), figsize=FIG_SIZE, sharey=True, layout="constrained")
    headline(fig, "Why the methods differ",
             "How each method bends the rod. The force weakens the bend in segment 1, but the flip at "
             "the junction stays; only Force as input gets both right")
    axes = np.atleast_1d(axes)
    for i, (ax, (fz, s)) in enumerate(zip(axes, loaded.items())):
        s_mm = s["s"].values[1:-1] * 1000.0
        for col, label, st in lines:
            ax.plot(s_mm, curvature_from_positions(s, col, a, b), label=label, **st)
        ax.axhline(0, color=INK_2, linewidth=0.6)
        if junction is not None:
            ax.axvline(junction * 1000.0, color=INK_2, linestyle=":", linewidth=0.8)
            ax.annotate("segment\njunction", (junction * 1000.0, 1), xycoords=("data", "axes fraction"),
                        xytext=(2, -2), textcoords="offset points", ha="left", va="top", fontsize=FS, color=INK_2)
        panel(ax, i, rf"Unknown tip force $F$ = {float(fz):.2f} N (one trial)")
        ax.set_xlabel("Arc length along the rod, base to tip [mm]")
    axes[0].set_ylabel("Bending curvature [1/m]\n(+ bends one way, − the other)")
    legend_below(fig, axes[0], ncol=3)
    save(fig, "fig5_curvature")


def fig_ranking(series):
    """Headline: which method recovers the shape best (mean error vs. ground truth)."""
    model = next(s for s in series if s[0] == "model")
    rows = []
    for key, label, _, st, d in series:
        if key == "model":
            continue
        name = label
        if key.startswith("m2_"):                   # "Strain as measurement, high trust (σ = 0.01)"
            head, tail = label.split(", ", 1)
            name = f"{head}\n{tail.replace(' trust', ' trust in model')}"
        rows.append((name, st["color"], d["rmse"].mean()))
    rows.sort(key=lambda r: r[2], reverse=True)     # barh draws bottom-up: best ends on top
    model_mean = model[4]["rmse"].mean()

    fig, ax = plt.subplots(figsize=FIG_SIZE, layout="constrained")
    headline(fig, "Which method recovers the rod shape best?",
             "Average distance between the estimated and the real rod, 350 simulated runs per method")
    y = np.arange(len(rows))
    ax.barh(y, [r[2] for r in rows], height=0.62, color=[r[1] for r in rows], edgecolor="white", linewidth=1.0)
    for yi, (_, _, v) in zip(y, rows):
        ax.annotate(f"{v:.1f} mm", (v, yi), xytext=(3, 0), textcoords="offset points",
                    va="center", fontsize=FS, color=INK)
    ax.axvline(model_mean, color=MODEL_GRAY, linestyle="--", linewidth=1.2, zorder=0)
    ax.annotate(f"model alone,\nno estimator:\n{model_mean:.1f} mm", (model_mean, len(rows) - 0.5),
                xytext=(3, 0), textcoords="offset points", ha="left", va="top", fontsize=FS, color=INK_2)
    ax.set_yticks(y, [r[0] for r in rows])
    best = ax.get_yticklabels()[-1]
    best.set_fontweight("bold")
    best.set_color(INK)
    ax.set_xlim(0, max(r[2] for r in rows) * 1.25)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel("Average shape error [mm]   (shorter bar = better)", loc="right")
    ax.grid(axis="y", visible=False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    save(fig, "fig3_ranking")


def fig_accuracy(series):
    """Median error vs. force; vertical bars span the 25th-75th percentile."""
    est = [s for s in series if s[0] != "model"]
    fig, ax = plt.subplots(figsize=FIG_SIZE, layout="constrained")
    headline(fig, "Accuracy as the unknown force grows",
             "Typical (median) distance from the real shape over 50 runs; lower is better")
    forces = np.sort(series[0][4]["fz"].unique())
    step = np.min(np.diff(forces)) if len(forces) > 1 else 0.05
    offsets = np.linspace(-0.18, 0.18, len(est)) * step         # dodge so bars do not overlap
    for key, label, _, st, d in series:
        g = d.groupby("fz")["rmse"]
        med, q1, q3 = g.median(), g.quantile(0.25), g.quantile(0.75)
        x = med.index.values
        if key == "model":
            ax.plot(x, med.values, color=st["color"], linestyle="--", linewidth=1.3, label=label, zorder=1)
            continue
        dx = offsets[[s[0] for s in est].index(key)]
        ax.errorbar(x + dx, med.values, yerr=[med.values - q1.values, q3.values - med.values],
                    capsize=1.5, elinewidth=0.7, label=label, **st)
    mark_ground_truth(ax)
    room_below_zero(ax)
    ax.set_xlabel(r"Unknown tip force $F$ [N]")
    ax.set_ylabel(ERROR_LABEL)
    legend_below(fig, ax, ncol=2)
    save(fig, "fig6_error_vs_force")


def fig_distribution(series, forces):
    methods = [s for s in series if s[0] != "model"]
    model = next(s for s in series if s[0] == "model")[4]
    fig, axes = plt.subplots(1, len(forces), figsize=FIG_SIZE, sharey=True, layout="constrained")
    headline(fig, "How consistent each method is",
             "Each dot is one run with a different sensor-noise draw; a tight cluster means a repeatable result")
    axes = np.atleast_1d(axes)
    rng = np.random.default_rng(0)
    for i, (ax, fz) in enumerate(zip(axes, forces)):
        data = [d.loc[d["fz"] == fz, "rmse"].values for *_, d in methods]
        pos = np.arange(1, len(methods) + 1)
        bp = ax.boxplot(data, positions=pos, widths=0.55, whis=(5, 95), showfliers=False,
                        patch_artist=True, medianprops=dict(color=INK, linewidth=1.0),
                        whiskerprops=dict(color=INK_2, linewidth=0.7), capprops=dict(color=INK_2, linewidth=0.7))
        for patch, (_, _, _, st, _) in zip(bp["boxes"], methods):
            patch.set(facecolor=st["color"], alpha=0.35, edgecolor=st["color"], linewidth=0.8)
        for p, vals, (_, _, _, st, _) in zip(pos, data, methods):
            ax.scatter(p + rng.uniform(-0.17, 0.17, len(vals)), vals, s=3, color=st["color"],
                       alpha=0.6, linewidths=0, zorder=3)
        m = model.loc[model["fz"] == fz, "rmse"]
        if len(m) and m.iloc[0] > 0:
            ax.axhline(m.iloc[0], color=MODEL_GRAY, linestyle="--", linewidth=1.0)
            ax.annotate("model only", (pos[-1] + 0.45, m.iloc[0]), xytext=(0, 2), textcoords="offset points",
                        ha="right", va="bottom", fontsize=FS, color=INK_2)
        ax.set_xticks(pos, [s[2] for s in methods])
        panel(ax, i, rf"Unknown tip force $F$ = {fz:.2f} N")
        ax.grid(axis="x", visible=False)
    log_axis(axes[0])
    axes[0].set_ylabel(ERROR_LABEL.replace("[mm]", "[mm], log scale"))
    save(fig, "fig8_error_distribution")


def fig_backbone_error(series, err_cols, forces, junction):
    fig, axes = plt.subplots(1, len(forces), figsize=FIG_SIZE, sharey=True, layout="constrained")
    headline(fig, "Where along the rod each method goes wrong", "Distance from the real shape at each point, base to tip")
    axes = np.atleast_1d(axes)
    for i, (ax, fz) in enumerate(zip(axes, forces)):
        for key, label, _, st, d in series:
            rows = d[d["fz"] == fz]
            if rows.empty:
                continue
            L = rows["L"].iloc[0]
            s_mm = np.linspace(0.0, L, len(err_cols)) * 1000.0
            e = rows[err_cols].values
            if key == "model":
                ax.plot(s_mm, e[0], color=st["color"], linestyle="--", label=label)
                continue
            ax.fill_between(s_mm, np.percentile(e, 25, axis=0), np.percentile(e, 75, axis=0),
                            color=st["color"], alpha=0.15, linewidth=0)
            ax.plot(s_mm, np.median(e, axis=0), label=label, **st)
        if junction is not None:
            ax.axvline(junction * 1000.0, color=INK_2, linestyle=":", linewidth=0.8)
            ax.annotate("segment\njunction", (junction * 1000.0, 1), xycoords=("data", "axes fraction"),
                        xytext=(2, -2), textcoords="offset points", ha="left", va="top", fontsize=FS, color=INK_2)
        mark_ground_truth(ax)
        panel(ax, i, rf"Unknown tip force $F$ = {fz:.2f} N")
        ax.set_xlabel("Position along the rod, base to tip [mm]")
    room_below_zero(axes[0])
    axes[0].set_ylabel("Distance from real shape [mm]")
    legend_below(fig, axes[0], ncol=3)
    save(fig, "fig7_error_along_backbone")


def force_in_plot_frame(series, fz):
    """Tip-force direction in the estimator/plot frame (assumes T_i0 = identity).

    The driver applies (fx, fy, fz) in the Cosserat frame (z = base axis); the
    estimator frame is the cyclic permutation P: (x, y, z)_est = (z, x, y)_cos.
    """
    d = series[0][4]
    row = d[np.isclose(d["fz"], float(fz))]
    if row.empty:
        return None
    fx, fy, fzc = row[["fx", "fy", "fz"]].iloc[0]
    return dict(x=fzc, y=fx, z=fy)


def fig_shapes(series, forces, shape_rv):
    files = {fz: DATA_DIR / f"shape_rv{shape_rv}_fz{fz}.csv" for fz in forces}
    missing = [str(p.name) for p in files.values() if not p.exists()]
    if missing:
        print(f"  (skipping fig2_shapes, missing {', '.join(missing)})")
        return
    # The driver appends one row with node = -1 holding the noisy tip measurement
    # in the m1_old columns.
    raw = {fz: pd.read_csv(p) for fz, p in files.items()}
    shapes = {fz: r[r["node"] >= 0] for fz, r in raw.items()}
    tips = {fz: r[r["node"] == -1] for fz, r in raw.items()}
    # Draw in the plane where the rod bends most: the two axes with the largest spread.
    allpts = pd.concat([s[[f"{w}_{c}" for w in ("prior", "truth") for c in "xyz"]] for s in shapes.values()])
    spread = {c: np.ptp(allpts[[f"prior_{c}", f"truth_{c}"]].values) for c in "xyz"}
    a, b = sorted(sorted(spread, key=spread.get)[-2:])
    by_key = {s[0]: s for s in series}
    m2 = by_key[f"m2_{shape_rv}"]
    lines = [
        ("truth", "Ground truth (real shape, with tip force)", dict(color=INK, linestyle="-", marker="o", linewidth=2.0, markersize=3)),
        ("prior", by_key["model"][1], dict(color=MODEL_GRAY, linestyle="--", marker=None)),
        ("m1_old", by_key["m1"][1], by_key["m1"][3]),
        ("m2_strain", m2[1], m2[3]),
        ("m3_force", by_key["m3"][1], by_key["m3"][3]),
    ]
    fig, axes = plt.subplots(1, len(forces), figsize=FIG_SIZE, layout="constrained")
    headline(fig, "What the estimates look like",
             "A method is good when its line lies on the black ground truth")
    axes = np.atleast_1d(axes)
    for i, (ax, (fz, s)) in enumerate(zip(axes, shapes.items())):
        for col, label, st in lines:
            ax.plot(s[f"{col}_{a}"] * 1000.0, s[f"{col}_{b}"] * 1000.0, label=label, **st)
        tip = tips[fz]
        if not tip.empty:
            ax.plot(tip[f"m1_old_{a}"] * 1000.0, tip[f"m1_old_{b}"] * 1000.0, marker="*", markersize=10,
                    markerfacecolor="white", markeredgecolor=INK, markeredgewidth=0.8, linestyle="none",
                    label="Noisy tip measurement (the only sensor)", zorder=5)
        # Arrow for the unknown tip force, drawn at the true tip.
        f = force_in_plot_frame(series, fz)
        if f is not None and np.hypot(f[a], f[b]) > 0:
            direction = np.array([f[a], f[b]]) / np.hypot(f[a], f[b])
            p_tip = np.array([s[f"truth_{a}"].iloc[-1], s[f"truth_{b}"].iloc[-1]]) * 1000.0
            end = p_tip + 18 * direction
            ax.plot(*end, alpha=0)                           # keeps the arrow inside the axes limits
            ax.annotate("", xy=end, xytext=p_tip,
                        arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=1.0, mutation_scale=8))
            ax.annotate(r"$F$", end, xytext=(-4, 4), textcoords="offset points", fontsize=FS)
        ax.plot(0, 0, marker="s", color=INK, markersize=5, linestyle="none", label="Fixed base")
        panel(ax, i, rf"Unknown tip force $F$ = {float(fz):.2f} N (one trial)")
        ax.set_xlabel(f"{a} [mm]")
    # Same limits in every panel; equal aspect then widens them to fill the panel.
    xs = [ax.get_xlim() for ax in axes]
    ys = [ax.get_ylim() for ax in axes]
    for i, ax in enumerate(axes):
        ax.set_xlim(min(x[0] for x in xs), max(x[1] for x in xs))
        ax.set_ylim(min(y[0] for y in ys), max(y[1] for y in ys))
        ax.set_aspect("equal", adjustable="datalim")
        if i > 0:
            ax.tick_params(labelleft=False)
    axes[0].set_ylabel(f"{b} [mm]")
    handles, labels = axes[0].get_legend_handles_labels()
    keep = [i for i, l in enumerate(labels) if not l.startswith("_")]
    fig.legend([handles[i] for i in keep], [labels[i] for i in keep], loc="outside lower center",
               ncol=3, handlelength=2.6, columnspacing=1.2)
    save(fig, "fig4_shapes")


def fig_runtime(series):
    methods = [s for s in series if s[0] != "model"][::-1]     # top-to-bottom = legend order
    fig, ax = plt.subplots(figsize=FIG_SIZE, layout="constrained")
    headline(fig, "Computation time per estimate", "Each dot is one run; lower is faster")
    data = [d["time_ms"].values for *_, d in methods]
    pos = np.arange(1, len(methods) + 1)
    bp = ax.boxplot(data, positions=pos, vert=False, widths=0.55, whis=(5, 95), showfliers=False,
                    patch_artist=True, medianprops=dict(color=INK, linewidth=1.0),
                    whiskerprops=dict(color=INK_2, linewidth=0.7), capprops=dict(color=INK_2, linewidth=0.7))
    for patch, (_, _, _, st, _) in zip(bp["boxes"], methods):
        patch.set(facecolor=st["color"], alpha=0.35, edgecolor=st["color"], linewidth=0.8)
    rng = np.random.default_rng(0)
    for p, vals, (_, _, _, st, _) in zip(pos, data, methods):
        ax.scatter(vals, p + rng.uniform(-0.2, 0.2, len(vals)), s=2, color=st["color"],
                   alpha=0.45, linewidths=0, zorder=3)
    log_axis(ax, "x")
    xmax = max(np.max(v) for v in data)
    ax.set_xlim(right=xmax * 5)
    for p, (_, _, _, _, d) in zip(pos, methods):
        ax.annotate(f"{100 * d['converged'].mean():.0f}% converged", (xmax * 1.4, p), va="center",
                    fontsize=FS, color=INK_2)
    ax.set_yticks(pos, [s[2].replace("\n", " ") for s in methods])
    ax.set_xlabel("Computation time per estimate [ms], log scale")
    ax.grid(axis="y", visible=False)
    save(fig, "fig9_runtime")


def cmd_plot(args):
    df, err_cols = load_trials()
    set_style()
    series = build_series(df)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Figures:")
    fig_problem(series, args.shape_forces[-1], args.shape_rv)
    fig_pipeline()
    fig_ranking(series)
    fig_shapes(series, args.shape_forces, args.shape_rv)
    fig_curvature(series, args.shape_forces, args.shape_rv, args.junction)
    fig_accuracy(series)
    fig_backbone_error(series, err_cols, [float(f) for f in args.box_forces], args.junction)
    fig_distribution(series, [float(f) for f in args.box_forces])
    fig_runtime(series)
    write_summary(series)


# --------------------------------------------------------------------------- CLI

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", choices=["run", "plot", "all"])
    p.add_argument("--rv", nargs="+", default=DEFAULT_RV, help="R_v = R_u levels for Method 2")
    p.add_argument("--forces", nargs="+", default=DEFAULT_FORCES, help="tip forces F_z [N]")
    p.add_argument("--trials", type=int, default=50, help="noise draws per force")
    p.add_argument("--sigma", type=float, default=0.001, help="tip noise std [m]")
    p.add_argument("--shape-forces", nargs="+", default=["0.10", "0.30"], help="forces for the shape overlays")
    p.add_argument("--box-forces", nargs="+", default=["0.00", "0.15", "0.30"], help="forces for figs 7 and 8")
    p.add_argument("--shape-rv", default="0.5", help="R_v of Method 2 shown in the shape overlays")
    p.add_argument("--junction", type=float, default=0.10, help="segment junction arc length [m]; marked in figs 5 and 7")
    p.add_argument("--paper", action="store_true", help="no figure titles (the paper caption carries them)")
    args = p.parse_args()
    global SHOW_HEADLINES
    SHOW_HEADLINES = not args.paper
    if args.command in ("run", "all"):
        cmd_run(args)
    if args.command in ("plot", "all"):
        cmd_plot(args)


if __name__ == "__main__":
    main()
