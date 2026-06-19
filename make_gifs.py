#!/usr/bin/env python3
"""
make_gifs.py  —  generate animated GIFs and copy figures for the GitHub Pages site.

Run from the repo root AFTER the analysis pipeline has been run:
    python3 make_gifs.py

Outputs → docs/images/
"""

import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from pathlib import Path

ROOT  = Path(".")
FEAT  = ROOT / "data" / "features"
FIGS  = ROOT / "output" / "figures"
DOCS  = ROOT / "docs" / "images"
DOCS.mkdir(parents=True, exist_ok=True)

Y_MID = 61.0

# ── colour palette (matches website) ──────────────────────────────────────────
C_BOTTOM  = "#e07a5f"   # warm orange-red  (bottom half)
C_TOP     = "#3d405b"   # dark slate-blue  (top half)
C_RUN1    = "#e07a5f"   # Run 1
C_FAM     = "#3d405b"   # Runs 2-6
C_BG      = "#f7f4ef"   # page background

FISH_COLORS = [
    "#e07a5f","#3d405b","#81b29a","#f2cc8f","#6a994e",
    "#e76f51","#264653","#2a9d8f","#a8dadc","#457b9d","#e63946","#588157"
]


# ══════════════════════════════════════════════════════════════════════════════
# GIF 1 — Single fish trajectory
# ══════════════════════════════════════════════════════════════════════════════
def make_trajectory_gif():
    print("GIF 1: fish trajectory …")
    df = pd.read_csv(FEAT / "tracking_interp.csv")
    track = (df[(df["fish_id"] == 21) & (df["run"] == 1)]
               .sort_values("frame").reset_index(drop=True))

    # Subsample to ~120 frames (20 min compressed to ~15 s at 8 fps)
    step  = max(1, len(track) // 120)
    samp  = track.iloc[::step].reset_index(drop=True)
    N     = len(samp)
    TRAIL = 12                          # number of trail points to show

    x_all = samp["x_interp"].values
    y_all = samp["y_interp"].values
    t_all = samp["time_min"].values

    xlo, xhi = 8.0,  250.0
    ylo, yhi = 0.0,  125.0

    fig, ax = plt.subplots(figsize=(3.2, 4.5))
    fig.patch.set_facecolor(C_BG)
    ax.set_facecolor(C_BG)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_aspect("equal")
    ax.axis("off")

    # Zone fills
    ax.fill_between([xlo, xhi], ylo, Y_MID, color=C_BOTTOM, alpha=0.18, zorder=0)
    ax.fill_between([xlo, xhi], Y_MID, yhi,  color=C_TOP,    alpha=0.18, zorder=0)

    # Tank border
    for spine in ["top","bottom","left","right"]:
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_color("#888")
        ax.spines[spine].set_linewidth(1)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    # Midline
    ax.axhline(Y_MID, color="#888", lw=1, ls="--", zorder=1)
    ax.text(xhi - 4, Y_MID + 2, "midline", fontsize=6, color="#888",
            ha="right", va="bottom")

    # Zone labels
    ax.text((xlo+xhi)/2, (Y_MID+yhi)/2 + 10, "TOP HALF",
            ha="center", va="center", fontsize=8, color=C_TOP,
            fontweight="bold", alpha=0.6, zorder=2)
    ax.text((xlo+xhi)/2, Y_MID/2 - 6, "BOTTOM HALF",
            ha="center", va="center", fontsize=8, color=C_BOTTOM,
            fontweight="bold", alpha=0.6, zorder=2)

    trail_line, = ax.plot([], [], lw=2.5, color=C_TOP, alpha=0.35, zorder=3)
    fish_dot,   = ax.plot([], [], "o", ms=9, color=C_TOP,
                          markeredgecolor="white", markeredgewidth=1.2, zorder=5)
    timer_text  = ax.text(xlo + 4, yhi - 6, "", fontsize=8,
                          color="#333", ha="left", va="top", zorder=6)

    def init():
        trail_line.set_data([], [])
        fish_dot.set_data([], [])
        timer_text.set_text("")
        return trail_line, fish_dot, timer_text

    def update(i):
        lo = max(0, i - TRAIL)
        xs = x_all[lo:i+1]
        ys = y_all[lo:i+1]
        trail_line.set_data(xs, ys)
        xi, yi = x_all[i], y_all[i]
        if np.isfinite(xi) and np.isfinite(yi):
            fish_dot.set_data([xi], [yi])
            fish_dot.set_color(C_TOP if yi >= Y_MID else C_BOTTOM)
        minutes = t_all[i]
        timer_text.set_text(f"min {minutes:.0f}")
        return trail_line, fish_dot, timer_text

    ani = animation.FuncAnimation(fig, update, frames=N, init_func=init,
                                  interval=125, blit=True)
    out = DOCS / "gif_trajectory.gif"
    ani.save(out, writer=PillowWriter(fps=8))
    plt.close(fig)
    print(f"  → {out}  ({out.stat().st_size//1024} KB)")


# ══════════════════════════════════════════════════════════════════════════════
# GIF 2 — Minute-by-minute curve building up (Run 1 vs Runs 2-6)
# ══════════════════════════════════════════════════════════════════════════════
def make_habituation_gif():
    print("GIF 2: habituation build-up …")
    df = pd.read_csv(FEAT / "features_by_minute.csv")
    df["pct_top"] = 1.0 - df["pct_bottom"]

    mins = sorted(df["minute"].unique())   # 0–20

    def curve(sub):
        g  = sub.groupby("minute")["pct_top"]
        m  = g.mean()
        se = g.std() / np.sqrt(g.count())
        return [m.get(mi, np.nan) for mi in mins], [se.get(mi, np.nan) for mi in mins]

    r1_means,  r1_se  = curve(df[df["run"] == 1])
    fam_means, fam_se = curve(df[df["run"].isin([2,3,4,5,6])])

    fig, ax = plt.subplots(figsize=(6, 3.6))
    fig.patch.set_facecolor(C_BG)
    ax.set_facecolor(C_BG)
    ax.set_xlim(-0.5, 20.5)
    ax.set_ylim(0.25, 0.85)
    ax.axhline(0.5, color="#bbb", lw=0.8, ls=":")
    ax.set_xlabel("Minute into session", fontsize=9, color="#333")
    ax.set_ylabel("Proportion of time in top half", fontsize=9, color="#333")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    l1_fill = ax.fill_between([], [], [], color=C_RUN1, alpha=0.15)
    l2_fill = ax.fill_between([], [], [], color=C_FAM,  alpha=0.15)
    l1_line, = ax.plot([], [], "-o", ms=4, lw=2, color=C_RUN1,
                       label="Run 1  (novel tank)")
    l2_line, = ax.plot([], [], "-o", ms=4, lw=2, color=C_FAM,
                       label="Runs 2–6  (familiar)")
    ax.legend(fontsize=8, loc="upper left", frameon=False)

    # Extra hold frames at the end
    HOLD = 10
    n_frames = len(mins) + HOLD

    def update(frame):
        nonlocal l1_fill, l2_fill
        fi = min(frame, len(mins) - 1)
        xd = [m for m in mins[:fi+1] if np.isfinite(r1_means[m])]
        y1 = [r1_means[m]  for m in mins[:fi+1] if np.isfinite(r1_means[m])]
        s1 = [r1_se[m]     for m in mins[:fi+1] if np.isfinite(r1_means[m])]
        y2 = [fam_means[m] for m in mins[:fi+1] if np.isfinite(fam_means[m])]
        s2 = [fam_se[m]    for m in mins[:fi+1] if np.isfinite(fam_means[m])]
        l1_line.set_data(xd, y1)
        l2_line.set_data(xd, y2)

        l1_fill.remove(); l2_fill.remove()
        l1_fill = ax.fill_between(
            xd, [a-b for a,b in zip(y1,s1)], [a+b for a,b in zip(y1,s1)],
            color=C_RUN1, alpha=0.15)
        l2_fill = ax.fill_between(
            xd, [a-b for a,b in zip(y2,s2)], [a+b for a,b in zip(y2,s2)],
            color=C_FAM, alpha=0.15)
        return l1_line, l2_line

    ani = animation.FuncAnimation(fig, update, frames=n_frames,
                                  interval=250, blit=False)
    out = DOCS / "gif_habituation.gif"
    ani.save(out, writer=PillowWriter(fps=4))
    plt.close(fig)
    print(f"  → {out}  ({out.stat().st_size//1024} KB)")


# ══════════════════════════════════════════════════════════════════════════════
# GIF 3 — Rank consistency across 6 runs (slope chart)
# ══════════════════════════════════════════════════════════════════════════════
def make_ranks_gif():
    print("GIF 3: rank stability …")
    df = pd.read_csv(FEAT / "features_by_minute.csv")
    df["pct_top"] = 1.0 - df["pct_bottom"]

    run_top = (df.groupby(["fish_id","run"])["pct_top"]
                 .mean().reset_index())

    # Pick 12 fish with all 6 runs
    have_all = run_top.groupby("fish_id")["run"].nunique()
    fish_ids = list(have_all[have_all == 6].index[:12])
    data = run_top[run_top["fish_id"].isin(fish_ids)]
    wide = data.pivot(index="fish_id", columns="run", values="pct_top")

    # Sort by Run-1 value for consistent ordering
    wide = wide.loc[wide[1].sort_values().index]
    fish_list = list(wide.index)
    n_fish = len(fish_list)

    RUNS   = [1, 2, 3, 4, 5, 6]
    HOLD   = 18   # frames held at each run
    TRANS  = 10   # frames for transition
    STAGES = len(RUNS)

    # Precompute all positions: shape (n_fish, n_runs)
    Y = wide[RUNS].values       # (n_fish, 6)

    # x positions for the 6 run columns
    X_POS = np.arange(1, 7, dtype=float)
    RUN_LABELS = ["Run 1\n(Day 1 AM)", "Run 2\n(Day 1 PM)",
                  "Run 3\n(Day 2 AM)", "Run 4\n(Day 2 PM)",
                  "Run 5\n(Day 3 AM)", "Run 6\n(Day 3 PM)"]

    # Build frame sequence: (current_stage, progress 0→1)
    frame_seq = []
    for s in range(STAGES):
        for _ in range(HOLD):
            frame_seq.append((s, 1.0))
        if s < STAGES - 1:
            for t in range(TRANS):
                frame_seq.append((s, (t + 1) / TRANS))
    # Extra hold at end
    for _ in range(25):
        frame_seq.append((STAGES - 1, 1.0))
    N = len(frame_seq)

    fig, ax = plt.subplots(figsize=(7, 4))
    fig.patch.set_facecolor(C_BG)
    ax.set_facecolor(C_BG)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlim(0.5, 6.5)
    ax.set_xticks(X_POS)
    ax.set_xticklabels(RUN_LABELS, fontsize=7.5, color="#444")
    ax.set_ylabel("Proportion of time\nin top half  (pct_top)", fontsize=9)
    ax.axhline(0.5, color="#ccc", lw=0.8, ls=":")
    ax.text(6.4, 0.52, "50%", fontsize=7, color="#aaa", ha="right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(bottom=False)

    # Pre-create one line + dot per fish
    lines = []
    dots  = []
    for fi, fid in enumerate(fish_list):
        c = FISH_COLORS[fi % len(FISH_COLORS)]
        ln, = ax.plot([], [], "-", lw=1.8, color=c, alpha=0.7, zorder=2)
        dt, = ax.plot([], [], "o", ms=8, color=c, zorder=3,
                      markeredgecolor="white", markeredgewidth=0.8)
        lines.append(ln); dots.append(dt)

    title_text = ax.set_title("", fontsize=10, pad=8, color="#222")

    run_labels_shown = []

    def update(frame_idx):
        stage, progress = frame_seq[frame_idx]

        # How many runs to draw lines for
        n_cols = stage + 1          # columns revealed so far
        x_visible = X_POS[:n_cols]

        # Current dot x-position (interpolating during transition)
        if progress < 1.0 and stage < STAGES - 1:
            dot_x = X_POS[stage] + progress * (X_POS[stage+1] - X_POS[stage])
        else:
            dot_x = X_POS[stage]

        for fi in range(n_fish):
            y_vals = Y[fi, :n_cols]
            lines[fi].set_data(x_visible, y_vals)

            # Interpolate dot y during transition
            if progress < 1.0 and stage < STAGES - 1:
                dot_y = Y[fi, stage] + progress * (Y[fi, stage+1] - Y[fi, stage])
            else:
                dot_y = Y[fi, stage]

            dots[fi].set_data([dot_x], [dot_y])

        title_text.set_text(f"Individual fish  ·  Run {stage+1} of 6")
        return lines + dots + [title_text]

    ani = animation.FuncAnimation(fig, update, frames=N,
                                  interval=80, blit=True)
    out = DOCS / "gif_ranks.gif"
    ani.save(out, writer=PillowWriter(fps=12))
    plt.close(fig)
    print(f"  → {out}  ({out.stat().st_size//1024} KB)")


# ══════════════════════════════════════════════════════════════════════════════
# Copy static figures
# ══════════════════════════════════════════════════════════════════════════════
def copy_static():
    print("Copying static figures …")
    for src, dst in [
        (FIGS / "04_reliability_matrix.png", DOCS / "reliability_matrix.png"),
        (FIGS / "run1_minute_curve.png",      DOCS / "run1_curve.png"),
    ]:
        if src.exists():
            shutil.copy(src, dst)
            print(f"  → {dst}")
        else:
            print(f"  MISSING: {src}")


# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    make_trajectory_gif()
    make_habituation_gif()
    make_ranks_gif()
    copy_static()
    print("\nAll done. Files written to docs/images/")
