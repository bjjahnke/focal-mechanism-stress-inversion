"""Lower-hemisphere stereonet plots of the results (equal-angle projection)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .geometry import plane_vectors, to_trend_plunge  # noqa: E402

SIGMA_COLORS = {"sigma1": "#c0392b", "sigma2": "#2e8b57", "sigma3": "#2c6fbb"}


def project(trend, plunge):
    """Equal-angle (stereographic) x, y for lines given as trend/plunge in degrees."""
    r = np.tan(np.radians(90 - np.asarray(plunge, float)) / 2)
    t = np.radians(np.asarray(trend, float))
    return r * np.sin(t), r * np.cos(t)


def draw_stereonet(ax):
    ax.add_patch(plt.Circle((0, 0), 1, facecolor="white", edgecolor="black", lw=1))
    for p in range(10, 90, 10):
        ax.add_patch(plt.Circle((0, 0), np.tan(np.radians(90 - p) / 2), fill=False, color="0.88", lw=0.6))
    for t in range(0, 180, 30):
        x, y = np.sin(np.radians(t)), np.cos(np.radians(t))
        ax.plot([-x, x], [-y, y], color="0.88", lw=0.6)
    for label, (x, y) in {"N": (0, 1.08), "E": (1.08, 0), "S": (0, -1.12), "W": (-1.12, 0)}.items():
        ax.text(x, y, label, ha="center", va="center", fontsize=9)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect("equal")
    ax.axis("off")


def _poles(strike, dip, rake):
    normals, _ = plane_vectors(strike, dip, rake)
    return to_trend_plunge(normals)


def plot_stress_solutions(solutions: pd.DataFrame, solution_planes: pd.DataFrame, path, max_panels=4):
    top = solutions.head(max_panels)
    fig, axes = plt.subplots(1, len(top), figsize=(4.2 * len(top), 4.8), squeeze=False)
    for ax, sol in zip(axes[0], top.itertuples(index=False)):
        draw_stereonet(ax)
        planes = solution_planes[solution_planes["solution_id"] == sol.solution_id]
        x, y = project(*_poles(planes["strike"], planes["dip"], planes["rake"]))
        ax.plot(x, y, "o", ms=4, color="0.35", label="Poles of chosen planes")
        for key, color in SIGMA_COLORS.items():
            sx, sy = project(getattr(sol, f"{key}_trend"), getattr(sol, f"{key}_plunge"))
            ax.plot(sx, sy, "s", ms=10, color=color, mec="black", label=f"σ{key[-1]}")
        ax.set_title(
            f"Solution {sol.solution_id}: {sol.fraction_of_runs:.0%} of runs\n"
            f"R = {sol.shape_ratio:.2f}, misfit = {sol.misfit_mean_deg:.0f}° ± {sol.misfit_std_deg:.0f}°",
            fontsize=10,
        )
    axes[0][0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_nodal_plane_frequency(frequency: pd.DataFrame, path):
    fig, ax = plt.subplots(figsize=(5.2, 5))
    draw_stereonet(ax)
    x, y = project(*_poles(frequency["strike"], frequency["dip"], frequency["rake"]))
    sc = ax.scatter(x, y, c=frequency["fraction_chosen"], cmap="Greys", vmin=0, vmax=1,
                    edgecolors="black", linewidths=0.4, s=36, zorder=3)
    fig.colorbar(sc, ax=ax, shrink=0.75, label="Fraction of runs the plane was chosen")
    ax.set_title("Poles of both nodal planes for every event", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def make_plots(summary, output_dir) -> list[Path]:
    plot_dir = Path(output_dir) / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    paths = [plot_dir / "stress_solutions.png", plot_dir / "nodal_plane_frequency.png"]
    plot_stress_solutions(summary.stress_solutions, summary.solution_planes, paths[0])
    plot_nodal_plane_frequency(summary.nodal_plane_frequency, paths[1])
    return paths
