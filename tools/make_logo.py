"""Book logo: hex sticker with a mini Manhattan plot in the course role colours. Output: book/images/logo.png."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.path import Path
from matplotlib.patches import PathPatch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "book", "images")
NAVY, PURP, TEAL, AMBER, GREEN, RED, DARK = "#3B5B92", "#6B3FA0", "#1F6F80", "#A85510", "#2E7D5B", "#B03A2E", "#1E2A44"
ROLE = [NAVY, PURP, TEAL, AMBER, GREEN, RED]
plt.rcParams.update({"font.family": "DejaVu Sans"})


def hexagon(r):
    a = np.deg2rad(np.arange(90, 450, 60))
    return np.c_[r * np.cos(a), r * np.sin(a)]


def frame(ax):
    ax.add_patch(Polygon(hexagon(1.0), closed=True, facecolor=DARK, edgecolor="none", zorder=1))
    ax.add_patch(Polygon(hexagon(0.93), closed=True, facecolor="white", edgecolor="none", zorder=2))
    ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.02); ax.set_aspect("equal"); ax.axis("off")
    return Polygon(hexagon(0.93), closed=True, transform=ax.transData)


def text(ax):
    ax.text(0, -0.43, "Genomics & Applied\nBioinformatics", ha="center", va="center", fontsize=16.5,
            fontweight="bold", color=DARK, zorder=6, linespacing=1.15)


# ---------- A: mini Manhattan plot in role colours ----------
rng = np.random.default_rng(4)
fig = plt.figure(figsize=(4, 4)); ax = fig.add_axes([0, 0, 1, 1])
clip = frame(ax)
x0, x1, base = -0.62, 0.62, -0.10
xs = np.linspace(x0, x1, 600)
seg = np.array_split(np.arange(xs.size), 6)
peaks = [0.30, 0.52, 0.22, 0.62, 0.34, 0.44]
for k, idx in enumerate(seg):
    xx = xs[idx]
    y = rng.exponential(0.05, xx.size)
    c = xx[len(xx) // 2]
    y += peaks[k] * np.exp(-((xx - c) / 0.012) ** 2) * rng.uniform(0.7, 1.0, xx.size)
    ax.scatter(xx, base + y, s=5, color=ROLE[k], lw=0, zorder=4, clip_path=clip)
ax.plot([x0, x1], [base + 0.36, base + 0.36], color=RED, lw=1.6, ls=(0, (3, 2)), zorder=5)
ax.plot([x0, x1], [base, base], color=DARK, lw=1.4, zorder=5)
text(ax)
fig.savefig(os.path.join(OUT, "logo.png"), dpi=200, transparent=True); plt.close(fig)
print("ok")
