"""Simplified tab icon: navy hexagon, four bold role-coloured peaks over a red threshold. Output favicon.png (256 px)."""
import os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "book", "images")
DARK, RED = "#1E2A44", "#B03A2E"
COLS = ["#3B5B92", "#6B3FA0", "#A85510", "#2E7D5B"]
a = np.deg2rad(np.arange(90, 450, 60))
fig = plt.figure(figsize=(2.56, 2.56)); ax = fig.add_axes([0, 0, 1, 1])
ax.add_patch(Polygon(np.c_[np.cos(a), np.sin(a)], facecolor=DARK, edgecolor="none"))
ax.add_patch(Polygon(0.86 * np.c_[np.cos(a), np.sin(a)], facecolor="white", edgecolor="none"))
base = -0.42
for x, h, c in zip([-0.45, -0.15, 0.15, 0.45], [0.55, 0.95, 0.70, 0.40], COLS):
    ax.plot([x, x], [base, base + h], color=c, lw=13, solid_capstyle="round")
ax.plot([-0.62, 0.62], [base + 0.5, base + 0.5], color=RED, lw=4, ls=(0, (2.2, 1.4)))
ax.plot([-0.62, 0.62], [base - 0.08, base - 0.08], color=DARK, lw=5)
ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.02); ax.set_aspect("equal"); ax.axis("off")
fig.savefig(os.path.join(OUT, "favicon.png"), dpi=100, transparent=True)
print("ok")
