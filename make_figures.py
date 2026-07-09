#!/usr/bin/env python3
"""Generate the three figures used in the paper, all as vector PDFs.
Figures 2 and 3 are built directly from /home/MSKHOKHAR/validate/report.json,
"""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
})

with open("/home/MSKHOKHAR/validate/report.json") as f:
    R = json.load(f)

# ----------------------------------------------------------------
# Figure 1: qubit encoding / oracle architecture schematic
# ----------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.9, 3.0))
ax.set_xlim(0, 10)
ax.set_ylim(0, 5)
ax.axis("off")

codon_color = "#3A6EA5"
anc_color = "#B85C38"
op_color = "#4B4B4B"

# input codon qubit lines q0..q5
for i in range(6):
    y = 4.6 - i * 0.5
    ax.plot([0.4, 3.0], [y, y], color=codon_color, lw=1.1, zorder=1)
    ax.text(0.05, y, f"$q_{{{i}}}$", fontsize=8, va="center", ha="left")
ax.text(0.4, 4.95, "codon register (6 qubits, 2 bits/base)", fontsize=7.3,
        color=codon_color, ha="left", va="bottom", style="italic")

# ancilla / output qubit lines q6..q10
for i in range(5):
    y = 1.9 - i * 0.5
    ax.plot([0.4, 3.0], [y, y], color=anc_color, lw=1.1, zorder=1)
    ax.text(0.05, y, f"$q_{{{6+i}}}$", fontsize=8, va="center", ha="left")
ax.text(0.4, -0.15, r"output/ancilla register (5 qubits, $|00000\rangle$ init.)",
        fontsize=7.3, color=anc_color, ha="left", va="top", style="italic")

# the U box
box = FancyBboxPatch((3.15, -0.1), 1.0, 4.9, boxstyle="round,pad=0.02",
                      linewidth=1.2, edgecolor=op_color, facecolor="#EDEDED", zorder=2)
ax.add_patch(box)
ax.text(3.65, 2.35, r"$U$", fontsize=15, ha="center", va="center", zorder=3)

# output lines continue
for i in range(11):
    y = 4.6 - i * 0.5
    ax.plot([4.15, 5.7], [y, y], color=(codon_color if i < 6 else anc_color), lw=1.1, zorder=1)

# measurement symbols on the 5 output qubits
for i in range(5):
    y = 1.9 - i * 0.5
    mbox = Rectangle((5.7, y - 0.18), 0.62, 0.36, linewidth=1.0,
                      edgecolor=op_color, facecolor="white", zorder=2)
    ax.add_patch(mbox)
    arc = mpatches.Arc((6.01, y - 0.03), 0.42, 0.32, angle=0, theta1=20, theta2=160, lw=1.0)
    ax.add_patch(arc)
    ax.annotate("", xy=(6.19, y + 0.1), xytext=(6.01, y - 0.05),
                arrowprops=dict(arrowstyle="->", lw=0.9))
    ax.plot([6.5, 7.15], [y, y], color=anc_color, lw=1.1, zorder=1)
ax.text(7.25, 1.9, r"$\longrightarrow \; y \in \{0,\dots,20\}$" "\n" r"$\rightarrow$ amino acid",
        fontsize=8, va="center", ha="left")

# dashed pass-through for codon qubits (unused after U for readout, shown faded)
for i in range(6):
    y = 4.6 - i * 0.5
    ax.plot([5.7, 7.15], [y, y], color=codon_color, lw=0.8, ls=(0, (2, 2)), alpha=0.5, zorder=1)

ax.text(1.7, -0.75,
        r"Input: $\,|x\rangle_{6}\,|00000\rangle_{5}$, $\;x=$codon$\to\{0,\ldots,63\}$   "
        r"$\quad\quad$ Output (intended): $\,|x\rangle_6\,|y\rangle_5,\;\; y=\mathrm{table}(x)$",
        fontsize=7.6, ha="left", va="top")

plt.tight_layout()
plt.savefig("/home/MSKHOKHAR/paper/figures/fig_architecture.pdf", bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------
# Figure 2: unitarity validation, original vs. corrected construction
# ----------------------------------------------------------------
orig = R["unitarity_original"]
fixed = R["unitarity_fixed"]
labels = [d["label"] for d in orig]
short = [l if len(l) < 18 else l[:16] + "\u2026" for l in labels]
dev_orig = [d["unitarity_frobenius_dev"] for d in orig]
dev_fixed = [max(d["unitarity_frobenius_dev"], 1e-16) for d in fixed]  # avoid log(0)

fig, ax = plt.subplots(figsize=(6.9, 3.1))
x = np.arange(len(labels))
w = 0.38
b1 = ax.bar(x - w/2, dev_orig, width=w, label="Original construction",
            color="#B85C38", edgecolor="black", linewidth=0.4)
b2 = ax.bar(x + w/2, dev_fixed, width=w, label="Corrected construction",
            color="#3A6EA5", edgecolor="black", linewidth=0.4)
ax.set_yscale("log")
ax.set_ylabel(r"$\|U^\dagger U - I\|_F$  (log scale)", fontsize=8.5)
ax.set_xticks(x)
ax.set_xticklabels(short, rotation=55, ha="right", fontsize=6.6)
ax.axhline(1e-10, color="gray", ls="--", lw=0.8)
ax.text(len(labels) - 1, 1.6e-10, "machine-precision threshold", fontsize=6.3,
        ha="right", va="bottom", color="gray")
ax.legend(fontsize=7.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=2)
ax.set_ylim(1e-17, 5e2)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
plt.tight_layout()
plt.savefig("/home/MSKHOKHAR/paper/figures/fig_unitarity.pdf", bbox_inches="tight")
plt.close()

# ----------------------------------------------------------------
# Figure 3: runtime benchmark, classical vs. quantum-simulated
# ----------------------------------------------------------------
timing = R["timing"]
n = [t["n_codons"] for t in timing]
tc = [t["classical_seconds"] for t in timing]
tq = [t["quantum_seconds"] for t in timing]

fig, ax = plt.subplots(figsize=(6.9, 3.0))
ax.plot(n, tc, "o-", color="#3A6EA5", label="Classical dict lookup", lw=1.3, ms=4)
ax.plot(n, tq, "s-", color="#B85C38", label="Quantum-simulated pathway\n(per-codon 11-qubit apply_gate)", lw=1.3, ms=4)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Sequence length (codons)", fontsize=8.5)
ax.set_ylabel("Wall-clock time (s), mean of 3 runs", fontsize=8.5)
ax.legend(fontsize=7.3, frameon=False, loc="upper left")
ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.6)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

# annotate slowdown factor at largest n
last = timing[-1]
ax.annotate(f"{last['slowdown_factor']:.0f}\u00d7 slower\nat n={last['n_codons']}",
            xy=(last["n_codons"], last["quantum_seconds"]),
            xytext=(last["n_codons"] * 0.12, last["quantum_seconds"] * 2.2),
            fontsize=7, ha="center",
            arrowprops=dict(arrowstyle="-", lw=0.6, color="gray"))
plt.tight_layout()
plt.savefig("/home/MSKHOKHAR/paper/figures/fig_benchmark.pdf", bbox_inches="tight")
plt.close()

print("Figures written:")
import os
for f in sorted(os.listdir("/home/MSKHOKHAR/paper/figures")):
    if f.endswith(".pdf"):
        print(" ", f)
