#!/usr/bin/env python3
"""Figure for the miRanda fourth arm (committed JSONs only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

M = json.load(open("results/miranda_benchmark.json"))
R = json.load(open("results/rnahybrid_benchmark.json"))
DB = json.load(open("results/mirdb_benchmark.json"))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.4))
tools = ["miRDB", "DuplexCNN", "miRanda", "RNAhybrid"]
vals = [DB["pooled_mirdb_auroc"], M["identical_row_three_way"]["pooled_auroc_cnn"],
        M["pooled_auroc_miranda"], R["pooled_auroc_rnahybrid"]]
cols = ["#27ae60", "#2980b9", "#e67e22", "#7f8c8d"]
b = ax1.bar(tools, vals, color=cols)
for r, v in zip(b, vals):
    ax1.text(r.get_x() + r.get_width()/2, v + 0.008, f"{v:.3f}", ha="center", fontsize=8)
ax1.set_ylim(0.45, 0.85); ax1.set_ylabel("pooled AUROC (miRTarBase rows)")
ax1.axhline(0.5, ls=":", c="k", lw=0.8)
ax1.set_title("Four-way ordering, same benchmark", fontsize=9)

pm = {p["mirna"]: p for p in M["per_mirna"]}
x = np.array([pm[k]["auroc_cnn"] for k in pm])
y = np.array([pm[k]["auroc_miranda"] for k in pm])
ax2.scatter(x, y, s=14, c=np.where(y > x, "#e67e22", "#2980b9"), alpha=0.75)
ax2.plot([0.3, 1], [0.3, 1], ls="--", c="k", lw=0.8)
ax2.set_xlabel("DuplexCNN per-miRNA AUROC"); ax2.set_ylabel("miRanda per-miRNA AUROC")
ax2.set_title(f"Per-miRNA: CNN wins {int((x > y).sum())}/{len(pm)}", fontsize=9)
fig.tight_layout()
fig.savefig("paper/figs/fig_miranda.pdf")
print("wrote paper/figs/fig_miranda.pdf")
