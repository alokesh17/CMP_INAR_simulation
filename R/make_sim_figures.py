"""
Generate the two simulation-illustration figures for cmp_inar_draft.tex,
Section 4.2 (parameter recovery) and Section 4.4 (model comparison).

Numbers are transcribed directly from Table~\\ref{tab:sim_placeholder} and
Table~\\ref{tab:model_comparison_winrate} already in the draft -- no new
simulation is run here, this only visualizes existing results.

Palette: dataviz skill's validated default categorical slots 1-4
(blue/orange/aqua/yellow), confirmed via scripts/validate_palette.js
(all-pairs, light mode) before use.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#52514e",
    "axes.labelcolor": "#0b0b0b",
    "xtick.color": "#52514e",
    "ytick.color": "#52514e",
    "axes.grid": True,
    "grid.color": "#e6e5e0",
    "grid.linewidth": 0.6,
})

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"

n_vals = [100, 200, 400, 600]

# Table tab:sim_placeholder, transcribed exactly.
# regime -> parameter -> (true_value, [mean at n=100,200,400,600], [sd ...])
data = {
    "Overdispersion ($\\nu=0.5$)": {
        "color": BLUE,
        "alpha": (0.30, [0.291, 0.312, 0.300, 0.301], [0.063, 0.043, 0.034, 0.022]),
        "lambda": (1.56, [2.435, 2.052, 1.943, 1.719], [0.803, 0.770, 0.480, 0.314]),
        "nu": (0.50, [0.698, 0.620, 0.605, 0.548], [0.205, 0.212, 0.148, 0.116]),
        "rho": (0.30, [0.313, 0.296, 0.314, 0.297], [0.075, 0.060, 0.063, 0.043]),
    },
    "Equidispersion ($\\nu=1.0$)": {
        "color": ORANGE,
        "alpha": (0.30, [0.323, 0.283, 0.295, 0.297], [0.067, 0.053, 0.036, 0.024]),
        "lambda": (3.00, [5.077, 4.314, 3.410, 3.480], [1.924, 1.474, 1.057, 1.114]),
        "nu": (1.00, [1.236, 1.150, 1.030, 1.056], [0.299, 0.240, 0.191, 0.208]),
        "rho": (0.30, [0.340, 0.307, 0.296, 0.293], [0.078, 0.059, 0.055, 0.039]),
    },
    "Underdispersion ($\\nu=1.5$)": {
        "color": AQUA,
        "alpha": (0.30, [0.313, 0.292, 0.307, 0.306], [0.070, 0.038, 0.032, 0.029]),
        "lambda": (5.66, [5.957, 6.670, 7.370, 6.430], [2.327, 2.315, 2.658, 1.607]),
        "nu": (1.50, [1.427, 1.514, 1.617, 1.547], [0.324, 0.210, 0.256, 0.187]),
        "rho": (0.30, [0.308, 0.284, 0.305, 0.297], [0.074, 0.049, 0.036, 0.032]),
    },
}

param_labels = {
    "alpha": r"$\hat\alpha$",
    "lambda": r"$\hat\lambda$",
    "nu": r"$\hat\nu$",
    "rho": r"$\hat\rho$",
}
panel_order = ["nu", "alpha", "lambda", "rho"]  # nu first: the parameter of interest

fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.0))
axes = axes.ravel()

jitter = {0: -6, 1: 0, 2: 6}  # horizontal offset in "n units" so error bars don't overlap

for p_idx, pname in enumerate(panel_order):
    ax = axes[p_idx]
    for r_idx, (regime, d) in enumerate(data.items()):
        true_val, means, sds = d[pname]
        color = d["color"]
        x = np.array(n_vals) + jitter[r_idx]
        ax.errorbar(
            x, means, yerr=sds, fmt="o-", color=color, ecolor=color,
            elinewidth=1.2, capsize=3, markersize=5, linewidth=1.4,
            alpha=0.95, label=regime if p_idx == 0 else None,
        )
        # true-value reference line, thin, regime-colored, recessive
        ax.axhline(true_val, color=color, linestyle=(0, (4, 3)), linewidth=1.0, alpha=0.45)
    ax.set_xticks(n_vals)
    ax.set_xlabel("Sample size $n$")
    ax.set_ylabel(param_labels[pname] + " (posterior mean $\\pm$ SD)")
    ax.set_title(param_labels[pname], fontsize=11, color="#0b0b0b")

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("Parameter recovery across dispersion regimes and sample sizes", y=1.00, fontsize=12)
fig.tight_layout(rect=[0, 0.04, 1, 0.97])
fig.savefig("fig_recovery.pdf", bbox_inches="tight")
fig.savefig("fig_recovery.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 2: model-comparison win rate (Table tab:model_comparison_winrate)
# ---------------------------------------------------------------------------
regimes = ["Overdispersion\n($\\nu=0.5$)", "Equidispersion\n($\\nu=1.0$)", "Underdispersion\n($\\nu=1.5$)"]
models = ["ZICMP", "ZIP", "ZINB", "ZIGP"]
model_colors = [BLUE, ORANGE, AQUA, YELLOW]
winrate = np.array([
    [0.70, 0.00, 0.30, 0.00],   # overdispersion
    [0.00, 0.90, 0.10, 0.00],   # equidispersion
    [0.90, 0.10, 0.00, 0.00],   # underdispersion
])

fig2, ax2 = plt.subplots(figsize=(6.5, 3.6))
x = np.arange(len(regimes))
width = 0.19
for m_idx, model in enumerate(models):
    offset = (m_idx - 1.5) * width
    bars = ax2.bar(x + offset, winrate[:, m_idx], width, label=model,
                    color=model_colors[m_idx], edgecolor="white", linewidth=0.6)
    for b, v in zip(bars, winrate[:, m_idx]):
        if v > 0:
            ax2.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.2f}",
                      ha="center", va="bottom", fontsize=7.5, color="#0b0b0b")

ax2.set_xticks(x)
ax2.set_xticklabels(regimes)
ax2.set_ylabel("Fraction of replicates with\nlowest WAIC$_2$ ($R=10$, $n=600$)")
ax2.set_ylim(0, 1.05)
fig2.suptitle("Which model wins, by dispersion regime", fontsize=12, y=1.02)
ax2.legend(loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.16))
fig2.tight_layout(rect=[0, 0, 1, 0.90])
fig2.savefig("fig_modelcomparison.pdf", bbox_inches="tight")
fig2.savefig("fig_modelcomparison.png", dpi=200, bbox_inches="tight")
plt.close(fig2)

# ---------------------------------------------------------------------------
# Figure 3: estimation bias (posterior mean - true value), same data as
# Figure 1 (Table~tab:sim_placeholder), isolating the direction and
# shrinkage-with-n of the nu attenuation referenced in Section 4.2.
# ---------------------------------------------------------------------------
fig3, axes3 = plt.subplots(2, 2, figsize=(7.5, 6.0))
axes3 = axes3.ravel()

for p_idx, pname in enumerate(panel_order):
    ax = axes3[p_idx]
    ax.axhline(0.0, color="#52514e", linestyle="-", linewidth=0.9, alpha=0.6)
    for r_idx, (regime, d) in enumerate(data.items()):
        true_val, means, sds = d[pname]
        color = d["color"]
        bias = np.array(means) - true_val
        x = np.array(n_vals) + jitter[r_idx]
        ax.plot(
            x, bias, "o-", color=color, markersize=5, linewidth=1.4,
            alpha=0.95, label=regime if p_idx == 0 else None,
        )
    ax.set_xticks(n_vals)
    ax.set_xlabel("Sample size $n$")
    ax.set_ylabel("Bias: " + param_labels[pname] + r"$- $ true value")
    ax.set_title(param_labels[pname], fontsize=11, color="#0b0b0b")

handles3, labels3 = axes3[0].get_legend_handles_labels()
fig3.legend(handles3, labels3, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig3.suptitle("Estimation bias (posterior mean $-$ true value) by dispersion regime and $n$", y=1.00, fontsize=12)
fig3.tight_layout(rect=[0, 0.04, 1, 0.97])
fig3.savefig("fig_bias.pdf", bbox_inches="tight")
fig3.savefig("fig_bias.png", dpi=200, bbox_inches="tight")
plt.close(fig3)

print("Wrote fig_recovery.{pdf,png}, fig_modelcomparison.{pdf,png}, and fig_bias.{pdf,png}")
