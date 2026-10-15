"""Forest plot of celltype-level bootstrap estimates."""

import matplotlib.pyplot as plt
import numpy as np


def plot_bootstrap_forest(bootstrap_df, save_path=None):
    """Draw bias-corrected estimates with 95% CI error bars."""
    df = bootstrap_df.drop(index="Intercept", errors="ignore")
    df = df.iloc[::-1]
    y = np.arange(len(df))
    est = df["bias_corrected_estimate"].to_numpy()
    lo = df["CI_lo"].to_numpy()
    hi = df["CI_hi"].to_numpy()

    fig, ax = plt.subplots(figsize=(6, max(3, 0.4 * len(df))))
    ax.errorbar(est, y, xerr=[est - lo, hi - est], fmt="o", color="#d57e7e",
                ecolor="#6d8299", capsize=3)
    ax.axvline(0, color="#444444", linestyle="--", alpha=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(df.index)
    ax.set_xlabel("bias_corrected_estimate (95% CI)")
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150)
    return fig
