"""Figures for docs/briefing/BRIEFING.md. Reads only committed result files (no model is re-run).

Inputs: results/test/test_results.json, results/baselines/baselines_dev.md, results/ground_truth/GT_REPORT.md,
results/ground_truth/harmony_bins.tsv, mvp/model/kinases.tsv, mvp/model/edges.tsv.gz (built locally; optional),
docs/briefing/dev_trajectory.tsv. Output: docs/briefing/figures/fig*.png
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/briefing/figures"
OUT.mkdir(parents=True, exist_ok=True)

# reference palette (dataviz skill): slots 1-3 validated all-pairs; gray for context series
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GRAY, LIGHT, INK, INK2, GRID = "#9b9a95", "#d6d5d0", "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.dpi": 200, "legend.frameon": False,
})

T = json.loads((ROOT / "results/test/test_results.json").read_text())
M = T["methods"]
SHORT = {"core-v3": "core-v3 (final)", "core-v1": "core-v1", "zscore curated_union min1": "z-score",
         "KSEA curated_union min1": "KSEA", "PhosX-like cantley site-top15 min3": "PhosX-like",
         "JUMPsem curated_union min3 (aligned)": "JUMPsem*", "zscore psp_all min3": "z-score (PSP only)",
         "prior-only": "prior-only control", "random": "random"}


def save(fig, name):
    fig.savefig(OUT / name, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / name)


def md_table(text: str, header_start: str) -> pd.DataFrame:
    lines = text.splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith(header_start))
    rows = []
    for l in lines[i:]:
        if not l.startswith("|"):
            break
        rows.append([c.strip() for c in l.strip().strip("|").split("|")])
    return pd.DataFrame(rows[2:], columns=rows[0])


# ------------------------------------------------------------------ Fig 1: the agentic workflow
def fig1():
    fig, ax = plt.subplots(figsize=(7.4, 3.6))
    ax.set_axis_off(); ax.set_xlim(-0.35, 10.35); ax.set_ylim(0, 5)
    steps = [("Brief\n(human)", "PROMPT.md"), ("Literature\nreview", "LITERATURE.md"),
             ("Design +\n3 reviews", "DESIGN.md"), ("Ground truth\n+ audit", "results/ground_truth"),
             ("Splits +\nharness", "G1 freeze"), ("Baselines", "14 methods"),
             ("Method\nbranches", "B1-B8"), ("Sealed test\n(once)", "G2 freeze"),
             ("Reproduce +\nfact-check", "independent agents"), ("Report", "REPORT.md")]
    xs = np.linspace(0.55, 9.45, 5)
    for i, (t, sub) in enumerate(steps):
        r, c = divmod(i, 5)
        x = xs[c] if r == 0 else xs[4 - c]
        y = 3.6 if r == 0 else 1.4
        col = ORANGE if i in (0,) else (BLUE if i in (4, 7) else "#f2f1ed")
        tc = "white" if i in (0, 4, 7) else INK
        ax.add_patch(FancyBboxPatch((x - 0.78, y - 0.55), 1.56, 1.1, boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc=col, ec=INK2 if col == "#f2f1ed" else col, lw=0.8))
        ax.text(x, y + 0.12, t, ha="center", va="center", fontsize=8.2, color=tc, weight="bold")
        ax.text(x, y - 0.38, sub, ha="center", va="center", fontsize=6.6, color=tc if tc == "white" else INK2)
        if i < len(steps) - 1:
            if r == 0 and c < 4:
                ax.annotate("", (xs[c + 1] - 0.8, y), (x + 0.8, y), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.9))
            elif r == 0 and c == 4:
                ax.annotate("", (x, 1.4 + 0.57), (x, y - 0.57), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.9))
            else:
                ax.annotate("", (xs[3 - c] + 0.8, y), (x - 0.8, y), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.9))
    ax.text(5, 4.75, "PI agent delegates to sub-agents, reviews, decides, and keeps the research record",
            ha="center", fontsize=8.5, color=INK)
    ax.text(5, 0.2, "Blue = freeze gates (hashes pinned before results are seen).  Every headline result is checked "
            "by an agent that did not produce it.", ha="center", fontsize=7.2, color=INK2)
    save(fig, "fig1_workflow.png")


# ------------------------------------------------------------------ Fig 2: benchmark data and study bias
def fig2():
    gt = (ROOT / "results/ground_truth/GT_REPORT.md").read_text()
    sec = gt[gt.index("### Per source x split"):]
    tab = md_table(sec, "| gt conditions")
    tab = tab.set_index(tab.columns[0]).astype(int)
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 2.9), gridspec_kw={"width_ratios": [1, 1.25]})
    names = {"H": "Hijazi 2020", "D": "decryptM", "O": "other studies"}
    cols = {"H": BLUE, "D": ORANGE, "O": AQUA}
    splits = ["dev_cv", "devconfirm", "test"]
    lab = ["Development\n(5 CV folds)", "Confirmation", "Sealed test"]
    bottom = np.zeros(3)
    for s in ["H", "D", "O"]:
        v = tab.loc[s, splits].to_numpy()
        a.bar(lab, v, bottom=bottom, color=cols[s], width=0.6, label=names[s], edgecolor="white", linewidth=1.5)
        bottom += v
    a.set_ylabel("ground-truth conditions"); a.set_title("a  Benchmark conditions per split", loc="left", fontsize=9)
    a.legend(fontsize=7.5, loc="upper right"); a.grid(axis="x", visible=False)

    kin = pd.read_csv(ROOT / "mvp/model/kinases.tsv", sep="\t", keep_default_na=False)
    e = ROOT / "mvp/model/edges.tsv.gz"
    if e.exists():
        ed = pd.read_csv(e, sep="\t")
        n = ed[ed["source"] == "cur"].groupby("kinase").size().reindex(kin["kinase"], fill_value=0)
        n = n.sort_values(ascending=False).to_numpy()
        top10 = n[:10].sum() / n.sum()
        b.bar(np.arange(len(n)), np.maximum(n, 0.8), width=1.0, color=BLUE, linewidth=0)
        b.set_yscale("log"); b.set_xlim(-2, len(n) + 2)
        b.set_xlabel("503 human protein kinases, ranked"); b.set_ylabel("curated substrate sites")
        b.text(0.97, 0.95, f"top 10 kinases: {top10:.0%} of all edges\n{(n == 0).sum()} kinases: no curated substrate",
               transform=b.transAxes, ha="right", va="top", fontsize=7.5, color=INK)
        b.grid(axis="x", visible=False)
    b.set_title("b  Study bias in kinase–substrate knowledge", loc="left", fontsize=9)
    fig.tight_layout()
    save(fig, "fig2_data.png")


# ------------------------------------------------------------------ Fig 3: baselines on DEV and the coverage penalty
def fig3():
    txt = (ROOT / "results/baselines/baselines_dev.md").read_text()
    best = md_table(txt[txt.index("## Best config per method"):], "| method")
    full = md_table(txt, "| # | method")
    full["key"] = full["method"] + "|" + full["library"].str.rstrip("*") + "|" + full["cfg"]
    rows = []
    for _, r in best.iterrows():
        lib, cfg = [x.strip() for x in r["best library/cfg"].split("/")]
        f = full[full["key"] == f"{r['method']}|{lib}|{cfg}"].iloc[0]
        rows.append((r["method"], float(f["M1-sel"]), float(f["M1-sel-scored"]), float(f["cov"])))
    d = pd.DataFrame(rows, columns=["method", "m1", "scored", "cov"])
    d = d[d["method"] != "jumpsem"].sort_values("m1")
    pretty = {"zscore": "z-score", "ksea": "KSEA", "phosx_like": "PhosX-like", "kl_enrich": "Kinase Library enrich.",
              "motif_zscore": "motif z-score", "viper": "VIPER", "norm_wmean": "norm. wmean", "ulm": "ULM",
              "wmean": "wmean", "gsea": "GSEA", "rokai": "RoKAI", "mlm": "MLM", "jumpsem_aligned": "JUMPsem*"}
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    y = np.arange(len(d))
    ax.hlines(y, d["m1"], d["scored"], color=LIGHT, lw=2.2, zorder=1)
    ax.scatter(d["scored"], y, s=34, color="white", edgecolor=GRAY, lw=1.4, zorder=2, label="accuracy on kinases it scores")
    ax.scatter(d["m1"], y, s=38, color=BLUE, zorder=3, label="primary metric (unscored targets count as chance)")
    ax.set_yticks(y, [f"{pretty.get(m, m)}  ({c:.0%})" for m, c in zip(d["method"], d["cov"])], fontsize=7.8)
    ax.axvline(0.5, color=INK2, lw=0.8, ls=":"); ax.text(0.5, len(d) - 0.4, " chance", fontsize=7, color=INK2)
    po = 0.687
    ax.axvline(po, color=ORANGE, lw=1.2); ax.text(po - 0.002, len(d) - 0.4, "prior-only\ncontrol 0.687 ", color=INK, fontsize=7, ha="right", va="top")
    ax.set_xlabel("development-set AUROC (M1-sel)"); ax.set_xlim(0.47, 0.72); ax.grid(axis="y", visible=False)
    ax.set_title("Published methods on the development set (best library each; % = target coverage)", loc="left", fontsize=9)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=7.2)
    save(fig, "fig3_baselines_dev.png")


# ------------------------------------------------------------------ Fig 4: development trajectory
def fig4():
    d = pd.read_csv(ROOT / "docs/briefing/dev_trajectory.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(7.4, 3.0))
    x = d["step"].to_numpy(); v = d["dev_m1sel"].to_numpy()
    ax.plot(x, v, color=BLUE, lw=2, marker="o", ms=7, zorder=3)
    for xi, vi in zip(x, v):
        ax.text(xi, vi - 0.016, f"{vi:.3f}", ha="center", fontsize=7.5, color=INK)
    ax.axhline(0.687, color=ORANGE, lw=1.2); ax.text(-0.3, 0.692, "prior-only control (0.687)", fontsize=7, color=INK)
    ax.axhline(0.5, color=INK2, lw=0.8, ls=":"); ax.text(-0.3, 0.503, "chance", fontsize=7, color=INK2)
    ax.set_xticks(x, [s.replace("\\n", "\n") for s in d["label"]], fontsize=7.2); ax.set_xlim(-0.5, 5.5)
    ax.set_ylim(0.48, 0.72); ax.set_ylabel("development AUROC (M1-sel)"); ax.grid(axis="x", visible=False)
    ax.set_title("Development: what each adopted component added (cross-validated on development data)", loc="left", fontsize=9)
    save(fig, "fig4_dev_trajectory.png")


# ------------------------------------------------------------------ Fig 5: sealed TEST, primary result
def fig5():
    order = ["core-v3", "prior-only", "core-v1", "zscore curated_union min1", "KSEA curated_union min1",
             "PhosX-like cantley site-top15 min3", "zscore psp_all min3", "JUMPsem curated_union min3 (aligned)", "random"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.6, 3.2), gridspec_kw={"width_ratios": [1.05, 1]})
    v = [M[m]["M1-sel"] for m in order]
    col = [BLUE if m == "core-v3" else ORANGE if m == "prior-only" else GRAY for m in order]
    y = np.arange(len(order))[::-1]
    a.barh(y, np.array(v) - 0.4, left=0.4, color=col, height=0.62)
    for yi, vi in zip(y, v):
        a.text(vi + 0.005, yi, f"{vi:.3f}", va="center", fontsize=7.3, color=INK)
    a.set_yticks(y, [SHORT[m] for m in order], fontsize=7.6)
    a.axvline(0.5, color=INK2, lw=0.8, ls=":"); a.set_xlim(0.4, 0.78); a.grid(axis="y", visible=False)
    a.set_xlabel("TEST AUROC (M1-sel)"); a.set_title("a  Sealed test set: primary metric", loc="left", fontsize=9)

    H = T["holm_family"]
    lab = {"zscore curated_union min1": "vs z-score", "KSEA curated_union min1": "vs KSEA",
           "JUMPsem curated_union min3 (aligned)": "vs JUMPsem*", "PhosX-like cantley site-top15 min3": "vs PhosX-like",
           "prior-only": "vs prior-only"}
    rows = [h for h in H if h["endpoint"] == "M1-sel" and h["comparator"] in lab]
    rows += [h for h in H if h["endpoint"] == "M2"]
    yy = np.arange(len(rows))[::-1]
    for yi, h in zip(yy, rows):
        sig = h["p_holm"] < 0.05
        c = BLUE if sig else GRAY
        b.hlines(yi, h["ci_low"], h["ci_high"], color=c, lw=2)
        b.plot(h["estimate"], yi, "o", color=c, ms=6)
        b.text(0.335, yi, f"p={h['p_holm']:.3f}", va="center", fontsize=7, color=INK2)
    b.set_yticks(yy, [lab.get(h["comparator"], "vs z-score") + (" (M2)" if h["endpoint"] == "M2" else "") for h in rows],
                 fontsize=7.6)
    b.axvline(0, color=INK2, lw=0.9); b.set_xlim(-0.06, 0.38); b.grid(axis="y", visible=False)
    b.set_xlabel("difference for core-v3 (95% CI)")
    b.set_title("b  Pre-registered tests (Holm-adjusted p)", loc="left", fontsize=9)
    fig.tight_layout()
    save(fig, "fig5_test_primary.png")


# ------------------------------------------------------------------ Fig 6: who is it good for? strata
def fig6():
    ms = [("core-v3", BLUE), ("prior-only", ORANGE), ("zscore curated_union min1", AQUA)]
    groups = [("0", "0"), ("1-4", "1–4"), ("5-19", "5–19"), (">=20", "≥20")]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.6, 2.9), gridspec_kw={"width_ratios": [1.5, 1]})
    w = 0.26
    for i, (m, c) in enumerate(ms):
        v = [M[m]["harmony_per_bin"][g] for g, _ in groups]
        a.bar(np.arange(4) + (i - 1) * w, v, width=w - 0.02, color=c, label=SHORT[m])
    a.set_xticks(range(4), [l for _, l in groups]); a.set_xlabel("known (PSP) substrates measured per kinase")
    a.set_ylim(0.4, 0.92); a.axhline(0.5, color=INK2, lw=0.8, ls=":"); a.set_ylabel("TEST AUROC")
    a.set_title("a  By how well-studied the kinase is", loc="left", fontsize=9); a.grid(axis="x", visible=False)
    a.legend(fontsize=7, loc="upper left", ncol=1)
    keys = [("T-seen", "seen in\ndevelopment"), ("T-unseen", "never seen in\ndevelopment")]
    for i, (m, c) in enumerate(ms):
        v = [M[m][k] for k, _ in keys]
        b.bar(np.arange(2) + (i - 1) * w, v, width=w - 0.02, color=c)
    b.set_xticks(range(2), [l for _, l in keys]); b.set_ylim(0.2, 0.92); b.axhline(0.5, color=INK2, lw=0.8, ls=":")
    b.set_title("b  Seen vs unseen target kinases", loc="left", fontsize=9); b.grid(axis="x", visible=False)
    fig.tight_layout()
    save(fig, "fig6_test_strata.png")


# ------------------------------------------------------------------ Fig 7: accuracy vs runtime (DEV)
def fig7():
    txt = (ROOT / "results/baselines/baselines_dev.md").read_text()
    full = md_table(txt, "| # | method")
    full["m1"] = full["M1-sel"].astype(float)
    full["ms"] = pd.to_numeric(full["ms/cond 1t"], errors="coerce")
    best = full.sort_values("m1", ascending=False).drop_duplicates("method")
    best = best[best["method"] != "jumpsem"]
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    ax.scatter(best["ms"], best["m1"], s=30, color=GRAY, zorder=2)
    for _, r in best.iterrows():
        nm = {"jumpsem_aligned": "JUMPsem*", "zscore": "z-score", "ksea": "KSEA", "phosx_like": "PhosX-like",
              "rokai": "RoKAI", "gsea": "GSEA", "viper": "VIPER"}.get(r["method"])
        off = {"KSEA": (0.9, -0.011, "right"), "z-score": (1.15, 0.004, "left"), "RoKAI": (0.9, -0.010, "right"),
               "GSEA": (1.15, 0.002, "left")}.get(nm, (1.15, -0.002, "left"))
        if nm:
            ax.text(r["ms"] * off[0], r["m1"] + off[1], nm, fontsize=7, color=INK2, ha=off[2])
    rt = T["runtime"]["cpu"]["threads_1"]["ms_per_condition"]
    ax.scatter([rt], [0.680], s=60, color=BLUE, zorder=3)
    ax.text(rt * 1.2, 0.680, "core-v3 (dev 0.680)", fontsize=7.5, color=INK, va="center")
    ax.axhline(0.687, color=ORANGE, lw=1.1); ax.text(0.6, 0.690, "prior-only control", fontsize=7, color=INK)
    ax.set_xscale("log"); ax.set_xlabel("milliseconds per condition (1 CPU thread, log scale)")
    ax.set_ylabel("development AUROC (M1-sel)"); ax.set_ylim(0.49, 0.71)
    ax.set_title("Accuracy versus speed (development set)", loc="left", fontsize=9)
    save(fig, "fig7_speed.png")


# ------------------------------------------------------------------ Fig 8: external sets
def fig8():
    R = T["external"]["results"]
    ms = [("core-v3", BLUE), ("core-v1", YELLOW), ("zscore curated_union min1", AQUA), ("KSEA curated_union min1", GRAY)]
    sets = [("tyrosine", "Tyrosine-kinase\nperturbations (16 kin.)"), ("hernandez_new", "Hernández new\nstudies (7 kin.)"),
            ("tumour", "CPTAC tumours\n(367 kin.)")]

    def val(m, s):
        r = R[m][s]["M1"]
        return r["M1-kinase"] if s == "tumour" else r["M1-sel"]["value"]
    fig, ax = plt.subplots(figsize=(6.8, 2.9))
    w = 0.2
    for i, (m, c) in enumerate(ms):
        v = [val(m, s) for s, _ in sets]
        ax.bar(np.arange(3) + (i - 1.5) * w, v, width=w - 0.02, color=c, label=SHORT.get(m, m))
    ax.set_xticks(range(3), [l for _, l in sets], fontsize=7.8); ax.set_ylim(0.4, 0.8)
    ax.axhline(0.5, color=INK2, lw=0.8, ls=":"); ax.set_ylabel("AUROC (coverage-aware)")
    ax.legend(fontsize=7, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.13)); ax.grid(axis="x", visible=False)
    save(fig, "fig8_external.png")


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8):
        f()
