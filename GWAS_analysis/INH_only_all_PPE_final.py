#!/usr/bin/env python3
"""
Manhattan plot with kinked Y-axis:
  PPE1 + key genes → ALWAYS LABELED
  Top 5 non-overlapping PPE genes (>275) → LABELED
  No label overlap
Input: eg_for_now.xlsx
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# ----------------------------------------------------------------------
# 1. CONFIG
# ----------------------------------------------------------------------
EXCEL_FILE   = "INH_R_higher_only.xlsx"
SHEET_NAME   = "Sheet1"
POS_COL      = "Position"
GENE_COL     = "Gene"
LOGP_COL     = "-log10(FDR_p)"

KINK_VALUE   = 100.0
Y_MAX        = 320.0

OUTPUT_PNG   = "manhattan_PPE_TOP5.png"
OUTPUT_PDF   = "manhattan_PPE_TOP5.pdf"
DPI          = 300
FIG_WIDTH    = 14
FIG_HEIGHT   = 8

TOP_RATIO    = 1.0
BOTTOM_RATIO = 3.0

PAD_LOW  = 2.0
PAD_HIGH = 20.0

# ----------------------------------------------------------------------
# GENE GROUPS
# ----------------------------------------------------------------------
PPE_GENES = {f"PPE{i}" for i in range(1, 70)}
PPE1_GENE = "PPE1"
KEY_GENES = {"PPE1", "katG", "inhA", "fabG1", "ahpC", "furA"}

# Colors
COLOR_PPE1      = "#d62728"  # RED
COLOR_KEY_OTHER = "#2ca02c"  # GREEN
COLOR_PPE_OTHER = "#1f77b4"  # BLUE
COLOR_OTHER     = "#bbbbbb"  # GREY

MIN_SIZE = 15
MAX_SIZE = 80
SIZE_SCALE_FACTOR = 0.8

# Label spacing: minimum distance (in Mb) between labels
MIN_LABEL_DISTANCE_MB = 0.3  # ~300 kb

# ----------------------------------------------------------------------
# 2. LOAD DATA
# ----------------------------------------------------------------------
df = pd.read_excel(EXCEL_FILE, sheet_name=SHEET_NAME)

required = [POS_COL, GENE_COL, LOGP_COL]
for c in required:
    if c not in df.columns:
        raise ValueError(f"Column '{c}' missing")

df["neglogp"] = df[LOGP_COL]
df["pos_mb"] = df[POS_COL] / 1e6

# ----------------------------------------------------------------------
# 3. ONE DOT PER GENE
# ----------------------------------------------------------------------
gene_stats = df.groupby(GENE_COL).agg(
    neglogp_max=("neglogp", "max"),
    pos_mb=("pos_mb", "first"),
    mutation_count=("neglogp", "count")
).reset_index()

df_max = df.loc[df.groupby(GENE_COL)["neglogp"].idxmax()].copy()
df_max = df_max.merge(gene_stats[[GENE_COL, "mutation_count"]], on=GENE_COL)

print(f"Total genes: {len(df_max)}")
print(f"Max -log10(p): {df_max['neglogp'].max():.1f}")

# ----------------------------------------------------------------------
# 4. COLOR & SIZE
# ----------------------------------------------------------------------
def assign_color(gene):
    if gene == PPE1_GENE:
        return COLOR_PPE1
    elif gene in KEY_GENES:
        return COLOR_KEY_OTHER
    elif gene in PPE_GENES:
        return COLOR_PPE_OTHER
    else:
        return COLOR_OTHER

df_max["color"] = df_max[GENE_COL].apply(assign_color)

sizes = df_max["mutation_count"]
sizes_scaled = MIN_SIZE + (sizes - sizes.min()) / (sizes.max() - sizes.min() + 1e-6) * (MAX_SIZE - MIN_SIZE)
df_max["size"] = np.clip(sizes_scaled * SIZE_SCALE_FACTOR, MIN_SIZE, MAX_SIZE)

# ----------------------------------------------------------------------
# 5. LABEL LOGIC
# ----------------------------------------------------------------------
df_max["label"] = df_max[GENE_COL].isin(KEY_GENES)

# ----------------------------------------------------------------------
# 6. FIGURE
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
gs = gridspec.GridSpec(2, 1, height_ratios=[TOP_RATIO, BOTTOM_RATIO], hspace=0.08)

ax_top    = fig.add_subplot(gs[0])
ax_bottom = fig.add_subplot(gs[1], sharex=ax_top)

# ----------------------------------------------------------------------
# 7. PLOT ALL NON-KEY GENES FIRST
# ----------------------------------------------------------------------
mask_high = df_max["neglogp"] > KINK_VALUE
mask_low  = df_max["neglogp"] <= KINK_VALUE

non_key_high = mask_high & ~df_max[GENE_COL].isin(KEY_GENES)
non_key_low  = mask_low  & ~df_max[GENE_COL].isin(KEY_GENES)

if non_key_high.any():
    ax_top.scatter(
        df_max.loc[non_key_high, "pos_mb"],
        df_max.loc[non_key_high, "neglogp"],
        c=df_max.loc[non_key_high, "color"],
        s=df_max.loc[non_key_high, "size"],
        edgecolor="none", alpha=0.9, clip_on=False, zorder=5
    )

if non_key_low.any():
    ax_bottom.scatter(
        df_max.loc[non_key_low, "pos_mb"],
        df_max.loc[non_key_low, "neglogp"],
        c=df_max.loc[non_key_low, "color"],
        s=df_max.loc[non_key_low, "size"],
        edgecolor="none", alpha=0.8, clip_on=True, zorder=5
    )

# ----------------------------------------------------------------------
# 8. PLOT KEY GENES (ON TOP)
# ----------------------------------------------------------------------
key_rows = df_max[df_max[GENE_COL].isin(KEY_GENES)]
for _, row in key_rows.iterrows():
    ax = ax_top if row["neglogp"] > KINK_VALUE else ax_bottom
    ax.scatter(row["pos_mb"], row["neglogp"], c=row["color"], s=row["size"],
               edgecolor="black", linewidth=1.2, alpha=1.0, clip_on=False, zorder=100)
    ax.annotate(row[GENE_COL], (row["pos_mb"], row["neglogp"]),
                xytext=(0, 15), textcoords='offset points',
                fontsize=10, fontweight='bold', color=row["color"],
                ha='center', va='bottom',
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="black", linewidth=1.0, alpha=0.9),
                zorder=101)

# ----------------------------------------------------------------------
# 9. SELECT TOP 5 NON-OVERLAPPING PPE GENES (>275)
# ----------------------------------------------------------------------
high_ppe = df_max[
    (df_max["neglogp"] > 275) &
    df_max[GENE_COL].isin(PPE_GENES) &
    ~df_max[GENE_COL].isin(KEY_GENES)
].copy()

if not high_ppe.empty:
    high_ppe = high_ppe.sort_values("neglogp", ascending=False)
    selected = []
    positions = []

    for _, row in high_ppe.iterrows():
        pos = row["pos_mb"]
        if not positions or min(abs(pos - p) for p in positions) >= MIN_LABEL_DISTANCE_MB:
            selected.append(row)
            positions.append(pos)
        if len(selected) >= 5:
            break

    # Plot selected PPEs
    for row in selected:
        ax = ax_top
        ax.scatter(row["pos_mb"], row["neglogp"], c=COLOR_PPE_OTHER, s=row["size"],
                   edgecolor="black", linewidth=1.0, alpha=1.0, clip_on=False, zorder=90)
        ax.annotate(row[GENE_COL], (row["pos_mb"], row["neglogp"]),
                    xytext=(0, 15), textcoords='offset points',
                    fontsize=10, fontweight='bold', color=COLOR_PPE_OTHER,
                    ha='center', va='bottom',
                    bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="none", alpha=0.9),
                    zorder=91)

# ----------------------------------------------------------------------
# 10. LABEL OTHER PPEs (≤ 275)
# ----------------------------------------------------------------------
other_ppe_labeled = df_max[
    df_max[GENE_COL].isin(PPE_GENES - KEY_GENES) &
    (df_max["neglogp"] <= 275)
]
def label_gene(ax, row):
    ax.annotate(row[GENE_COL], (row["pos_mb"], row["neglogp"]),
                xytext=(0, 15), textcoords='offset points',
                fontsize=10, fontweight='bold', color=row["color"],
                ha='center', va='bottom',
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="none", alpha=0.9),
                zorder=10)

for _, row in other_ppe_labeled[other_ppe_labeled["neglogp"] > KINK_VALUE].iterrows():
    label_gene(ax_top, row)
for _, row in other_ppe_labeled[other_ppe_labeled["neglogp"] <= KINK_VALUE].iterrows():
    label_gene(ax_bottom, row)

# ----------------------------------------------------------------------
# 11. AXES
# ----------------------------------------------------------------------
ax_top.set_ylim(KINK_VALUE - PAD_HIGH, Y_MAX)
ax_top.set_yticks([100, 150, 200, 250, 300])
ax_top.spines["bottom"].set_visible(False)
ax_top.tick_params(axis='y', length=4, labelsize=10)
ax_top.tick_params(axis='x', bottom=False, labelbottom=False)

ax_bottom.set_ylim(-PAD_LOW, KINK_VALUE + PAD_LOW)
ax_bottom.set_yticks(np.arange(0, 101, 10))
ax_bottom.spines["top"].set_visible(False)
ax_bottom.tick_params(axis='y', length=4, labelsize=10)

# ----------------------------------------------------------------------
# 12. KINK SYMBOLS
# ----------------------------------------------------------------------
d = 0.015
kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False, lw=1.5)
ax_top.plot((-d, +d), (-d, +d), **kwargs)
ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)

kwargs = dict(transform=ax_bottom.transAxes, color='k', clip_on=False, lw=1.5)
ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)
ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

# ----------------------------------------------------------------------
# 13. LABELS & TITLE
# ----------------------------------------------------------------------
fig.text(0.03, 0.5, r'$-\log_{10}(P_{\text{FDR}})$',
         va='center', ha='center', rotation='vertical', fontsize=13, fontweight='medium')

ax_bottom.set_xlabel("Genomic position (Mb)", fontsize=12)
fig.suptitle("Manhattan Plot – PPE1 (Red) + Key Genes (Green) + Top 5 PPE (>275, non-overlap)", 
             fontsize=14, y=0.96)

# ----------------------------------------------------------------------
# 14. SAVE
# ----------------------------------------------------------------------
plt.tight_layout(rect=[0.05, 0.02, 1, 0.94])
plt.savefig(OUTPUT_PNG, dpi=DPI, bbox_inches='tight')
plt.savefig(OUTPUT_PDF, dpi=DPI, bbox_inches='tight')
print(f"SAVED: {OUTPUT_PNG}  &  {OUTPUT_PDF}")
plt.show()
