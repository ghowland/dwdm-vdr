# ================================================================
# STANDALONE HELPERS — paste at top of script, replacing the
# data_5_diagram_lib import and set_outdir call
# ================================================================

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os

# --- Output directory ---
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures')
os.makedirs(OUTDIR, exist_ok=True)

# --- Colors ---
BG      = '#0a0a12'
PAN     = '#12121f'
GOLD    = '#d4a843'
SILVER  = '#a0a8b8'
CYAN    = '#4ecdc4'
MAG     = '#c74b7a'
BLUE    = '#5b8def'
GREEN   = '#6bcf7f'
RED     = '#e05555'
ORANGE  = '#e8944a'
WHITE   = '#e8e8f0'
DIM     = '#555570'
PURPLE  = '#9b7bd4'

# --- Provenance log ---
_PROV = []
def prov(name, value, source):
    _PROV.append((name, value, source))

# --- Figure creators ---
def _style_ax(ax):
    ax.set_facecolor(BG)
    for spine in ax.spines.values():
        spine.set_color(DIM)
        spine.set_linewidth(0.5)
    ax.tick_params(colors=DIM, labelsize=9)

def dark_fig(title, xlabel, ylabel, size=(16, 10)):
    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor(BG)
    _style_ax(ax)
    ax.set_title(title, color=GOLD, fontsize=15, fontweight='bold', pad=14)
    ax.set_xlabel(xlabel, color=SILVER, fontsize=11)
    ax.set_ylabel(ylabel, color=SILVER, fontsize=11)
    ax.grid(True, color=DIM, alpha=0.15, linewidth=0.5)
    return fig, ax

def dark_canvas(title, size=(16, 10)):
    fig, ax = plt.subplots(figsize=size)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    ax.axis('off')
    ax.set_title(title, color=GOLD, fontsize=15, fontweight='bold', pad=14)
    return fig, ax

def dark_fig_dual(title_l, title_r, size=(18, 9), wspace=0.30):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=size,
                                    gridspec_kw={'wspace': wspace})
    fig.patch.set_facecolor(BG)
    _style_ax(ax1)
    _style_ax(ax2)
    ax1.set_title(title_l, color=GOLD, fontsize=14, fontweight='bold', pad=12)
    ax2.set_title(title_r, color=GOLD, fontsize=14, fontweight='bold', pad=12)
    ax1.grid(True, color=DIM, alpha=0.15, linewidth=0.5)
    ax2.grid(True, color=DIM, alpha=0.15, linewidth=0.5)
    return fig, ax1, ax2

# --- Save ---
def save_fig(fig, filename):
    path = os.path.join(OUTDIR, filename)
    fig.savefig(path, dpi=180, facecolor=BG, bbox_inches='tight', pad_inches=0.3)
    plt.close(fig)
    print("  Saved: %s" % filename)

# --- Data points ---
def data_point(ax, x, y, label, color, size=200):
    ax.scatter([x], [y], s=size, c=color, edgecolors=WHITE,
               linewidth=1.5, zorder=5, label=label)

def data_point_err(ax, x, y, yerr, label, color, size=200):
    ax.errorbar(x, y, yerr=yerr, fmt='none', ecolor=color,
                elinewidth=1.5, capsize=4, capthick=1.5, zorder=4)
    data_point(ax, x, y, label, color, size)

def measured_diamond(ax, x, y, label, color, size=200):
    ax.scatter([x], [y], s=size, c=color, edgecolors=WHITE,
               linewidth=1.5, zorder=5, label=label, marker='D')

# --- Curves ---
def curve(ax, x, y, label, color, width=2.5, style='-', alpha=1.0):
    ax.plot(x, y, color=color, linewidth=width, linestyle=style,
            alpha=alpha, label=label, zorder=3)

# --- Regions ---
def shaded_region(ax, x0, x1, color, alpha=0.06, label=None):
    ax.axvspan(x0, x1, color=color, alpha=alpha, label=label, zorder=1)

def shaded_region_h(ax, y0, y1, color, alpha=0.06, label=None):
    ax.axhspan(y0, y1, color=color, alpha=alpha, label=label, zorder=1)

# --- Bands ---
def measurement_band(ax, value, unc, label, color):
    ax.axhspan(value - unc, value + unc, color=color, alpha=0.15, zorder=1)
    ax.axhspan(value - 3 * unc, value + 3 * unc, color=color, alpha=0.05,
               zorder=0, label=label)
    ax.axhline(y=value, color=color, lw=1.2, ls='--', alpha=0.7, zorder=2)

def measurement_band_v(ax, value, unc, label, color):
    ax.axvspan(value - unc, value + unc, color=color, alpha=0.15, zorder=1)
    ax.axvspan(value - 3 * unc, value + 3 * unc, color=color, alpha=0.05,
               zorder=0, label=label)
    ax.axvline(x=value, color=color, lw=1.2, ls='--', alpha=0.7, zorder=2)

# --- Lines ---
def threshold_line(ax, x, label, color, vertical=True):
    if x is None:
        return
    if vertical:
        ax.axvline(x=x, color=color, lw=1.2, ls='--', alpha=0.7, zorder=2,
                   label=label)
    else:
        ax.axhline(y=x, color=color, lw=1.2, ls='--', alpha=0.7, zorder=2,
                   label=label)

# --- Text ---
def result_box(ax, x, y, text, color, fontsize=10):
    ax.text(x, y, text, color=color, fontsize=fontsize,
            va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.4', facecolor=BG,
                      edgecolor=color, alpha=0.9),
            zorder=6)

def note(ax, x, y, text, color, fontsize=9):
    ax.text(x, y, text, color=color, fontsize=fontsize,
            ha='center', va='center', zorder=6)

def arrow_label(ax, x_data, y_data, x_text, y_text, text, color):
    ax.annotate(text, xy=(x_data, y_data), xytext=(x_text, y_text),
                color=color, fontsize=9, ha='center', va='center',
                arrowprops=dict(arrowstyle='->', color=color, lw=1.2),
                zorder=6)

# --- Bar chart ---
def bar_chart(ax, labels, values, colors, width=0.6):
    x = np.arange(len(labels))
    bars = ax.bar(x, values, width=width, color=colors, edgecolor=colors,
                  alpha=0.75, linewidth=1.5, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, color=SILVER, fontsize=9)
    for i in range(len(bars)):
        ax.text(x[i], values[i] + 0.02 * max(values), "%.4g" % values[i],
                ha='center', va='bottom', color=WHITE, fontsize=9)
    return bars

# --- Concentric shells ---
def concentric_shells(ax, shells, center=(0, 0)):
    for i in range(len(shells)):
        radius, color, label, alpha = shells[i]
        circle = plt.Circle(center, radius, facecolor=color, alpha=alpha,
                            edgecolor=color, linewidth=1.2, zorder=2 + i)
        ax.add_patch(circle)
        ax.text(center[0] + radius + 0.02 * radius, center[1],
                label, color=color, fontsize=9, va='center')

# --- Running curves ---
def running_curves(ax, log_mu, inv_alphas, labels, colors):
    for i in range(len(inv_alphas)):
        curve(ax, log_mu, inv_alphas[i], labels[i], colors[i])

# --- Legend ---
def legend(ax, loc='lower right'):
    leg = ax.legend(loc=loc, facecolor=PAN, edgecolor=DIM,
                    labelcolor=WHITE, fontsize=9)
    if leg:
        leg.set_zorder(7)


# ================================================================
# PHYSICAL CONSTANTS FROM THE PAPER
# ================================================================
prov("beta2", -20.41e-24, "README: -20.41 ps^2/km = -20.41e-24 s^2/m")
prov("gamma_nl", 1.3e-3, "README: 1.3 /W/km = 1.3e-3 /W/m")
prov("alpha_db_per_km", 0.2, "README: 0.2 dB/km")
prov("anchor_freq_thz", 193.1, "README: ITU anchor 193.1 THz")
prov("datacenter_drift", 9.62e-2, "Output: VDR-Float drift at 5 steps")
prov("datacenter_nl_phase", 0.016622, "Output: NL phase per channel")
prov("datacenter_pout_dbm", -0.7997, "Output: per-channel output power")
prov("launch_power_mw", 1.0, "Output: 0 dBm = 1 mW")
prov("q335_digits", 101, "2^335 ~ 10^101")
prov("dc_step_size_m", 400, "2 km / 5 steps = 400 m")

BETA2 = -20.41e-24       # s^2/m
GAMMA = 1.3e-3            # /W/m
DZ_DC = 400.0             # datacenter step size in meters
P_LAUNCH = 1.0e-3         # 1 mW = 0 dBm
ANCHOR = 193.1e12         # Hz

# Datacenter channel frequencies
DC_FREQS_THZ = [193.0625, 193.0875, 193.1125, 193.1375]
DC_OFFSETS_GHZ = [(f - 193.1) * 1000 for f in DC_FREQS_THZ]


# ================================================================
# FIG 1: FWM PHASE MATCHING EFFICIENCY VS CHANNEL SPACING
# Type: Curve
# Shows: Lorentzian rolloff explains why FWM matters at 25 GHz
#        but not at 100 GHz — shape IS the physics
# ================================================================

fig, ax = dark_fig(
    "FWM Phase Matching Efficiency vs Channel Spacing",
    "Channel Spacing (GHz)",
    "Phase Matching Efficiency  \u03b7",
    size=(16, 10)
)

spacings_ghz = np.linspace(5, 200, 500)
spacings_hz = spacings_ghz * 1e9

# Phase mismatch for degenerate FWM: delta_beta = beta2 * (2*pi*delta_f)^2
# Efficiency: eta = 1 / (1 + (delta_beta * dz)^2)
delta_beta = BETA2 * (2.0 * np.pi * spacings_hz) ** 2
eta = 1.0 / (1.0 + (delta_beta * DZ_DC) ** 2)

curve(ax, spacings_ghz, eta, "Lorentzian efficiency", CYAN, width=2.5)

# Mark standard spacings
standard_spacings = [12.5, 25, 50, 100]
spacing_labels = ["12.5 GHz", "25 GHz", "50 GHz", "100 GHz"]
spacing_colors = [RED, ORANGE, GREEN, BLUE]

# Precompute eta at each standard spacing for annotation
for i in range(len(standard_spacings)):
    sp = standard_spacings[i]
    sp_hz = sp * 1e9
    db = BETA2 * (2.0 * np.pi * sp_hz) ** 2
    eta_val = 1.0 / (1.0 + (db * DZ_DC) ** 2)
    data_point(ax, sp, eta_val, spacing_labels[i], spacing_colors[i], size=220)

# Annotations with generous offsets — staggered vertically
annotation_offsets = [
    (12.5, 38, 0.78),
    (25, 58, 0.62),
    (50, 90, 0.46),
    (100, 145, 0.30),
]
for i in range(len(annotation_offsets)):
    sp, xt, yt = annotation_offsets[i]
    sp_hz = sp * 1e9
    db = BETA2 * (2.0 * np.pi * sp_hz) ** 2
    eta_val = 1.0 / (1.0 + (db * DZ_DC) ** 2)
    arrow_label(ax, sp, eta_val, xt, yt,
                "%s\n\u03b7 = %.4f" % (spacing_labels[i], eta_val),
                spacing_colors[i])

# Shaded region: high-efficiency FWM zone
shaded_region(ax, 5, 35, RED, alpha=0.06, label="Strong FWM zone")

ax.set_xlim(-5, 215)
ax.set_ylim(-0.08, 1.12)

result_box(ax, 130, 0.95,
           "dz = 400 m (datacenter)\n\u03b2\u2082 = \u221220.41 ps\u00b2/km",
           SILVER, fontsize=9)

legend(ax, loc="upper right")
save_fig(fig, "dwdm_01_fwm_efficiency.png")


# ================================================================
# FIG 2: FOUR-WAVE MIXING TRIPLET GEOMETRY
# Type: Geometric diagram
# Shows: Three channels generating FWM product at fourth frequency,
#        with actual datacenter channel frequencies
# ================================================================

fig, ax = dark_canvas("Four-Wave Mixing: Frequency Conservation", size=(18, 12))
ax.set_xlim(-1, 17)
ax.set_ylim(-1, 11)
ax.set_aspect("equal")

# Frequency axis at y=2
ax.annotate("", xy=(16, 2), xytext=(0, 2),
            arrowprops=dict(arrowstyle="->", color=DIM, lw=1.5))
note(ax, 16.0, 1.4, "Frequency", DIM, fontsize=10)

# Channel positions mapped to x-coords
# 193.0625, 193.0875, 193.1125, 193.1375 THz
# Offsets from 193.0625: 0, 25, 50, 75 GHz
# Map to x: 2, 5.5, 9, 12.5
ch_x = [2.0, 5.5, 9.0, 12.5]
ch_labels = ["f\u2081\n193.0625", "f\u2082\n193.0875", "f\u2083\n193.1125", "f\u2084\n193.1375"]
ch_colors = [CYAN, GREEN, ORANGE, RED]

for i in range(4):
    # Vertical lines at channel positions
    ax.plot([ch_x[i], ch_x[i]], [2, 4.0], color=ch_colors[i], lw=2.5, solid_capstyle="round")
    data_point(ax, ch_x[i], 4.0, None, ch_colors[i], size=180)
    note(ax, ch_x[i], 0.6, ch_labels[i], ch_colors[i], fontsize=9)

# Show degenerate FWM triplet: f1 + f2 - f3 -> f_product
# Example: f1 + f2 - f1 = f2 (trivial), skip
# Interesting: 2*f1 - f2 = 2*193.0625 - 193.0875 = 193.0375 (outside grid)
# Show: f1 + f3 - f2 = 193.0625 + 193.1125 - 193.0875 = 193.0875 = f2 (lands on channel!)
# Show: f1 + f4 - f3 = 193.0625 + 193.1375 - 193.1125 = 193.0875 = f2

# Triplet 1: f1 + f3 -> f2 (absorb from f2)
# Draw arcs from f1 and f3 converging
y_arc1 = 6.5
note(ax, 5.5, 7.6,
     "Triplet A:  f\u2081 + f\u2083 \u2212 f\u2082 = f\u2082",
     GOLD, fontsize=11)
note(ax, 5.5, 7.0,
     "193.0625 + 193.1125 \u2212 193.0875 = 193.0875 THz",
     SILVER, fontsize=9)

# Arrows from f1 and f3 to midpoint above f2
ax.annotate("", xy=(5.5, y_arc1), xytext=(2.0, 4.3),
            arrowprops=dict(arrowstyle="->, head_width=0.3", color=CYAN,
                            lw=1.8, connectionstyle="arc3,rad=-0.2"))
ax.annotate("", xy=(5.5, y_arc1), xytext=(9.0, 4.3),
            arrowprops=dict(arrowstyle="->, head_width=0.3", color=ORANGE,
                            lw=1.8, connectionstyle="arc3,rad=0.2"))
# Arrow down to f2 position
ax.annotate("", xy=(5.5, 4.3), xytext=(5.5, y_arc1),
            arrowprops=dict(arrowstyle="->, head_width=0.3", color=GOLD,
                            lw=2.0))

# Triplet 2: 2*f1 - f2 -> f_out (outside grid)
y_arc2 = 9.5
f_out = 2 * 193.0625 - 193.0875  # = 193.0375
# Map f_out to x: 193.0375 is 25 GHz below f1, so x = 2.0 - 3.5 = -1.5
# Clamp to visible range
x_out = 0.2

note(ax, 4.5, 10.2,
     "Triplet B:  2f\u2081 \u2212 f\u2082 = f_out",
     MAG, fontsize=11)
note(ax, 4.5, 9.6,
     "2(193.0625) \u2212 193.0875 = 193.0375 THz (outside grid)",
     SILVER, fontsize=9)

ax.annotate("", xy=(1.1, y_arc2 - 1.5), xytext=(2.0, 4.3),
            arrowprops=dict(arrowstyle="->, head_width=0.3", color=CYAN,
                            lw=1.5, connectionstyle="arc3,rad=-0.15"))
ax.annotate("", xy=(1.1, y_arc2 - 1.5), xytext=(5.5, 4.3),
            arrowprops=dict(arrowstyle="->, head_width=0.3", color=GREEN,
                            lw=1.5, linestyle="dashed",
                            connectionstyle="arc3,rad=-0.25"))
# Mark the out-of-grid product
ax.plot([x_out, x_out], [2, 3.5], color=MAG, lw=2.0, linestyle="dashed")
data_point(ax, x_out, 3.5, None, MAG, size=150)
note(ax, x_out, 0.6, "f_out\n193.0375", MAG, fontsize=9)

# Phase matching condition box
result_box(ax, 8.5, 9.5,
           "Phase matching:\n\u0394\u03b2 = \u03b2\u2082(\u03c9\u2081\u00b2 + \u03c9\u2082\u00b2 \u2212 \u03c9\u2083\u00b2 \u2212 \u03c9\u2084\u00b2)\n\u03b7 = 1/(1 + (\u0394\u03b2\u00b7dz)\u00b2)",
           CYAN, fontsize=9)

# Spacing annotation
ax.annotate("", xy=(5.5, 1.6), xytext=(2.0, 1.6),
            arrowprops=dict(arrowstyle="<->", color=DIM, lw=1.2))
note(ax, 3.75, 1.2, "25 GHz", DIM, fontsize=9)

save_fig(fig, "dwdm_02_fwm_triplet.png")


# ================================================================
# FIG 3: MULTI-CHANNEL CROSSTALK MATRIX
# Type: Heatmap
# Shows: Channel-to-channel interaction strength: diagonal = SPM,
#        off-diagonal = XPM (2x), FWM at phase-matched triplets
# ================================================================

fig, ax = dark_canvas("Channel Crosstalk Interaction Matrix  (4-ch Datacenter)",
                      size=(16, 14))
ax.set_xlim(-2.5, 10)
ax.set_ylim(-3.5, 9.5)
ax.set_aspect("equal")

# 4x4 grid, cells at integer positions
# SPM diagonal: gamma * P * L_eff
# XPM off-diagonal: 2 * gamma * P * L_eff (factor of 2)
# FWM: at phase-matched triplets, additional contribution

n_ch = 4
cell_size = 1.8
gap = 0.15
grid_x0 = 1.0
grid_y0 = 0.5

ch_names = ["Ch 0\n193.0625", "Ch 1\n193.0875", "Ch 2\n193.1125", "Ch 3\n193.1375"]

# Interaction strengths (relative to SPM = 1)
# SPM = gamma * P * Leff -> 1x
# XPM = 2 * gamma * P * Leff -> 2x
# FWM adds small perturbation at phase-matched pairs
# Nearest neighbors: strongest XPM + some FWM
# Next-nearest: XPM only, weaker FWM

matrix = np.array([
    [1.0, 2.0, 2.0, 2.0],
    [2.0, 1.0, 2.0, 2.0],
    [2.0, 2.0, 1.0, 2.0],
    [2.0, 2.0, 2.0, 1.0],
])

# FWM bonus for degenerate triplets (nearest-neighbor pairs)
fwm_pairs = [(0, 1), (1, 0), (1, 2), (2, 1), (2, 3), (3, 2)]
for pi, pj in fwm_pairs:
    matrix[pi][pj] += 0.3

# Color mapping
max_val = np.max(matrix)
min_val = np.min(matrix)

for row in range(n_ch):
    for col in range(n_ch):
        x = grid_x0 + col * (cell_size + gap)
        y = grid_y0 + (n_ch - 1 - row) * (cell_size + gap)
        val = matrix[row][col]

        # Color intensity
        frac = (val - min_val) / (max_val - min_val)

        if row == col:
            # SPM: blue tones
            c = BLUE
            alpha = 0.3 + 0.4 * frac
            label = "SPM\n\u03b3P\u00b7L_eff"
            label_val = "1.0\u00d7"
        else:
            is_fwm = (row, col) in fwm_pairs
            if is_fwm:
                c = RED
                alpha = 0.3 + 0.5 * frac
                label = "XPM+FWM\n2\u03b3P + FWM"
                label_val = "%.1f\u00d7" % val
            else:
                c = ORANGE
                alpha = 0.2 + 0.4 * frac
                label = "XPM\n2\u03b3P\u00b7L_eff"
                label_val = "%.1f\u00d7" % val

        rect = FancyBboxPatch(
            (x, y), cell_size, cell_size,
            boxstyle="round,pad=0.1",
            facecolor=c, alpha=alpha,
            edgecolor=WHITE, linewidth=1.2
        )
        ax.add_patch(rect)

        # Value text
        cx = x + cell_size / 2
        cy = y + cell_size / 2 + 0.25
        ax.text(cx, cy, label_val, color=WHITE, fontsize=12,
                fontweight="bold", ha="center", va="center")
        ax.text(cx, cy - 0.45, label.split("\n")[0], color=SILVER, fontsize=8,
                ha="center", va="center")

# Row and column labels
for i in range(n_ch):
    x = grid_x0 + i * (cell_size + gap) + cell_size / 2
    y_top = grid_y0 + n_ch * (cell_size + gap) + 0.15
    ax.text(x, y_top, ch_names[i], color=GOLD, fontsize=9,
            ha="center", va="bottom")

    y_cell = grid_y0 + (n_ch - 1 - i) * (cell_size + gap) + cell_size / 2
    x_left = grid_x0 - 0.4
    ax.text(x_left, y_cell, ch_names[i], color=GOLD, fontsize=9,
            ha="right", va="center")

# Legend box
result_box(ax, 6.0, 8.5,
           "Relative NL coupling strength\n"
           "SPM: \u03b3\u00b7P\u00b7L_eff  (1.0\u00d7)\n"
           "XPM: 2\u03b3\u00b7P\u00b7L_eff  (2.0\u00d7)\n"
           "FWM: adds ~0.3\u00d7 at nearest neighbors",
           CYAN, fontsize=9)

# Axis labels
note(ax, 4.7, -2.5, "Target Channel  \u2192", SILVER, fontsize=11)
ax.text(-2.0, 4.5, "Source Channel  \u2192", color=SILVER, fontsize=11,
        ha="center", va="center", rotation=90)

save_fig(fig, "dwdm_03_crosstalk_matrix.png")


# ================================================================
# FIG 4: Q335 PRECISION FLOOR VS PHYSICAL SCALES
# Type: Scale/Landscape (log axis)
# Shows: 2^335 ~ 10^101 precision digits dwarfs every physical
#        scale in the problem — arithmetic is never the limit
# ================================================================

fig, ax = dark_fig(
    "Q335 Precision Floor vs Physical Scales",
    "Scale  (meters, seconds, or dimensionless ratio)",
    "",
    size=(18, 10)
)

# Log scale landmarks
landmarks = [
    (-101, "Q335 precision\nfloor (~10\u207b\u00b9\u2070\u00b9)", MAG, 14),
    (-35, "Planck length\n(1.6\u00d710\u207b\u00b3\u2075 m)", PURPLE, 10),
    (-24, "\u03b2\u2082 dispersion\n(10\u207b\u00b2\u2074 s\u00b2/m)", CYAN, 10),
    (-15, "Wavelength\n(1.55 \u03bcm = 10\u207b\u2076 m)", GREEN, 10),
    (-6, "Fiber core\n(80 \u03bcm\u00b2 = 10\u207b\u2075 m)", ORANGE, 10),
    (-3, "NL coefficient \u03b3\n(1.3\u00d710\u207b\u00b3 /W/m)", RED, 10),
    (0, "Channel power\n(1 mW = 10\u207b\u00b3 W)", GOLD, 10),
    (3, "Link length\n(2 km = 10\u00b3\u00b7\u2075 m)", BLUE, 10),
    (7, "Transoceanic\n(10\u2074 km = 10\u2077 m)", SILVER, 10),
    (14, "Frequency\n(193.1 THz = 10\u00b9\u2074 Hz)", WHITE, 10),
]

# Plot on a linear axis representing log10 scale
y_base = 0.5
positions = [lm[0] for lm in landmarks]
ax.set_xlim(-115, 25)
ax.set_ylim(-0.8, 2.5)

# Horizontal baseline
ax.plot([-110, 20], [y_base, y_base], color=DIM, lw=1.0, alpha=0.5)

# Tick marks and labels
for i in range(len(landmarks)):
    pos, label, color, fsize = landmarks[i]
    # Vertical tick
    ax.plot([pos, pos], [y_base - 0.08, y_base + 0.08], color=color, lw=2.0)
    data_point(ax, pos, y_base, None, color, size=80)

    # Alternate labels above and below to avoid overlap
    if i % 2 == 0:
        y_text = y_base + 0.45 + (i % 3) * 0.15
        va = "bottom"
    else:
        y_text = y_base - 0.45 - (i % 3) * 0.15
        va = "top"

    ax.text(pos, y_text, label, color=color, fontsize=fsize,
            ha="center", va=va, fontweight="bold" if i == 0 else "normal")

# Brace / region showing the gap
shaded_region(ax, -101, -35, MAG, alpha=0.04, label="66 orders of magnitude below Planck")

# Arrow showing the huge gap
ax.annotate("", xy=(-35, 1.9), xytext=(-101, 1.9),
            arrowprops=dict(arrowstyle="<->", color=MAG, lw=2.0))
ax.text(-68, 2.05, "66 orders of magnitude", color=MAG, fontsize=12,
        ha="center", va="bottom", fontweight="bold")

result_box(ax, -55, -0.45,
           "2\u00b3\u00b3\u2075 \u2248 10\u00b9\u2070\u00b9  \u2014  101 decimal digits of exact precision\n"
           "Arithmetic is never the limiting factor",
           GOLD, fontsize=10)

ax.set_yticks([])
legend(ax, loc="upper right")
save_fig(fig, "dwdm_04_q335_landscape.png")


# ================================================================
# FIG 5: EXACT VS APPROXIMATE PHASE ROTATION ERROR
# Type: Curve
# Shows: Error of sin(phi) ~ phi approximation vs phase magnitude.
#        Operating point (0.0166 rad) sits deep in valid regime.
# ================================================================

fig, ax = dark_fig(
    "First-Order Phase Rotation Error",
    "Phase magnitude  \u03c6  (radians)",
    "Relative error  |sin(\u03c6) \u2212 \u03c6| / |sin(\u03c6)|",
    size=(16, 10)
)

phi = np.linspace(1e-4, 1.0, 1000)
sin_phi = np.sin(phi)
rel_error = np.abs(sin_phi - phi) / np.abs(sin_phi)

curve(ax, phi, rel_error, "Relative error of sin(\u03c6) \u2248 \u03c6", CYAN, width=2.5)

# Operating point
phi_op = 0.016622
err_op = abs(np.sin(phi_op) - phi_op) / abs(np.sin(phi_op))
data_point(ax, phi_op, err_op, "Datacenter operating point", GOLD, size=250)
arrow_label(ax, phi_op, err_op, 0.20, 0.025,
            "Datacenter NL phase\n\u03c6 = 0.0166 rad\nerror = %.1e" % err_op, GOLD)

# Mark some reference thresholds
threshold_line(ax, None, "1% error", RED, vertical=False)
ax.axhline(y=0.01, color=RED, lw=1.2, linestyle="--", alpha=0.7)
note(ax, 0.82, 0.013, "1% error", RED, fontsize=9)

threshold_line(ax, None, "0.1% error", ORANGE, vertical=False)
ax.axhline(y=0.001, color=ORANGE, lw=1.2, linestyle="--", alpha=0.7)
note(ax, 0.82, 0.0013, "0.1% error", ORANGE, fontsize=9)

# Shaded valid region
shaded_region(ax, 0, 0.25, GREEN, alpha=0.05,
              label="High-accuracy regime (<0.1% error)")

ax.set_xlim(-0.05, 1.08)
ax.set_ylim(-0.01, 0.18)
ax.set_yscale("log")
ax.set_ylim(1e-8, 0.2)

result_box(ax, 0.55, 1e-7,
           "At \u03c6 = 0.017 rad, approximation\n"
           "error is ~5\u00d710\u207b\u2078 (negligible)\n"
           "No transcendental calls needed\n"
           "in the inner propagation loop",
           SILVER, fontsize=9)

legend(ax, loc="upper left")
save_fig(fig, "dwdm_05_phase_rotation_error.png")


# ================================================================
# FIG 6: CHANNEL GRID ON ITU FREQUENCY AXIS
# Type: Scale diagram
# Shows: Channel positions at 100/50/25/12.5 GHz spacings
#        on the ITU grid — how tight the frontier spacings are
# ================================================================

fig, ax = dark_canvas("ITU-T G.694.1 Channel Grid: Four Standard Spacings",
                      size=(18, 10))
ax.set_xlim(192.85, 193.45)
ax.set_ylim(-0.5, 10.5)

# Frequency axis at bottom
ax.plot([192.88, 193.42], [0.3, 0.3], color=DIM, lw=1.5)

# Tick marks every 25 GHz
for f in np.arange(192.9, 193.4, 0.025):
    ax.plot([f, f], [0.15, 0.45], color=DIM, lw=0.8, alpha=0.4)

note(ax, 193.15, -0.2, "Frequency (THz)", SILVER, fontsize=11)

# Anchor
ax.plot([193.1, 193.1], [0.0, 10.0], color=GOLD, lw=1.0, linestyle=":", alpha=0.4)
note(ax, 193.1, 10.2, "193.1 THz\n(ITU anchor)", GOLD, fontsize=9)

# Four spacing tiers, each row showing channels around anchor
tiers = [
    (8.5, "100 GHz", 0.1, BLUE, 8),
    (6.2, "50 GHz", 0.05, GREEN, 8),
    (3.9, "25 GHz", 0.025, ORANGE, 12),
    (1.6, "12.5 GHz", 0.0125, RED, 16),
]

for y_row, label, spacing, color, n_ch in tiers:
    # Center channels around anchor
    freqs = []
    for ci in range(n_ch):
        f = ANCHOR / 1e12 - (n_ch // 2 - 0.5) * spacing + ci * spacing
        if 192.88 < f < 193.42:
            freqs.append(f)

    # Channel markers
    for f in freqs:
        ax.plot([f, f], [y_row - 0.35, y_row + 0.35], color=color, lw=2.2)
        ax.plot(f, y_row + 0.35, "o", color=color, markersize=4)

    # Label
    ax.text(192.895, y_row, label, color=color, fontsize=12,
            fontweight="bold", ha="left", va="center")

    # Spacing bracket between first two channels
    if len(freqs) >= 2:
        f1, f2 = freqs[0], freqs[1]
        y_brace = y_row - 0.55
        ax.annotate("", xy=(f2, y_brace), xytext=(f1, y_brace),
                    arrowprops=dict(arrowstyle="<->", color=color, lw=1.0))
        ax.text((f1 + f2) / 2, y_brace - 0.25,
                "%.1f GHz" % (spacing * 1000),
                color=color, fontsize=8, ha="center", va="top")

# Capacity note
result_box(ax, 193.28, 9.2,
           "Tighter spacing = more channels\n"
           "12.5 GHz: 2\u00d7 capacity of 25 GHz\n"
           "But float drift worst at tight spacing",
           CYAN, fontsize=9)

save_fig(fig, "dwdm_06_itu_channel_grid.png")


# ================================================================
# FIG 7: DISPERSION PHASE VS FREQUENCY OFFSET
# Type: Curve (parabola)
# Shows: beta2 * omega^2 dispersion grows quadratically with
#        channel offset — outer channels see more dispersion
# ================================================================

fig, ax = dark_fig(
    "Chromatic Dispersion Phase vs Channel Offset",
    "Channel offset from 193.1 THz  (GHz)",
    "Dispersion phase per step  |\u03b2\u2082\u03c9\u00b2\u00b7dz|  (radians)",
    size=(16, 10)
)

offsets_ghz = np.linspace(-200, 200, 1000)
offsets_hz = offsets_ghz * 1e9
omega = 2.0 * np.pi * offsets_hz
disp_phase = np.abs(BETA2 * omega ** 2 * DZ_DC)

curve(ax, offsets_ghz, disp_phase, "|\u03b2\u2082|\u00b7(2\u03c0\u0394f)\u00b2\u00b7dz", CYAN, width=2.5)

# Mark datacenter channel offsets
dc_offsets = [-37.5, -12.5, 12.5, 37.5]
dc_colors = [RED, ORANGE, GREEN, BLUE]
dc_labels = ["Ch 0", "Ch 1", "Ch 2", "Ch 3"]

for i in range(4):
    off = dc_offsets[i]
    off_hz = off * 1e9
    phase = abs(BETA2 * (2.0 * np.pi * off_hz) ** 2 * DZ_DC)
    data_point(ax, off, phase, dc_labels[i], dc_colors[i], size=220)

# Annotations with staggered offsets
annot_positions = [
    (-37.5, -120, 3.8e-4),
    (-12.5, -90, 3.0e-4),
    (12.5, 90, 3.0e-4),
    (37.5, 120, 3.8e-4),
]
for i in range(4):
    off_data = dc_offsets[i]
    off_hz = off_data * 1e9
    phase = abs(BETA2 * (2.0 * np.pi * off_hz) ** 2 * DZ_DC)
    x_text, y_text = annot_positions[i][1], annot_positions[i][2]
    arrow_label(ax, off_data, phase, x_text, y_text,
                "%s\n\u0394f = %+.1f GHz\n\u03c6 = %.2e rad" % (dc_labels[i], off_data, phase),
                dc_colors[i])

# Symmetry annotation
ax.annotate("", xy=(37.5, 4.3e-4), xytext=(-37.5, 4.3e-4),
            arrowprops=dict(arrowstyle="<->", color=GOLD, lw=1.5))
ax.text(0, 4.5e-4, "Symmetric: Ch0 = Ch3, Ch1 = Ch2", color=GOLD,
        fontsize=10, ha="center", va="bottom")

ax.set_xlim(-220, 220)
ax.set_ylim(-2e-5, 5.5e-4)

result_box(ax, 100, 1.5e-4,
           "\u03b2\u2082 = \u221220.41 ps\u00b2/km\n"
           "dz = 400 m\n"
           "Parabolic \u2192 quadratic growth\n"
           "with channel offset",
           SILVER, fontsize=9)

legend(ax, loc="upper center")
save_fig(fig, "dwdm_07_dispersion_phase.png")


# ================================================================
# FIG 8: DWDM DEPLOYMENT LANDSCAPE
# Type: Scale/Landscape diagram
# Shows: Datacenter / Metro / Transoceanic on log scale with
#        step counts, channel counts, and key parameters
# ================================================================

fig, ax = dark_canvas("DWDM Deployment Landscape: Three Regimes",
                      size=(18, 10))
ax.set_xlim(-1, 15)
ax.set_ylim(-1.5, 7.5)

# Log-distance axis: 2 km -> 10,000 km
# Map log10(km) to x: log10(2)=0.3, log10(160)=2.2, log10(10000)=4
# Scale: x = 2 + (log10(km) - 0.3) * 3.0
def km_to_x(km):
    return 2.0 + (np.log10(km) - 0.3) * 3.0

# Axis line
ax.plot([0.5, 14.5], [1.0, 1.0], color=DIM, lw=1.5)
ax.annotate("", xy=(14.5, 1.0), xytext=(14.0, 1.0),
            arrowprops=dict(arrowstyle="->", color=DIM, lw=1.5))

# Distance tick marks
dist_ticks = [1, 2, 5, 10, 50, 100, 500, 1000, 5000, 10000]
for d in dist_ticks:
    x = km_to_x(d)
    if 0.5 < x < 14.5:
        ax.plot([x, x], [0.85, 1.15], color=DIM, lw=0.8)
        ax.text(x, 0.55, "%s km" % d, color=DIM, fontsize=7,
                ha="center", va="top")

note(ax, 7.5, 0.1, "Propagation Distance (log scale)", SILVER, fontsize=10)

# Three deployment boxes
deployments = [
    {
        "name": "Datacenter",
        "km": 2,
        "color": GREEN,
        "channels": "4\u201316",
        "spacing": "25 GHz",
        "spans": "1",
        "steps": "5\u201320",
        "status": "WORKING",
        "drift": "9.62%",
    },
    {
        "name": "Metro",
        "km": 160,
        "color": ORANGE,
        "channels": "4\u201340",
        "spacing": "50\u2013100 GHz",
        "spans": "2\u20136",
        "steps": "10\u2013600",
        "status": "PROPAGATION OK\nAMPLIFICATION DEFERRED",
        "drift": "\u221e (model issue)",
    },
    {
        "name": "Transoceanic",
        "km": 10000,
        "color": CYAN,
        "channels": "96\u2013192",
        "spacing": "50 GHz",
        "spans": "125",
        "steps": "12,500",
        "status": "DEFERRED",
        "drift": "Projected: catastrophic\nfor float64",
    },
]

box_y_positions = [3.5, 3.5, 3.5]
box_height = 3.2
box_width = 3.8

for i in range(len(deployments)):
    dep = deployments[i]
    x_center = km_to_x(dep["km"])
    x_box = x_center - box_width / 2

    # Clamp boxes to not overlap
    if i == 0:
        x_box = 0.8
    elif i == 1:
        x_box = 5.2
    elif i == 2:
        x_box = 9.8

    y_box = box_y_positions[i]

    # Box
    rect = FancyBboxPatch(
        (x_box, y_box), box_width, box_height,
        boxstyle="round,pad=0.25",
        facecolor=dep["color"], alpha=0.08,
        edgecolor=dep["color"], linewidth=1.5
    )
    ax.add_patch(rect)

    # Title
    ax.text(x_box + box_width / 2, y_box + box_height - 0.3,
            dep["name"], color=dep["color"], fontsize=13,
            fontweight="bold", ha="center", va="top")

    # Info lines with generous internal margin
    info_lines = [
        "Distance: %s km" % dep["km"],
        "Channels: %s" % dep["channels"],
        "Spacing: %s" % dep["spacing"],
        "Spans: %s" % dep["spans"],
        "Steps: %s" % dep["steps"],
        "Drift: %s" % dep["drift"],
    ]
    for li in range(len(info_lines)):
        ax.text(x_box + 0.35, y_box + box_height - 0.85 - li * 0.38,
                info_lines[li], color=SILVER, fontsize=8,
                ha="left", va="top")

    # Connector line down to axis
    x_tick = km_to_x(dep["km"])
    ax.plot([x_tick, x_tick], [1.15, y_box], color=dep["color"],
            lw=1.2, linestyle=":", alpha=0.5)

# Status badges
status_info = [
    (2.7, 3.15, "VERIFIED EXACT", GREEN),
    (7.1, 3.15, "ENGINE OK / MODEL DEFERRED", ORANGE),
    (11.7, 3.15, "TARGET", CYAN),
]
for sx, sy, slabel, scolor in status_info:
    ax.text(sx, sy, slabel, color=scolor, fontsize=8,
            fontweight="bold", ha="center", va="top",
            bbox=dict(boxstyle="round,pad=0.2", facecolor=BG,
                      edgecolor=scolor, alpha=0.8))

save_fig(fig, "dwdm_08_deployment_landscape.png")


# ================================================================
# SUMMARY
# ================================================================
filenames = [
    "dwdm_01_fwm_efficiency.png",
    "dwdm_02_fwm_triplet.png",
    "dwdm_03_crosstalk_matrix.png",
    "dwdm_04_q335_landscape.png",
    "dwdm_05_phase_rotation_error.png",
    "dwdm_06_itu_channel_grid.png",
    "dwdm_07_dispersion_phase.png",
    "dwdm_08_deployment_landscape.png",
]

print("\n" + "=" * 72)
print("DWDM-VDR DIAGRAMS COMPLETE")
print("=" * 72)
for i in range(len(filenames)):
    print("  Fig %d: %s" % (i + 1, filenames[i]))
print("=" * 72)
