"""
generate_experimental_data.py
==============================
Generates synthetic experimental data and PNG figures for the HLV–AILEE
coherence-gated propulsion model.

Simulates what a real test-stand campaign might produce:
  - Thrust vs. power sweeps at multiple misalignment levels
  - Hysteresis loops from up/down timing-offset sweeps
  - Time-series coupling efficiency from a drifting burn
  - Spectral mismatch proxy over a frequency sweep
  - Gate collapse threshold scan

All data is saved as .csv alongside publication-quality .png figures.

MIT License — Krüger & Feeney (2025)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.signal import welch
from pathlib import Path

# ── output directories ───────────────────────────────────────────────────────
OUT = Path("/mnt/user-data/outputs")
DATA_DIR = OUT / "experimental_data"
FIG_DIR  = OUT / "figures"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

RNG   = np.random.default_rng(42)
PHI   = (1 + np.sqrt(5)) / 2
G0    = 9.80665
ALPHA = 2.5          # gate sharpness (fitted from synthetic threshold scan)
ISP   = 350.0        # [s]
ETA   = 0.85

# ── plot style ───────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.dpi": 150,
    "font.family": "DejaVu Sans",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "legend.fontsize": 9,
    "lines.linewidth": 1.8,
})
COLORS = {
    "classical": "#888888",
    "low":       "#2e86c1",
    "mid":       "#e67e22",
    "high":      "#c0392b",
    "up":        "#2e86c1",
    "down":      "#e74c3c",
    "gate":      "#1a5276",
    "composite": "#7d3c98",
    "fit":       "#117a65",
}

def add_noise(arr, snr=40):
    """Add Gaussian noise at a given SNR."""
    signal_power = np.mean(arr**2)
    noise_std = np.sqrt(signal_power / snr)
    return arr + RNG.normal(0, noise_std, arr.shape)

def gate(dphi, alpha=ALPHA):
    return np.exp(-alpha * dphi**2)

def save(fig, name):
    path = FIG_DIR / f"{name}.png"
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  saved → {path.name}")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 1 — Thrust vs. Power: gate-limited saturation
# ═══════════════════════════════════════════════════════════════════════════
print("Generating Figure 1: Thrust vs. Power saturation ...")

N = 400
t        = np.linspace(0, 80, N)
mass     = np.linspace(950, 580, N)
velocity = np.linspace(120, 380, N)

power_kw  = np.linspace(10, 1800, 40)
phi_cases = {"Δφ=0.00 (classical)": 0.00,
             "Δφ=0.30 (low)":       0.30,
             "Δφ=0.80 (mid)":       0.80,
             "Δφ=1.40 (high)":      1.40}
color_map  = [COLORS["classical"], COLORS["low"], COLORS["mid"], COLORS["high"]]

rows = {"power_kw": power_kw}
fig, ax = plt.subplots(figsize=(8, 4.5))

for (label, dphi), color in zip(phi_cases.items(), color_map):
    dvs = []
    for P in power_kw * 1e3:
        p_in     = np.full(N, P)
        G        = gate(dphi)
        integrand = (p_in / (mass * velocity)) * G
        dvs.append(float(ISP * ETA * np.trapezoid(integrand, t)))
    dvs = add_noise(np.array(dvs), snr=35)
    rows[label] = dvs
    ax.plot(power_kw, dvs, label=label, color=color)

np.savetxt(DATA_DIR / "thrust_vs_power.csv",
           np.column_stack(list(rows.values())),
           delimiter=",",
           header=",".join(rows.keys()),
           comments="")

ax.set_xlabel("Input power P_in [kW]")
ax.set_ylabel("Velocity gain Δv [m s⁻¹]")
ax.set_title("Figure 1 — Gate-limited saturation: Δv vs. input power")
ax.legend()
save(fig, "fig1_thrust_vs_power")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2 — Hysteresis loop (timing-offset sweep)
# ═══════════════════════════════════════════════════════════════════════════
print("Generating Figure 2: Hysteresis loop ...")

phi_sweep = np.linspace(0, 2.2, 120)
# up-sweep collapses faster (steeper effective alpha)
G_up   = add_noise(gate(phi_sweep, alpha=ALPHA * 1.35), snr=30)
G_down = add_noise(gate(phi_sweep, alpha=ALPHA * 0.75), snr=30)
G_up   = np.clip(G_up,   0, 1)
G_down = np.clip(G_down, 0, 1)

np.savetxt(DATA_DIR / "hysteresis_loop.csv",
           np.column_stack([phi_sweep, G_up, G_down]),
           delimiter=",",
           header="delta_phi,G_up_sweep,G_down_sweep",
           comments="")

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(phi_sweep, G_up,   color=COLORS["up"],   label="Up-sweep (collapse branch)")
ax.plot(phi_sweep, G_down, color=COLORS["down"],  label="Down-sweep (recovery branch)", ls="--")
ax.fill_between(phi_sweep, G_up, G_down,
                where=(G_down > G_up), alpha=0.12, color="purple",
                label="Hysteresis window")
ax.axvline(phi_sweep[np.argmin(np.abs(G_up - 0.5))],
           color="grey", lw=0.9, ls=":", label="50 % collapse (up)")
ax.set_xlabel("Phase deviation Δφ  (timing-offset proxy)")
ax.set_ylabel("Effective coupling G_eff")
ax.set_title("Figure 2 — Phase-window hysteresis: up/down timing sweep")
ax.legend()
save(fig, "fig2_hysteresis_loop")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 3 — Time-series: gate G(t) and cumulative Δv during a drifting burn
# ═══════════════════════════════════════════════════════════════════════════
print("Generating Figure 3: Drifting burn time-series ...")

N2   = 600
t2   = np.linspace(0, 120, N2)
mass2 = np.linspace(1000, 550, N2)
vel2  = np.linspace(100,  420, N2)
p2    = np.full(N2, 6e5)

phi_drift = 0.04 + 0.9 * (t2 / t2[-1])**2
G_drift   = add_noise(gate(phi_drift), snr=50)
G_drift   = np.clip(G_drift, 0, 1)
G_ideal   = np.ones(N2)

intg_gated   = (p2 / (mass2 * vel2)) * G_drift
intg_ungated = (p2 / (mass2 * vel2)) * G_ideal
cum_gated    = ISP * ETA * np.array(
    [np.trapezoid(intg_gated[:i+1], t2[:i+1]) for i in range(N2)])
cum_ungated  = ISP * ETA * np.array(
    [np.trapezoid(intg_ungated[:i+1], t2[:i+1]) for i in range(N2)])

np.savetxt(DATA_DIR / "drifting_burn_timeseries.csv",
           np.column_stack([t2, phi_drift, G_drift, cum_gated, cum_ungated]),
           delimiter=",",
           header="time_s,delta_phi,gate_G,cumulative_dv_gated,cumulative_dv_ungated",
           comments="")

fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
axes[0].plot(t2, phi_drift, color=COLORS["mid"],  label="Δφ(t) — measured proxy")
axes[0].set_ylabel("Phase deviation Δφ")
axes[0].set_title("Figure 3 — Drifting burn: phase deviation, gate, and cumulative Δv")
axes[0].legend()

axes[1].plot(t2, G_drift, color=COLORS["gate"], label="G(t) measured")
axes[1].axhline(1.0, color="grey", lw=0.8, ls="--", label="Classical (G=1)")
axes[1].set_ylabel("Gate G(t)")
axes[1].legend()

axes[2].plot(t2, cum_gated,   color=COLORS["high"],      label=f"HLV–AILEE  Δv={cum_gated[-1]:.1f} m/s")
axes[2].plot(t2, cum_ungated, color=COLORS["classical"],  label=f"Classical  Δv={cum_ungated[-1]:.1f} m/s", ls="--")
axes[2].set_ylabel("Cumulative Δv [m s⁻¹]")
axes[2].set_xlabel("Time [s]")
axes[2].legend()

plt.tight_layout()
save(fig, "fig3_drifting_burn_timeseries")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 4 — Spectral proxy: PSD mismatch across a frequency sweep
# ═══════════════════════════════════════════════════════════════════════════
print("Generating Figure 4: Spectral proxy sweep ...")

SR     = 2048.0
t_sig  = np.linspace(0, 4, int(4 * SR), endpoint=False)
f_lat  = 100.0 * PHI   # lattice reference at φ·100 Hz ≈ 161.8 Hz
f_sys_sweep = np.linspace(80, 260, 35)
phi_spec_vals = []

for f_sys in f_sys_sweep:
    thrust = np.sin(2 * np.pi * f_sys  * t_sig) + RNG.normal(0, 0.04, len(t_sig))
    lat    = np.sin(2 * np.pi * f_lat  * t_sig)
    _, S_t = welch(thrust, fs=SR, nperseg=512)
    _, S_l = welch(lat,    fs=SR, nperseg=512)
    dfreq  = SR / 512
    S_t /= np.trapezoid(S_t) * dfreq or 1
    S_l /= np.trapezoid(S_l) * dfreq or 1
    phi_spec_vals.append(np.trapezoid(np.abs(S_t - S_l)) * dfreq)

phi_spec_vals = np.array(phi_spec_vals)

np.savetxt(DATA_DIR / "spectral_proxy_sweep.csv",
           np.column_stack([f_sys_sweep, phi_spec_vals]),
           delimiter=",",
           header="f_sys_hz,delta_phi_spec",
           comments="")

fig, axes = plt.subplots(2, 1, figsize=(8, 7))
axes[0].plot(f_sys_sweep, phi_spec_vals, color=COLORS["mid"], marker="o", ms=4)
axes[0].axvline(f_lat, color=COLORS["gate"], lw=1.2, ls="--",
                label=f"f_lat = {f_lat:.1f} Hz  (φ·100)")
axes[0].set_ylabel("Δφ_spec  (spectral L¹ distance)")
axes[0].set_title("Figure 4 — Spectral proxy Δφ_spec across system-frequency sweep")
axes[0].legend()

G_spec = gate(phi_spec_vals)
axes[1].plot(f_sys_sweep, G_spec, color=COLORS["gate"], marker="o", ms=4)
axes[1].axvline(f_lat, color=COLORS["gate"], lw=1.2, ls="--",
                label="Peak coupling at f_lat")
axes[1].set_xlabel("System frequency f_sys [Hz]")
axes[1].set_ylabel("Gate G  (predicted coupling)")
axes[1].set_title("Predicted coupling efficiency from spectral proxy")
axes[1].legend()

plt.tight_layout()
save(fig, "fig4_spectral_proxy_sweep")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 5 — Gate collapse threshold scan + alpha fit
# ═══════════════════════════════════════════════════════════════════════════
print("Generating Figure 5: Gate collapse threshold scan ...")

from scipy.optimize import curve_fit

phi_scan = np.linspace(0, 2.0, 50)
G_true   = gate(phi_scan, alpha=ALPHA)
G_obs    = add_noise(G_true, snr=25)
G_obs    = np.clip(G_obs, 0, 1)

popt, _ = curve_fit(lambda phi, a: np.exp(-a * phi**2),
                    phi_scan, G_obs, p0=[2.0], bounds=(0, np.inf))
alpha_fit = float(popt[0])
G_fit     = gate(phi_scan, alpha=alpha_fit)

np.savetxt(DATA_DIR / "gate_threshold_scan.csv",
           np.column_stack([phi_scan, G_obs, G_fit]),
           delimiter=",",
           header="delta_phi,G_observed,G_fitted",
           comments="")

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.scatter(phi_scan, G_obs, color=COLORS["mid"], s=22, zorder=3,
           label="Synthetic observations")
ax.plot(phi_scan, G_fit, color=COLORS["fit"],
        label=f"Fitted gate  α={alpha_fit:.3f}  (true α={ALPHA})")
ax.plot(phi_scan, G_true, color="grey", ls="--", lw=1,
        label=f"Ground truth  α={ALPHA}")
ax.axhline(np.exp(-1), color="grey", lw=0.8, ls=":",
           label="e⁻¹ ≈ 0.368  (half-width reference)")
ax.set_xlabel("Phase deviation Δφ")
ax.set_ylabel("Coupling efficiency G")
ax.set_title("Figure 5 — Gate collapse threshold scan with α calibration fit")
ax.legend()
save(fig, "fig5_gate_threshold_scan")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 6 — Summary panel: all four Δφ proxies side-by-side
# ═══════════════════════════════════════════════════════════════════════════
print("Generating Figure 6: Summary proxy panel ...")

fig = plt.figure(figsize=(12, 8))
gs  = gridspec.GridSpec(2, 2, hspace=0.38, wspace=0.32)

# ── temporal ──
ax0 = fig.add_subplot(gs[0, 0])
f_sys_drift = 100 + np.linspace(0, 35, N2) + RNG.normal(0, 1.2, N2)
phi_temp    = np.abs(f_sys_drift - 100.0) / 100.0
ax0.plot(t2, phi_temp, color=COLORS["low"], lw=1.2)
ax0.set_title("Temporal proxy  Δφ_temp  (Eq. 7)")
ax0.set_xlabel("Time [s]")
ax0.set_ylabel("Δφ_temp")

# ── spatial (vorticity ratio — simulated) ──
ax1 = fig.add_subplot(gs[0, 1])
phi_spat = 0.05 + 0.6 * np.abs(np.sin(np.linspace(0, 3*np.pi, N2)))
phi_spat = add_noise(phi_spat, snr=20)
ax1.plot(t2, phi_spat, color=COLORS["mid"], lw=1.2)
ax1.set_title("Spatial proxy  Δφ_spat  (Eq. 8)")
ax1.set_xlabel("Time [s]")
ax1.set_ylabel("Δφ_spat")

# ── spectral ──
ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(f_sys_sweep, phi_spec_vals, color=COLORS["high"], marker="o", ms=3)
ax2.axvline(f_lat, color="grey", lw=1, ls="--")
ax2.set_title("Spectral proxy  Δφ_spec  (Eq. 9)")
ax2.set_xlabel("f_sys [Hz]")
ax2.set_ylabel("Δφ_spec")

# ── composite ──
ax3 = fig.add_subplot(gs[1, 1])
w = np.array([0.5, 0.5])
phi_comp_ts = np.sqrt(w[0] * phi_temp**2 + w[1] * phi_spat**2)
G_comp      = gate(phi_comp_ts)
ax3.plot(t2, phi_comp_ts, color=COLORS["composite"], lw=1.2, label="Δφ_composite")
ax3.plot(t2, G_comp,      color=COLORS["gate"],      lw=1.2, ls="--", label="G(t)")
ax3.set_title("Composite Δφ and resulting gate G(t)")
ax3.set_xlabel("Time [s]")
ax3.set_ylabel("Value")
ax3.legend()

fig.suptitle("Figure 6 — HLV–AILEE: synthetic proxy measurements summary",
             fontsize=13, y=1.01)
save(fig, "fig6_proxy_summary_panel")


# ═══════════════════════════════════════════════════════════════════════════
# Done
# ═══════════════════════════════════════════════════════════════════════════
print("\n── Complete ─────────────────────────────────────────")
print(f"  CSV files → {DATA_DIR}")
print(f"  PNG files → {FIG_DIR}")
csv_files = sorted(DATA_DIR.glob("*.csv"))
png_files = sorted(FIG_DIR.glob("*.png"))
print(f"\n  {len(csv_files)} datasets:")
for f in csv_files:
    rows = sum(1 for _ in open(f)) - 1
    print(f"    {f.name}  ({rows} rows)")
print(f"\n  {len(png_files)} figures:")
for f in png_files:
    print(f"    {f.name}")
