# HLV–AILEE · Geometry-Limited Propulsion

> *Coherence-Gated Momentum Transfer — simulation code for*  
> **"Geometry-Limited Propulsion: An HLV–AILEE Formulation of Coherence-Gated Momentum Transfer"**  
> Krüger & Feeney (2025)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![NumPy](https://img.shields.io/badge/NumPy-required-013243?logo=numpy)](https://numpy.org/)
[![SciPy](https://img.shields.io/badge/SciPy-optional-8CAAE6?logo=scipy)](https://scipy.org/)
[![Jupyter](https://img.shields.io/badge/Jupyter-notebook-F37626?logo=jupyter)](https://jupyter.org/)
[![CI](https://github.com/dfeen87/Geometry-Limited-Propulsion/workflows/CI/badge.svg)](https://github.com/dfeen87/Geometry-Limited-Propulsion/actions)

---

## Overview

Classical propulsion models absorb saturation losses into empirical efficiency terms. This repository implements the **HLV–AILEE** (Helix–Light–Vortex / Admissible-Interaction Lattice Energy–Exchange) formulation, which proposes a deeper explanation: momentum transfer is gated by the *phase and geometric alignment* between propulsion dynamics and an admissible host-lattice structure.

The central claim is that a coherence gate

$$G(t) = e^{-\alpha\,\Delta\phi(t)^2}$$

controls realizable energy–momentum coupling. When alignment is good (Δφ ≈ 0), coupling proceeds normally. When alignment degrades, the gate collapses — and additional input power *cannot compensate*.

The model makes **no reactionless or superluminal claims**. All predictions are falsifiable. See §6 of the paper for the full falsification protocol.

> A note on current status: The simulation code implements the HLV–AILEE mathematical framework faithfully and produces clean, falsifiable predictions. However, the bridge between the measurable Δφ proxies and the underlying lattice geometry remains theoretical — no empirical validation has been performed yet. What this repository provides is a numerically precise, testable formulation: the experiments needed to confirm or falsify it are the logical next step, and the measurement protocol in simulation_stub.ipynb and §7 of the paper is designed with exactly that in mind.

---

## Repository Structure

```
hlv-ailee/
├── README.md
├── LICENSE
├── requirements.txt            # Runtime dependencies
├── .gitignore                  # Common Python/Jupyter ignore rules
├── CITATION.cff                # Citation metadata for GitHub/Zenodo
├── CONTRIBUTING.md             # Contribution guidelines
├── CODE_OF_CONDUCT.md          # Community standards
├── src/
│   ├── coherence_gate.py           # Core gate function and velocity-gain law
│   ├── phase_alignment_metrics.py  # Measurable Δφ proxies (Eqs. 7–9)
│   └── simulation_stub.ipynb       # End-to-end simulation notebook
└── paper/
    └── propulsion_hlv_ailee.md     # Paper in Markdown (Unicode math)
```

---

## Quickstart

**1. Clone and install dependencies**

```bash
git clone https://github.com/dfeen87/Geometry-Limited-Propulsion.git
cd Geometry-Limited-Propulsion
pip install -r requirements.txt
# Optional (needed only for delta_phi_spectral):
pip install scipy
```

SciPy is optional — it is only required for the spectral proxy (`delta_phi_spectral`). All other functionality runs on NumPy alone.

**2. Run the notebook**

```bash
cd src
jupyter notebook simulation_stub.ipynb
```

**3. Use the library directly**

```python
import numpy as np
from coherence_gate import coherence_gate, hlv_ailee_delta_v

# Gate value at a given phase deviation
G = coherence_gate(delta_phi=0.5, alpha=2.0)
print(f"G = {G:.4f}")   # → 0.7788

# Gated velocity gain over a burn
N = 500
t        = np.linspace(0, 100, N)
p_in     = np.full(N, 5e5)           # 500 kW constant
mass     = np.linspace(1000, 600, N) # propellant burn [kg]
velocity = np.linspace(100, 400, N)  # reference trajectory [m/s]
delta_phi = 0.05 + 0.8 * (t / t[-1])**2  # slowly drifting misalignment

dv = hlv_ailee_delta_v(
    t, p_in, mass, velocity,
    delta_phi=delta_phi,
    isp=350.0, eta=0.85, alpha=2.0,
)
print(f"Δv = {dv:.2f} m/s")
```

---

## Module Reference

### `coherence_gate.py`

| Function | Description |
|---|---|
| `coherence_gate(delta_phi, alpha)` | Gate function G(t) = exp(−α·Δφ²); returns values in (0, 1] |
| `classical_delta_v(isp, m0, mf)` | Tsiolkovsky baseline Δv = Isp·g₀·ln(M₀/Mf) |
| `hlv_ailee_delta_v(t, p_in, mass, velocity, delta_phi, isp, eta, alpha)` | Full numerical integration of the gated velocity-gain law (Eq. 6) |
| `gate_vs_phase_sweep(delta_phi_values, alpha)` | Sweep G over a range of Δφ — useful for plotting acceptance-window shape |
| `hysteresis_sweep(phi_up, phi_down, alpha_up, alpha_down)` | Asymmetric up/down sweep for hysteresis loop visualization (Figure 2) |

### `phase_alignment_metrics.py`

| Function | Proxy | Paper equation |
|---|---|---|
| `delta_phi_temporal(f_sys, f_lat)` | Fractional frequency mismatch | Eq. 7 |
| `delta_phi_spatial(vflow, dx, dy)` | Vorticity-to-speed ratio via 2-D curl | Eq. 8 |
| `delta_phi_spectral(thrust, lat, sr)` | L¹ distance between normalised PSDs | Eq. 9 |
| `composite_delta_phi(*proxies, weights)` | Weighted RMS combination of any proxies | — |

---

## Key Predictions

The table below summarises how HLV–AILEE predictions differ from classical expectations, and how each can be tested.

| Observable | Classical expectation | HLV–AILEE prediction | Test |
|---|---|---|---|
| Thrust vs. power | Smooth saturation | **Threshold-like collapse** beyond Δφ window | Sweep power at fixed geometry; look for non-linearity |
| Timing offsets | Small, continuous degradation | **Sharp non-linear drop** sensitive to phase drift | Introduce controlled timing offsets |
| Ignition / shutdown history | Mostly irrelevant | **Hysteresis / path dependence** | Up/down sweeps of timing and geometry controls |
| Reduced-power regime | Lower Δv | Potential **"sweet spot"** under improved alignment | Reduce power while improving phase alignment |

The model is **falsified** if:

1. No measurable correlation exists between Δφ proxies and coupling efficiency,
2. Increased power fully compensates induced misalignment (no gate-collapse regime), or
3. Controlled phase drift produces no threshold-like degradation or hysteresis.

---

## Notebook Walkthrough

`simulation_stub.ipynb` covers six sections:

1. **Gate shape** — G(Δφ) plotted for α ∈ {0.5, 1.0, 2.0, 5.0}; shows narrowing acceptance window with increasing sharpness.
2. **Hysteresis loop** — Asymmetric α values on up/down sweep branches reproduce Figure 2 of the paper, with shaded hysteresis window.
3. **Δv vs. power** — Three-way comparison (Δφ = 0, 0.3, 1.2) demonstrates gate-limited saturation vs. classical linear scaling.
4. **Temporal proxy** — Simulated frequency drift mapped through Eq. 7; plots f_sys(t) alongside Δφ_temp(t).
5. **Spectral proxy** — Welch PSD comparison of a 50 Hz thrust tone against aligned and misaligned lattice references; Δφ_spec printed for each case.
6. **Full integration** — Three-panel plot: (i) time-varying Δφ and composite, (ii) gate G(t) over the burn, (iii) cumulative Δv gated vs. ungated with final suppression factor.

---

## Measurement Protocol

For physical experiments, Section 7 of the paper outlines a four-phase protocol:

```
Phase 1 — Baseline
    Measure thrust, Isp, input power, spectra.
    Log geometry and environment.

Phase 2 — Perturbations
    Introduce timing offsets.
    Modulate geometry.
    Induce controlled phase drift.
    Record ignition/shutdown sequences.

Phase 3 — Analysis
    Fit α from observed threshold data.
    Plot hysteresis loops (up vs. down sweeps).
    Compare against null/surrogate baselines.

Phase 4 — Replication
    Pre-register all metrics and thresholds.
    Release raw data and configuration.
    Independent lab replication.
```

---

## Citation

If you use this code in your work, please cite the paper:

```bibtex
@misc{kruger_feeney_2025_hlv_ailee,
  author       = {Krüger, Marcel and Feeney, Don},
  title        = {Geometry-Limited Propulsion: An {HLV–AILEE} Formulation
                  of Coherence-Gated Momentum Transfer},
  year         = {2025},
  note         = {Manuscript draft},
}
```

The HLV mathematical framework this builds on:

```bibtex
@misc{kruger_2025_hlv,
  author       = {Krüger, Marcel},
  title        = {A Mathematical Unification of the Helix–Light–Vortex ({HLV})
                  Framework: Discrete Geometry, Spiral Time, Unified
                  Lagrangians, and {TOE}-Level Structure},
  year         = {2025},
  publisher    = {Zenodo},
}
```

---

## Authors

**Marcel Krüger** · Independent Researcher, HLV–RFP (Meiningen, Germany)  
✉ marcelkrueger092@gmail.com · ORCID [0009-0002-5709-9729](https://orcid.org/0009-0002-5709-9729)

**Don Feeney** · Independent Researcher (Pennsylvania, USA)  
✉ dfeen87@gmail.com · ORCID [0009-0003-1350-4160](https://orcid.org/0009-0003-1350-4160)

---

## License

This project is released under the [MIT License](LICENSE).

---

## Acknowledgement

Thank you to Claude (Anthropic) for assisting with the simulation code and PAPER.md. Thank you to Codex (OpenAI) for the final code review. All scientific content, theoretical framework, and experimental claims originate solely with the authors.
