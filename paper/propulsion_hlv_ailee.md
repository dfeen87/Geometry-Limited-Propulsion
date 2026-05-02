# Geometry-Limited Propulsion:
## An HLV–AILEE Formulation of Coherence-Gated Momentum Transfer

**Marcel Krüger**
Independent Researcher, HLV–RFP (Meiningen, Germany)
✉ marcelkrueger092@gmail.com · ORCID: 0009-0002-5709-9729

**Don Feeney**
Independent Researcher (Pennsylvania, USA)
✉ dfeen87@gmail.com · ORCID: 0009-0003-1350-4160

*May 2, 2026*

---

## Abstract

Classical propulsion models encode efficiency losses through empirical terms,
implicitly assuming near-optimal coupling between injected power and momentum
transfer. Yet many systems display sharp saturation behavior where additional
energy fails to produce proportional velocity gain. We propose that such
saturation reflects a deeper constraint: phase and geometric alignment between
propulsion dynamics and an admissible interaction structure of an emergent
substrate. Using the Helix–Light–Vortex framework — in which spacetime is
modeled as a discrete golden-ratio quasicrystal lattice governed by triadic
spiral time — we formalize a coherence gate that exponentially suppresses
energy–momentum coupling under phase deviation. The resulting HLV–AILEE
velocity-gain law preserves conservation principles and makes no reactionless
or superluminal claims. Instead, it yields falsifiable predictions:

1. threshold-like collapse of effective performance beyond a phase-alignment window,
2. hysteresis/path dependence under timing and geometric drift, and
3. counterintuitive optimal regimes at reduced power under improved coherence.

We outline measurement protocols and computational tests accessible with
existing instrumentation.

---

## 1  Introduction: Saturation in Classical Propulsion

Classical propulsion formulations typically assume that increasing input power
and improving engineering efficiency lead to proportional gains in thrust or
velocity, modulo known loss channels (thermal, turbulent, nozzle, and material
constraints). In practice, many systems exhibit sharp saturation: beyond an
operating regime, additional injected energy fails to yield proportional
momentum gain and instead dissipates into heat, noise, or unstructured flow.

This paper proposes that such saturation can be reinterpreted as a **constraint
problem**: momentum transfer is realizable only when the propulsion dynamics
remain sufficiently aligned with an admissible phase/geometry window. We
formalize this within the Helix–Light–Vortex framework and derive a
gate-controlled coupling law (HLV–AILEE) that is empirically distinguishable
from generic "efficiency tuning."

---

## 2  The Missing Variable: Geometry and Phase Alignment

### 2.1  From Empirical Efficiency to an Explicit Coupling Constraint

Consider the standard rocket equation:

> **Δv = Isp · g₀ · ln(M₀ / Mf)**  ·················· (1)

Here *I*sp is an engineering parameter. In the HLV–AILEE view, realizable
coupling also depends on whether the system remains within an admissible phase
window relative to the host geometry. This motivates introducing a first-class
alignment variable rather than absorbing everything into *η*.

---

## 3  HLV Host Geometry (Minimal Ingredients)

### 3.1  Discrete Quasicrystal Host (𝒢_φ)

Within HLV, physical spacetime is modeled as a 3D projection of a
higher-dimensional quasicrystalline lattice (cut-and-project). We denote the
host schematically by 𝒢_φ with golden ratio φ = (1 + √5) / 2 and acceptance
window W_φ:

> **𝒢_φ = { π‖(x) ∈ E‖ : x ∈ ℤ⁶, π⊥(x) ∈ W_φ }**  ·············· (2)

### 3.2  Triadic Spiral Time

HLV extends time to a triadic operator:

> **ψ(t) = t + i·φ(t) + j·χ(t)**  ·················· (3)

where *t* is the U₁ component, φ(t) a coherence mode (U₂), and χ(t) a
memory/retrocausal mode (U₃). For propulsion modeling we only require that
this structure can induce small controlled deformations of propagation:

> **A(t) = 1 + ε(t),   |ε(t)| ≪ 1**  ················ (4)

without violating conservation laws.

---

## 4  Phase Deviation as a Control Variable Δφ

### 4.1  Definition (Operational)

We define a phase/geometric deviation metric that compares system dynamics to
the host admissible structure:

> **Δφ(t) = |∂_ψ ψ_sys(t) − ψ_lat(t)|**  ················ (5)

### 4.2  Coherence-Gate Behavior

Δφ does not introduce a new force. It gates coupling efficiency:

- **Δφ ≃ 0** → coherent coupling into momentum transfer
- **Δφ outside acceptance window** → rapid collapse of realizable coupling

---

## 5  HLV–AILEE Law: Geometry-Gated Velocity Gain

### 5.1  Core Equations

The HLV–AILEE gate is defined as:

> **G(t) = exp(−α · Δφ(t)²)**  ·················· (6)

where α > 0 controls the width of the admissible phase window.

For **mass-flow propulsion**, the primary rocket-equation form is:

> **dv = g₀ · Isp · G(t) · (−dM / M)**  ················ (7)

For G(t) = 1 this reduces exactly to the classical Tsiolkovsky equation.

For **power-limited electric or ion propulsion**, a separate power-based
variant may be written as:

> **Δv = η · ∫₀^tf [ P_in(t) / (M(t) · v(t)) ] · G(t) dt**  ·········· (8)

The two forms describe different propulsion regimes and should not be
interpreted as the same physical law.

### 5.2  Interpretation

The gate factor G(t) does not create additional momentum and does not modify
conservation laws. It only reduces the fraction of otherwise available momentum
transfer that is admitted by the chosen propulsion channel.

In the mass-flow formulation, G(t) suppresses the effective contribution of
propellant exchange to dv. In the power-limited formulation, G(t) suppresses
the conversion of injected power into velocity gain. Thus, once the system
leaves the admissible phase/geometric window, neither additional propellant
exchange nor additional input power can fully compensate for the loss of
coherent coupling.

---

## 6  Falsifiability and Predictions

### 6.1  Primary Falsification Criteria

HLV–AILEE is **falsified** if:

1. no measurable correlation exists between Δφ-like metrics and coupling efficiency,
2. increased power fully compensates induced misalignment (no gate-collapse regime),
3. controlled phase drift does not produce threshold-like degradation or hysteresis.

### 6.2  Positive Predictions

Conversely, the model predicts:

1. threshold-like drops in effective performance beyond a critical window,
2. a possible "sweet spot" at reduced power if alignment improves,
3. hysteresis/path dependence under up/down sweeps of timing and geometry controls.

---

## 7  Measurement Protocol and Implementation

### 7.1  Inferring Δφ from Measurable Proxies

**Temporal proxy:**
> Δφ_temp ∼ |f_sys − f_lat| / f_lat  ················ (9)

**Spatial/vorticity proxy:**
> Δφ_spat ∼ ⟨|∇ × v_flow|⟩ / ⟨|v_flow|⟩  ············· (10)

**Spectral proxy:**
> Δφ_spec ∼ ∫ dω |S_thrust(ω) − S_lat(ω)|  ············· (11)

### 7.2  Protocol Phases

| Phase | Actions |
|-------|---------|
| **1 — Baseline** | Measure thrust, Isp, power, spectra; log geometry & environment |
| **2 — Perturbations** | Apply timing offsets; geometry modulation; induced phase drift; ignition/shutdown history |
| **3 — Analysis** | Identify thresholds, plot hysteresis loops, compare matched baselines and surrogate controls |
| **4 — Replication** | Independent labs; pre-registered metrics; public release of data and settings |

---

## 8  Figures

### Figure 1 — Conceptual Comparison

```
╔══════════════════════════════════╗   ╔══════════════════════════════════════════════╗
║       CLASSICAL MODEL            ║   ║            HLV–AILEE GATE                    ║
╠══════════════════════════════════╣   ╠══════════════════════════════════════════════╣
║                                  ║   ║                                              ║
║  Input power / mass-flow         ║   ║  System phase state  ψ_sys(t)               ║
║  P_in(t),  ṁ(t)                  ║   ║            │                                 ║
║         │                        ║   ║            ▼                                 ║
║         ▼                        ║   ║  Host phase structure  ψ_lat(t)             ║
║  Engine / nozzle efficiency η    ║   ║            │                                 ║
║  (empirical)                     ║   ║            ▼                                 ║
║         │                        ║   ║  Δφ(t) = |∂_ψ ψ_sys(t) − ψ_lat(t)|        ║
║         ▼                        ║   ║  Gate:  G(t) = exp(−α · Δφ(t)²)            ║
║  Momentum transfer  Δv           ║   ║            │                                 ║
║  (rocket eq. + losses)           ║   ║            ▼                                 ║
║         │                        ║   ║  Realizable coupling                         ║
║         ▼                        ║   ║  rocket: dv = g₀·Isp·G(t)·(−dM/M)         ║
║  Losses absorbed into            ║   ║  power:  Δv = η·∫[P_in/(M·v)]·G(t) dt     ║
║  empirical efficiency terms      ║   ║                                              ║
║                                  ║   ║  ⚠ No new force.                             ║
║                                  ║   ║    Gate encodes admissibility.               ║
╚══════════════════════════════════╝   ╚══════════════════════════════════════════════╝
```

*Conceptual comparison between a classical propulsion model and the HLV–AILEE
coherence-gated formulation. Classical descriptions absorb coupling losses into
empirical efficiency terms, whereas HLV–AILEE introduces an explicit
phase/geometric deviation variable and a bounded gate factor G(t). The gate may
be applied either to the mass-flow rocket-equation form or to the power-limited
electric/ion-propulsion form, depending on propulsion regime.*

---

### Figure 2 — Phase Window with Hysteresis

```
Effective coupling G_eff
  1.0 │╲   ←  Up-sweep (collapse branch)
      │  ╲
  0.8 │    ╲_ _
      │   ╱     ╲
  0.6 │  ╱  hyst- ╲
      │ ╱   esis   ╲
  0.4 │╱   window   ╲
      │               ╲_ _
  0.2 │   ← Down-sweep      ╲
      │     (recovery branch) ╲___
  0.0 ┼────┬────┬────┬────┬────┬────▶  Phase deviation Δφ
      0   0.5   1.0  1.5  2.0  2.5
```

*Phase-window with hysteresis: coupling collapses at smaller Δφ during an
up-sweep, but recovers only after returning deeper into the window
(down-sweep). Parameters are illustrative and should be fitted to data.*

---

## 8  Synthetic Prediction Templates for Future Experimental Tests

This section translates the proposed HLV–AILEE coherence-gate mechanism into
synthetic, model-generated prediction templates for future experimental tests.
The figures shown here are **not experimental measurements, not thrust data,
and not empirical validation of propulsion**. No physical test stand, thrust
balance, or laboratory instrument produced these curves. Their purpose is only
to define qualitative signal shapes, calibration targets, and falsification
structures that a future experimental campaign could test.

The central model assumption is that the effective coupling gate decreases when
a phase/geometric deviation proxy Δφ increases. A minimal bounded gate model is:

> **G(Δφ) = exp(−α · Δφ²)**  ·················· (12)

where α > 0 controls the sharpness of the gate collapse. Equation (12) is not
inferred from propulsion data in this work. It is used only as a model-defined
prediction template.

Figures 3–6 define four separate synthetic prediction templates:

- The **hysteresis template** tests whether collapse and recovery follow different branches under an up/down sweep of Δφ.
- The **spectral template** tests whether a localized gate response appears near a predefined frequency benchmark.
- The **calibration template** tests whether the gate sharpness parameter α can be recovered from controlled threshold-scan data.
- The **multi-proxy template** defines how temporal, spatial, and spectral deviations could be combined into a single experimental readout.

A future falsification test would require replacing the synthetic signals in
these figures with independently acquired laboratory data. The HLV–AILEE gate
mechanism would be disfavored if controlled measurements fail to show any
reproducible relationship between Δφ-type phase-deviation proxies and the
predicted gate response, or if matched null models explain the data equally
well without the proposed gate structure.

---

## Code Availability

A reference implementation of the HLV–AILEE coherence-gated propulsion model,
including the rocket-equation and power-limited variants, simulation notebooks,
and unit tests, is publicly available at:

> **https://github.com/dfeen87/Geometry-Limited-Propulsion**

The repository is provided for reproducibility and computational inspection of
the model. It does not constitute experimental validation of the HLV–AILEE
framework.

> *"No empirical propulsion data are claimed or included."*

---

## Appendix A — Summary Tables

### Table 1 — Symbol Definitions

| Symbol | Meaning |
|--------|---------|
| *I*sp | Specific impulse (standard engineering definition) |
| *η* | Aggregate classical efficiency (thermal / nozzle / conversion losses) |
| P_in(t) | Input power as a function of time |
| M(t) | Time-dependent system mass |
| v(t) | Instantaneous velocity |
| Δφ(t) | Phase/geometric deviation between system and lattice |
| *α* | Acceptance-window sharpness parameter |
| G(t) | Coherence gate  G(t) = exp(−α · Δφ(t)²) |

---

### Table 2 — Falsifiable Predictions vs. Classical Expectation

| Observable Test | Classical Expectation | HLV–AILEE Prediction |
|---|---|---|
| Thrust vs. power | Smooth saturation (empirical) | Threshold-like collapse beyond Δφ window |
| Timing offsets | Small, continuous degradation | Sharp non-linear drop; high sensitivity to phase drift |
| Ignition / shutdown history | Mostly irrelevant | Hysteresis / path dependence |
| Reduced-power regime | Lower Δv | Potential "sweet spot" under improved alignment |

---

## References

[1] M. Krüger, *A Mathematical Unification of the Helix–Light–Vortex (HLV)
Framework: Discrete Geometry, Spiral Time, Unified Lagrangians, and
TOE-Level Structure* (Zenodo, 2025).

[2] M. Krüger and D. Feeney, *Geometry-Limited Propulsion: An HLV–AILEE
Formulation of Coherence-Gated Momentum Transfer* (manuscript draft, 2025).
