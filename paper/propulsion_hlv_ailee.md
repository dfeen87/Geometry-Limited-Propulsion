# Geometry-Limited Propulsion:
## An HLV–AILEE Formulation of Coherence-Gated Momentum Transfer

**Marcel Krüger**
Independent Researcher, HLV–RFP (Meiningen, Germany)
marcelkrueger092@gmail.com · ORCID: 0009-0002-5709-9729

**Don Feeney**
Independent Researcher (Pennsylvania, USA)
dfeen87@gmail.com · ORCID: 0009-0003-1350-4160

*December 26, 2025*

---

## Abstract

Classical propulsion models encode efficiency losses through empirical terms, implicitly assuming near-optimal coupling between injected power and momentum transfer. Yet many systems display sharp saturation behavior where additional energy fails to produce proportional velocity gain. We propose that such saturation reflects a deeper constraint: phase and geometric alignment between propulsion dynamics and an admissible interaction structure of an emergent substrate. Using the Helix–Light–Vortex framework — in which spacetime is modeled as a discrete golden-ratio quasicrystal lattice governed by triadic spiral time — we formalize a coherence gate that exponentially suppresses energy–momentum coupling under phase deviation. The resulting HLV–AILEE velocity-gain law preserves conservation principles and makes no reactionless or superluminal claims. Instead, it yields falsifiable predictions:

1. threshold-like collapse of effective performance beyond a phase-alignment window,
2. hysteresis/path dependence under timing and geometric drift, and
3. counterintuitive optimal regimes at reduced power under improved coherence.

We outline measurement protocols and computational tests accessible with existing instrumentation.

---

## 1  Introduction: Saturation in Classical Propulsion

Classical propulsion formulations typically assume that increasing input power and improving engineering efficiency lead to proportional gains in thrust or velocity, modulo known loss channels (thermal, turbulent, nozzle, and material constraints). In practice, many systems exhibit sharp saturation: beyond an operating regime, additional injected energy fails to yield proportional momentum gain and instead dissipates into heat, noise, or unstructured flow.

This paper proposes that such saturation can be reinterpreted as a **constraint problem**: momentum transfer is realizable only when the propulsion dynamics remain sufficiently aligned with an admissible phase/geometry window. We formalize this within the Helix–Light–Vortex framework and derive a gate-controlled coupling law (HLV–AILEE) that is empirically distinguishable from generic "efficiency tuning."

---

## 2  The Missing Variable: Geometry and Phase Alignment

### 2.1  From Empirical Efficiency to an Explicit Coupling Constraint

Consider the standard rocket equation:

> **Δv = Isp · g₀ · ln(M₀ / Mf)**  ·················· (1)

Here *I*sp is an engineering parameter. In the HLV–AILEE view, realizable coupling also depends on whether the system remains within an admissible phase window relative to the host geometry. This motivates introducing a first-class alignment variable rather than absorbing everything into *η*.

---

## 3  HLV Host Geometry (Minimal Ingredients)

### 3.1  Discrete Quasicrystal Host (𝒢_φ)

Within HLV, physical spacetime is modeled as a 3D projection of a higher-dimensional quasicrystalline lattice (cut-and-project). We denote the host schematically by 𝒢_φ with golden ratio φ = (1 + √5) / 2 and acceptance window W_φ:

> **𝒢_φ = { π‖(x) ∈ E‖ : x ∈ ℤ⁶, π⊥(x) ∈ W_φ }**  ·············· (2)

### 3.2  Triadic Spiral Time

HLV extends time to a triadic operator:

> **ψ(t) = t + i·φ(t) + j·χ(t)**  ·················· (3)

where *t* is the U₁ component, φ(t) a coherence mode (U₂), and χ(t) a memory/retrocausal mode (U₃). For propulsion modeling we only require that this structure can induce small controlled deformations of propagation:

> **A(t) = 1 + ε(t),   |ε(t)| ≪ 1**  ················ (4)

without violating conservation laws.

---

## 4  Phase Deviation as a Control Variable Δφ

### 4.1  Definition (Operational)

We define a phase/geometric deviation metric that compares system dynamics to the host admissible structure:

> **Δφ(t) = |∂_ψ ψ_sys(t) − ψ_lat(t)|**  ················ (5)

### 4.2  Coherence-Gate Behavior

Δφ does not introduce a new force. It gates coupling efficiency:

- **Δφ ≃ 0** → coherent coupling into momentum transfer
- **Δφ outside acceptance window** → rapid collapse of realizable coupling

---

## 5  HLV–AILEE Law: Geometry-Gated Velocity Gain

### 5.1  Core Equation

We encode the constraint by an exponential coherence gate:

> **Δv = Isp · η · ∫₀^tf [ P_in(t) / (M(t)·v(t)) ] · exp(−α·Δφ(t)²) dt**  ···· (6)

### 5.2  Interpretation

Increasing input power **cannot** compensate for large phase/geometry mismatch once the gate suppresses coupling.

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
> Δφ_temp ∼ |f_sys − f_lat| / f_lat  ················ (7)

**Spatial/vorticity proxy:**
> Δφ_spat ∼ ⟨|∇ × v_flow|⟩ / ⟨|v_flow|⟩  ············· (8)

**Spectral proxy:**
> Δφ_spec ∼ ∫ dω |S_thrust(ω) − S_lat(ω)|  ············· (9)

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
╔══════════════════════════════════╗   ╔══════════════════════════════════════╗
║       CLASSICAL MODEL            ║   ║         HLV–AILEE GATE               ║
╠══════════════════════════════════╣   ╠══════════════════════════════════════╣
║                                  ║   ║                                      ║
║  Input power / mass-flow         ║   ║  System phase state  ψ_sys(t)        ║
║  P_in(t),  ṁ(t)                  ║   ║            │                         ║
║         │                        ║   ║            ▼                         ║
║         ▼                        ║   ║  Host phase structure  ψ_lat(t)      ║
║  Engine / nozzle efficiency η    ║   ║            │                         ║
║  (empirical)                     ║   ║            ▼                         ║
║         │                        ║   ║  Δφ(t) = |∂_ψ ψ_sys − ψ_lat|       ║
║         ▼                        ║   ║            │                         ║
║  Momentum transfer  Δv           ║   ║            ▼                         ║
║  (rocket eq. + losses)           ║   ║  Gate  G(t) = exp(−α·Δφ²)           ║
║         │                        ║   ║            │                         ║
║         ▼                        ║   ║            ▼                         ║
║  Losses absorbed into            ║   ║  Realizable coupling                 ║
║  empirical efficiency terms      ║   ║  Δv ∝ ∫ ··· G(t) dt                 ║
║                                  ║   ║                                      ║
║                                  ║   ║  ⚠ No new force.                     ║
║                                  ║   ║    Gate encodes admissibility.       ║
╚══════════════════════════════════╝   ╚══════════════════════════════════════╝
```

*Classical propulsion absorbs losses into empirical efficiency. HLV–AILEE introduces an alignment variable and a coherence gate that suppresses realizable coupling outside an acceptance window.*

---

### Figure 2 — Phase Window with Hysteresis

```
Effective coupling G_eff
  1.0 │╲
      │  ╲   ← Down-sweep (recovers later)
  0.8 │    ╲
      │      ╲
  0.6 │        ╲_ _ _
      │        ╱      ╲
  0.4 │      ╱   hysteresis╲
      │    ╱      window    ╲
  0.2 │  ╱  ← Up-sweep       ╲
      │╱   (collapses earlier)  ╲___
  0.0 ┼────┬────┬────┬────┬────┬────▶  Phase deviation Δφ
      0   0.5   1.0  1.5  2.0  2.5
```

*Coupling collapses at smaller Δφ during an up-sweep, but recovers only after returning deeper into the window (down-sweep). Parameters are illustrative and should be fitted to data.*

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

[1] M. Krüger, *A Mathematical Unification of the Helix–Light–Vortex (HLV) Framework: Discrete Geometry, Spiral Time, Unified Lagrangians, and TOE-Level Structure* (Zenodo, 2025).

[2] M. Krüger and D. Feeney, *Geometry-Limited Propulsion: An HLV–AILEE Formulation of Coherence-Gated Momentum Transfer* (manuscript draft, 2025).
