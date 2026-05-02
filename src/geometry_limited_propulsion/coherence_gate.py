"""
coherence_gate.py
=================
Core HLV–AILEE formulation: coherence gate and geometry-gated velocity-gain law.

References
----------
Krüger & Feeney (2025), "Geometry-Limited Propulsion: An HLV–AILEE
Formulation of Coherence-Gated Momentum Transfer."

Equations from the paper
------------------------
Gate function (Eq. 6):
    G(t) = exp(-α · Δφ(t)²)

Velocity-gain law (Eq. 6):
    Δv = Isp · η · ∫₀^tf  [P_in(t) / (M(t) · v(t))]  · G(t)  dt
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

# Backward-compatible trapezoid integration (np.trapz removed in NumPy 2.0)
try:
    _trapezoid = np.trapezoid
except AttributeError:
    _trapezoid = np.trapz  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

G0: float = 9.80665  # standard gravity [m s⁻²]


# ---------------------------------------------------------------------------
# Gate function
# ---------------------------------------------------------------------------

def coherence_gate(delta_phi: ArrayLike, alpha: float) -> np.ndarray:
    """Compute the coherence gate G(t) = exp(-α · Δφ(t)²).

    The gate returns values in [0, 1].  At perfect alignment (Δφ = 0) it
    equals 1 (full coupling); it decays exponentially as phase/geometry
    deviation grows.

    Parameters
    ----------
    delta_phi : array-like
        Phase/geometric deviation Δφ(t) at each time step.  Dimensionless;
        see ``phase_alignment_metrics`` for proxy definitions.
    alpha : float
        Acceptance-window sharpness parameter α > 0.  Larger values yield
        a narrower gate (steeper collapse).

    Returns
    -------
    np.ndarray
        Gate values G(t) ∈ (0, 1] with the same shape as *delta_phi*.

    Raises
    ------
    ValueError
        If *alpha* is not strictly positive.

    Examples
    --------
    >>> import numpy as np
    >>> from coherence_gate import coherence_gate
    >>> phi = np.linspace(0, 2, 9)
    >>> G = coherence_gate(phi, alpha=2.0)
    >>> G.round(4)
    array([1.    , 0.9394, 0.7788, 0.5698, 0.3679, 0.2096, 0.1054, 0.0467, 0.0183])
    """
    if alpha <= 0:
        raise ValueError(f"alpha must be strictly positive, got {alpha!r}.")
    delta_phi = np.asarray(delta_phi, dtype=float)
    return np.exp(-alpha * delta_phi**2)


# ---------------------------------------------------------------------------
# Classical (Tsiolkovsky) Δv — baseline comparison
# ---------------------------------------------------------------------------

def classical_delta_v(
    isp: float,
    mass_initial: float,
    mass_final: float,
) -> float:
    """Tsiolkovsky rocket equation (Eq. 1 of paper).

    Δv = Isp · g₀ · ln(M₀ / Mf)

    Parameters
    ----------
    isp : float
        Specific impulse [s].
    mass_initial : float
        Initial mass M₀ [kg].
    mass_final : float
        Final mass Mf [kg].

    Returns
    -------
    float
        Classical velocity gain Δv [m s⁻¹].
    """
    if isp <= 0:
        raise ValueError(f"isp must be strictly positive, got {isp!r}.")
    if mass_final <= 0 or mass_initial <= mass_final:
        raise ValueError(
            "Require mass_initial > mass_final > 0; "
            f"got M₀={mass_initial}, Mf={mass_final}."
        )
    return isp * G0 * np.log(mass_initial / mass_final)


# ---------------------------------------------------------------------------
# HLV–AILEE gated velocity gain (numerical integration)
# ---------------------------------------------------------------------------

# Valid mode identifiers for hlv_ailee_delta_v()
MODES = frozenset({"power", "rocket"})


def hlv_ailee_delta_v(
    t: ArrayLike,
    mass: ArrayLike,
    delta_phi: ArrayLike,
    isp: float,
    eta: float,
    alpha: float,
    *,
    mode: str = "power",
    # power-mode arguments
    p_in: ArrayLike | None = None,
    velocity: ArrayLike | None = None,
    # rocket-mode arguments
    m_dot: ArrayLike | None = None,
) -> float:
    """Geometry-gated velocity gain with selectable physical framing.

    Two dimensionally consistent formulations are supported via the *mode*
    toggle, addressing a design question raised during model review:

    **mode='power'** — Power-limited framing (original formulation, default).
    Best suited to electric / ion propulsion where thrust derives from input
    power rather than direct propellant expulsion:

        Δv = Isp · η · ∫₀^tf  [P_in(t) / (M(t) · v(t))]  · G(t)  dt

    **mode='rocket'** — Rocket-equation framing (Marcel Krüger's preferred
    formulation for chemical propulsion). Places the gate directly on each
    mass-expulsion event, recovering Tsiolkovsky exactly when G = 1:

        Δv = g₀ · Isp · ∫₀^tf  G(t) · (−Ṁ(t) / M(t))  dt

    In both cases G(t) = exp(−α · Δφ(t)²) is the coherence gate.

    Choosing between modes
    ----------------------
    * Chemical / bipropellant / solid rocket → ``mode='rocket'``
      Supply *m_dot*; *p_in* and *velocity* are unused.
    * Electric / ion / Hall thruster        → ``mode='power'``
      Supply *p_in* and *velocity*; *m_dot* is unused.

    Parameters shared by both modes
    --------------------------------
    t : array-like, shape (N,)
        Time grid [s].  Need not be uniform; trapezoid rule is used.
    mass : array-like, shape (N,)
        System mass M(t) [kg].  Must be strictly positive everywhere.
    delta_phi : array-like, shape (N,)
        Phase/geometric deviation Δφ(t).  Dimensionless.
    isp : float
        Specific impulse [s].
    eta : float
        Aggregate classical efficiency η ∈ (0, 1].
        Used only in ``mode='power'``; in ``mode='rocket'`` efficiency
        is already encoded in Isp and Ṁ.
    alpha : float
        Gate sharpness parameter α > 0.
    mode : str
        ``'power'`` (default) or ``'rocket'``.

    Power-mode only (mode='power')
    --------------------------------
    p_in : array-like, shape (N,)
        Input power P_in(t) [W].
    velocity : array-like, shape (N,)
        Instantaneous velocity v(t) [m s⁻¹].  Must be strictly positive.

    Rocket-mode only (mode='rocket')
    ---------------------------------
    m_dot : array-like, shape (N,)
        Mass-flow rate Ṁ(t) [kg s⁻¹].  Positive values indicate propellant
        consumption (i.e. dM/dt = −Ṁ, so Ṁ > 0 means mass is decreasing).

    Returns
    -------
    float
        Gated velocity gain Δv [m s⁻¹].

    Raises
    ------
    ValueError
        For unrecognised *mode*, missing required arguments, or invalid inputs.

    Examples
    --------
    >>> import numpy as np
    >>> N = 200
    >>> t = np.linspace(0, 60, N)
    >>> mass = np.linspace(800, 500, N)

    Power-limited (electric propulsion):

    >>> dv = hlv_ailee_delta_v(
    ...     t, mass, delta_phi=np.full(N, 0.2),
    ...     isp=3000.0, eta=0.65, alpha=2.0,
    ...     mode='power',
    ...     p_in=np.full(N, 2e5),
    ...     velocity=np.linspace(200, 800, N),
    ... )

    Rocket-equation (chemical propulsion):

    >>> dv = hlv_ailee_delta_v(
    ...     t, mass, delta_phi=np.full(N, 0.2),
    ...     isp=350.0, eta=0.85, alpha=2.0,
    ...     mode='rocket',
    ...     m_dot=np.full(N, 5.0),
    ... )

    Notes
    -----
    * No new forces are introduced.  The gate encodes *admissibility* of
      coupling, not an additional thrust mechanism.
    * When Δφ ≡ 0, mode='rocket' recovers the Tsiolkovsky equation exactly
      (up to numerical integration error).
    * When Δφ ≡ 0, mode='power' reduces to a power-weighted velocity integral.
    """
    if mode not in MODES:
        raise ValueError(
            f"mode must be one of {sorted(MODES)!r}, got {mode!r}."
        )
    if isp <= 0:
        raise ValueError(f"isp must be strictly positive, got {isp!r}.")

    t    = np.asarray(t,         dtype=float)
    mass = np.asarray(mass,      dtype=float)
    dphi = np.asarray(delta_phi, dtype=float)

    if not (t.shape == mass.shape == dphi.shape):
        raise ValueError("t, mass, and delta_phi must have the same shape.")
    if np.any(mass <= 0):
        raise ValueError("mass must be strictly positive at all time steps.")
    if not (0 < eta <= 1):
        raise ValueError(f"eta must be in (0, 1], got {eta!r}.")

    gate = coherence_gate(dphi, alpha)

    # ------------------------------------------------------------------
    if mode == "power":
        if p_in is None or velocity is None:
            raise ValueError(
                "mode='power' requires both p_in and velocity arrays."
            )
        p_in_arr = np.asarray(p_in,     dtype=float)
        vel_arr  = np.asarray(velocity, dtype=float)
        if not (t.shape == p_in_arr.shape == vel_arr.shape):
            raise ValueError("p_in and velocity must match the shape of t.")
        if np.any(vel_arr <= 0):
            raise ValueError("velocity must be strictly positive at all time steps.")
        integrand = (p_in_arr / (mass * vel_arr)) * gate
        return float(isp * eta * _trapezoid(integrand, t))

    # ------------------------------------------------------------------
    else:  # mode == "rocket"
        if m_dot is None:
            raise ValueError(
                "mode='rocket' requires m_dot (mass-flow rate) array."
            )
        m_dot_arr = np.asarray(m_dot, dtype=float)
        if not (t.shape == m_dot_arr.shape):
            raise ValueError("m_dot must match the shape of t.")
        if np.any(m_dot_arr < 0):
            raise ValueError(
                "m_dot must be non-negative (positive = propellant consumed)."
            )
        # dv = g₀ · Isp · G(t) · (−dM/M) = g₀ · Isp · G(t) · (Ṁ/M) dt
        integrand = gate * (m_dot_arr / mass)
        return float(G0 * isp * _trapezoid(integrand, t))


# ---------------------------------------------------------------------------
# Sweep helpers for prediction / falsifiability analysis
# ---------------------------------------------------------------------------

def gate_vs_phase_sweep(
    delta_phi_values: ArrayLike,
    alpha: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the gate across a range of Δφ values.

    Useful for plotting the acceptance-window shape (Figure 2 of paper).

    Parameters
    ----------
    delta_phi_values : array-like
        1-D array of Δφ values to sweep.
    alpha : float
        Gate sharpness parameter.

    Returns
    -------
    delta_phi : np.ndarray
        Echo of the input sweep values.
    gate : np.ndarray
        Corresponding gate values G = exp(-α · Δφ²).
    """
    delta_phi_values = np.asarray(delta_phi_values, dtype=float)
    return delta_phi_values, coherence_gate(delta_phi_values, alpha)


def hysteresis_sweep(
    delta_phi_up: ArrayLike,
    delta_phi_down: ArrayLike,
    alpha_up: float,
    alpha_down: float,
) -> dict[str, np.ndarray]:
    """Simulate a simple hysteresis loop (Figure 2 of paper).

    In the paper, up-sweeps (increasing Δφ) collapse coupling earlier than
    down-sweeps (decreasing Δφ) recover it.  This is modelled here by
    allowing different effective sharpness values for the two branches —
    a minimal phenomenological stand-in until a dynamical hysteresis model
    is fitted to data.

    Parameters
    ----------
    delta_phi_up : array-like
        Δφ values for the up-sweep branch.
    delta_phi_down : array-like
        Δφ values for the down-sweep branch.
    alpha_up : float
        Gate sharpness for the up-sweep (collapse branch).
    alpha_down : float
        Gate sharpness for the down-sweep (recovery branch); typically
        alpha_down < alpha_up to reflect delayed recovery.

    Returns
    -------
    dict with keys
        ``"phi_up"``, ``"gate_up"`` — up-sweep branch.
        ``"phi_down"``, ``"gate_down"`` — down-sweep branch.
    """
    phi_up = np.asarray(delta_phi_up, dtype=float)
    phi_down = np.asarray(delta_phi_down, dtype=float)
    return {
        "phi_up": phi_up,
        "gate_up": coherence_gate(phi_up, alpha_up),
        "phi_down": phi_down,
        "gate_down": coherence_gate(phi_down, alpha_down),
    }
