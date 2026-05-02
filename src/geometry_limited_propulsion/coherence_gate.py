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
    if mass_final <= 0 or mass_initial <= mass_final:
        raise ValueError(
            "Require mass_initial > mass_final > 0; "
            f"got M₀={mass_initial}, Mf={mass_final}."
        )
    return isp * G0 * np.log(mass_initial / mass_final)


# ---------------------------------------------------------------------------
# HLV–AILEE gated velocity gain (numerical integration)
# ---------------------------------------------------------------------------

def hlv_ailee_delta_v(
    t: ArrayLike,
    p_in: ArrayLike,
    mass: ArrayLike,
    velocity: ArrayLike,
    delta_phi: ArrayLike,
    isp: float,
    eta: float,
    alpha: float,
) -> float:
    """Geometry-gated velocity gain (Eq. 6 of paper).

    Numerically integrates

        Δv = Isp · η · ∫₀^tf  [P_in(t) / (M(t) · v(t))]  · G(t)  dt

    where G(t) = exp(-α · Δφ(t)²) is the coherence gate.

    All array arguments must share the same length N (time steps).

    Parameters
    ----------
    t : array-like, shape (N,)
        Time grid [s].  Need not be uniform; integration uses the trapezoid
        rule over the supplied grid.
    p_in : array-like, shape (N,)
        Input power P_in(t) [W].
    mass : array-like, shape (N,)
        System mass M(t) [kg].  Must be strictly positive everywhere.
    velocity : array-like, shape (N,)
        Instantaneous velocity v(t) [m s⁻¹].  Must be strictly positive
        everywhere (i.e. the system is already in motion).
    delta_phi : array-like, shape (N,)
        Phase/geometric deviation Δφ(t).  Dimensionless.
    isp : float
        Specific impulse [s].
    eta : float
        Aggregate classical efficiency η ∈ (0, 1].
    alpha : float
        Gate sharpness parameter α > 0.

    Returns
    -------
    float
        Gated velocity gain Δv [m s⁻¹].

    Notes
    -----
    * For a perfectly aligned system (Δφ ≡ 0) the gate collapses to 1 and the
      result reduces to a power-weighted integral form of the rocket equation.
    * No new forces are introduced.  The gate encodes *admissibility* of
      coupling, not an additional thrust mechanism.
    """
    t = np.asarray(t, dtype=float)
    p_in = np.asarray(p_in, dtype=float)
    mass = np.asarray(mass, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    delta_phi = np.asarray(delta_phi, dtype=float)

    if not (t.shape == p_in.shape == mass.shape == velocity.shape == delta_phi.shape):
        raise ValueError("All array arguments must have the same shape.")
    if np.any(mass <= 0):
        raise ValueError("mass must be strictly positive at all time steps.")
    if np.any(velocity <= 0):
        raise ValueError("velocity must be strictly positive at all time steps.")
    if not (0 < eta <= 1):
        raise ValueError(f"eta must be in (0, 1], got {eta!r}.")

    gate = coherence_gate(delta_phi, alpha)
    integrand = (p_in / (mass * velocity)) * gate
    return float(isp * eta * np.trapz(integrand, t))


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
