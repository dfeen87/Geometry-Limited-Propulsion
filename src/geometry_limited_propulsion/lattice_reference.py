# MIT License
# Copyright (c) 2025 Marcel Krüger, Don Feeney
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""
lattice_reference.py
====================
Derives lattice reference quantities directly from the HLV golden-ratio
quasicrystal geometry rather than treating them as free user parameters.

This module closes the gap between the measurable Δφ proxies in
``phase_alignment_metrics.py`` and the underlying HLV host-lattice structure.
Previously, quantities like ``f_lat`` and ``ψ_lat`` had to be chosen
arbitrarily.  Here they are computed from first principles via the φ-harmonic
frequency series implied by the quasicrystal acceptance window W_φ.

The derivation chain this module makes explicit:

    HLV geometry (φ)  →  φ-harmonic series  →  f_lat
                                             →  ψ_lat(t)
    Empirical threshold data                 →  α  (via AlphaCalibrator)

With this module imported, the full pipeline becomes deterministic:

    HLV geometry (φ) → f_lat → ψ_lat → Δφ proxies → G(t) → Δv

References
----------
Krüger (2025), "A Mathematical Unification of the Helix–Light–Vortex (HLV)
Framework" (Zenodo).

Krüger & Feeney (2025), "Geometry-Limited Propulsion: An HLV–AILEE
Formulation of Coherence-Gated Momentum Transfer."

Notes on the quasicrystal structure
------------------------------------
The HLV host lattice 𝒢_φ is defined (Eq. 2 of the paper) as a cut-and-project
of ℤ⁶ onto a 3-D physical subspace E‖, gated by an acceptance window W_φ
whose width is governed by the golden ratio φ = (1+√5)/2.

A key property of Penrose / icosahedral quasicrystals is that their
diffraction spectrum is *discrete but dense*: admissible wavevectors form the
set { m + n·φ : m,n ∈ ℤ } (the ring ℤ[φ]).  Translated to frequencies, the
admissible lattice frequencies are

    f_lat(m, n) = f₀ · φⁿ · m,    m, n ∈ ℤ,  m ≥ 1

where f₀ is a base (ground-mode) frequency fixed by the physical scale of the
system.  For a given system frequency f_sys the nearest admissible frequency
is found by minimising |f_sys − f_lat(m, n)| over a finite search grid.

The triadic time operator ψ(t) = t + i·φ(t) + j·χ(t) (Eq. 3) contributes
two oscillatory modes beyond the U₁ real-time component.  For the reference
signal ψ_lat we model the U₂ coherence mode as a φ-modulated cosine and the
U₃ memory mode as a sub-harmonic at f₀/φ, both with small amplitude ε ≪ 1
consistent with Eq. 4.
"""

from __future__ import annotations

import warnings
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike


# ---------------------------------------------------------------------------
# Module-level constant
# ---------------------------------------------------------------------------

PHI: float = (1.0 + np.sqrt(5.0)) / 2.0  # golden ratio φ ≈ 1.6180339887


# ---------------------------------------------------------------------------
# φ-harmonic frequency series
# ---------------------------------------------------------------------------

def phi_harmonic_series(
    f0: float,
    n_range: tuple[int, int] = (-4, 4),
    m_range: tuple[int, int] = (1, 8),
) -> np.ndarray:
    """Return the φ-harmonic frequency lattice for a given base frequency.

    Admissible lattice frequencies are

        f(m, n) = f₀ · φⁿ · m,    m ∈ [m_min, m_max],  n ∈ [n_min, n_max]

    arising from the ℤ[φ] structure of the quasicrystal diffraction spectrum.

    Parameters
    ----------
    f0 : float
        Base (ground-mode) frequency [Hz].  Sets the physical scale; in
        practice this might be the dominant combustion or oscillation
        frequency of the propulsion system at a reference operating point.
    n_range : tuple (n_min, n_max)
        Inclusive range of the φ-exponent n.  Negative n gives sub-harmonics
        (f < f₀), positive n gives super-harmonics.
    m_range : tuple (m_min, m_max)
        Inclusive range of the integer multiplier m ≥ 1.

    Returns
    -------
    np.ndarray
        1-D array of admissible frequencies, sorted ascending and deduplicated
        to a relative tolerance of 1 × 10⁻⁹.

    Examples
    --------
    >>> series = phi_harmonic_series(f0=100.0, n_range=(-1, 1), m_range=(1, 2))
    >>> series.round(4)
    array([ 61.8034,  100.    , 123.6068,  161.8034,  200.    ,  323.6068])
    """
    if f0 <= 0:
        raise ValueError(f"f0 must be strictly positive, got {f0!r}.")

    n_min, n_max = n_range
    m_min, m_max = m_range
    if m_min < 1:
        raise ValueError("m_min must be ≥ 1 (m is a positive integer multiplier).")

    freqs = []
    for n in range(n_min, n_max + 1):
        for m in range(m_min, m_max + 1):
            freqs.append(f0 * (PHI**n) * m)

    freqs_arr = np.array(sorted(set(freqs)))

    # Deduplicate within floating-point tolerance
    if len(freqs_arr) > 1:
        gaps = np.diff(freqs_arr)
        mask = np.concatenate(([True], gaps > 1e-9 * freqs_arr[1:]))
        freqs_arr = freqs_arr[mask]

    return freqs_arr


def nearest_lattice_frequency(
    f_sys: float | ArrayLike,
    f0: float,
    n_range: tuple[int, int] = (-4, 4),
    m_range: tuple[int, int] = (1, 8),
) -> tuple[np.ndarray, np.ndarray]:
    """Find the nearest φ-harmonic lattice frequency for each system frequency.

    For each value in *f_sys* the function searches the φ-harmonic series
    (see :func:`phi_harmonic_series`) and returns the closest admissible
    frequency.  This replaces the ad-hoc choice of ``f_lat`` in
    ``phase_alignment_metrics.delta_phi_temporal``.

    Parameters
    ----------
    f_sys : float or array-like, shape (N,)
        System frequency or array of system frequencies [Hz].
    f0 : float
        Base frequency for the φ-harmonic series [Hz].
    n_range, m_range :
        Passed through to :func:`phi_harmonic_series`.

    Returns
    -------
    f_lat : np.ndarray
        Nearest admissible lattice frequency for each element of *f_sys*,
        same shape as *f_sys*.
    delta_phi_temp : np.ndarray
        Corresponding temporal proxy Δφ_temp = |f_sys − f_lat| / f_lat,
        derived from lattice-grounded f_lat rather than a free parameter.

    Examples
    --------
    >>> f_lat, dphi = nearest_lattice_frequency(105.0, f0=100.0)
    >>> print(f"f_lat = {f_lat[0]:.4f} Hz,  Δφ_temp = {dphi[0]:.4f}")
    f_lat = 100.0000 Hz,  Δφ_temp = 0.0500
    """
    series = phi_harmonic_series(f0, n_range=n_range, m_range=m_range)
    f_sys_arr = np.atleast_1d(np.asarray(f_sys, dtype=float))

    # Vectorised nearest-neighbour search over the series
    f_lat = np.empty_like(f_sys_arr)
    for i, fs in enumerate(f_sys_arr.flat):
        idx = np.argmin(np.abs(series - fs))
        f_lat.flat[i] = series[idx]

    delta_phi_temp = np.abs(f_sys_arr - f_lat) / f_lat
    return f_lat, delta_phi_temp


# ---------------------------------------------------------------------------
# Lattice reference signal  ψ_lat(t)
# ---------------------------------------------------------------------------

def psi_lat(
    t: ArrayLike,
    f0: float,
    epsilon: float = 0.05,
) -> dict[str, np.ndarray]:
    """Generate the triadic lattice reference signal ψ_lat(t).

    Models the three components of the HLV triadic time operator (Eq. 3):

        ψ(t) = t  +  i·φ(t)  +  j·χ(t)

    as concrete oscillatory signals:

    * **U₁ — real time** ``u1``: a unit-amplitude cosine at f₀ (the
      ground-mode lattice frequency).
    * **U₂ — coherence mode** ``u2``: a cosine at f₀·φ, scaled by ε.
      The φ-shift places it at the nearest super-harmonic in ℤ[φ].
    * **U₃ — memory/retrocausal mode** ``u3``: a cosine at f₀/φ, scaled
      by ε.  The sub-harmonic captures the retrocausal texture of the U₃
      mode while keeping |ε| ≪ 1 (Eq. 4).

    The composite signal

        s_lat(t) = u1 + ε·u2 + ε·u3

    serves as the reference waveform for ``delta_phi_spectral`` in
    ``phase_alignment_metrics.py``, replacing the user-supplied
    ``lat_signal`` with a theoretically grounded one.

    Parameters
    ----------
    t : array-like, shape (N,)
        Time grid [s].
    f0 : float
        Ground-mode lattice frequency [Hz].
    epsilon : float
        Amplitude of the U₂ and U₃ perturbation modes.  Must satisfy
        |ε| ≪ 1 (paper Eq. 4).  Default 0.05.

    Returns
    -------
    dict with keys
        ``"t"``       — time array (echo of input).
        ``"u1"``      — U₁ real-time component.
        ``"u2"``      — U₂ coherence mode (φ·f₀).
        ``"u3"``      — U₃ memory mode (f₀/φ).
        ``"signal"``  — composite s_lat(t) = u1 + ε·(u2 + u3).
        ``"f0"``      — ground-mode frequency used.
        ``"f_phi"``   — U₂ frequency f₀·φ.
        ``"f_sub"``   — U₃ frequency f₀/φ.

    Raises
    ------
    ValueError
        If epsilon ≥ 0.5 (violates the |ε| ≪ 1 constraint of Eq. 4).

    Examples
    --------
    >>> import numpy as np
    >>> t = np.linspace(0, 1, 1024)
    >>> ref = psi_lat(t, f0=100.0, epsilon=0.05)
    >>> ref['signal'].shape
    (1024,)
    >>> round(ref['f_phi'], 4)
    161.8034
    """
    if abs(epsilon) >= 0.5:
        raise ValueError(
            f"epsilon={epsilon!r} violates |ε| ≪ 1 (Eq. 4).  Use a value < 0.5."
        )
    if epsilon >= 0.1:
        warnings.warn(
            f"epsilon={epsilon!r} is approaching the ≪ 1 boundary; "
            "consider using a smaller value (e.g. 0.01–0.05).",
            UserWarning,
            stacklevel=2,
        )

    t_arr = np.asarray(t, dtype=float)
    f_phi = f0 * PHI          # U₂ super-harmonic: f₀·φ
    f_sub = f0 / PHI          # U₃ sub-harmonic:   f₀/φ

    u1 = np.cos(2 * np.pi * f0    * t_arr)   # real-time component
    u2 = np.cos(2 * np.pi * f_phi * t_arr)   # coherence mode
    u3 = np.cos(2 * np.pi * f_sub * t_arr)   # memory mode

    signal = u1 + epsilon * (u2 + u3)

    return {
        "t":      t_arr,
        "u1":     u1,
        "u2":     u2,
        "u3":     u3,
        "signal": signal,
        "f0":     f0,
        "f_phi":  f_phi,
        "f_sub":  f_sub,
    }


# ---------------------------------------------------------------------------
# Acceptance-window width from quasicrystal geometry
# ---------------------------------------------------------------------------

def acceptance_window_width(n_range: tuple[int, int] = (-4, 4)) -> float:
    """Estimate the natural acceptance-window half-width from W_φ geometry.

    The quasicrystal acceptance window W_φ (Eq. 2) has a width in the
    perpendicular space E⊥ proportional to 1/φ² (the square of the reciprocal
    golden ratio).  Projected onto the frequency axis this implies a
    dimensionless half-width

        Δφ_window ≈ 1 / φ^|n_max|

    where |n_max| is the largest φ-exponent used in the harmonic series.
    This provides a geometry-derived estimate for the acceptance window
    boundary beyond which the gate G(t) is expected to collapse.

    In practice this should be treated as a *prior* estimate; the
    :class:`AlphaCalibrator` refines it from empirical threshold data.

    Parameters
    ----------
    n_range : tuple (n_min, n_max)
        Range of φ-exponents used in the harmonic series.

    Returns
    -------
    float
        Estimated dimensionless half-width of the acceptance window.

    Examples
    --------
    >>> acceptance_window_width(n_range=(-4, 4))
    0.14589803375031546
    """
    n_max = max(abs(n_range[0]), abs(n_range[1]))
    return float(1.0 / PHI**n_max)


def alpha_from_window_width(window_half_width: float) -> float:
    """Convert an acceptance-window half-width to the gate sharpness α.

    By convention we define the half-width as the Δφ at which G drops to
    e⁻¹ ≈ 0.368 (the 1/e point):

        G(Δφ_hw) = exp(−α · Δφ_hw²) = e⁻¹
        ⟹  α = 1 / Δφ_hw²

    This gives a geometry-derived prior for α rather than an arbitrary choice.

    Parameters
    ----------
    window_half_width : float
        Acceptance-window half-width Δφ at which coupling drops to e⁻¹.

    Returns
    -------
    float
        Corresponding gate sharpness parameter α.

    Examples
    --------
    >>> hw = acceptance_window_width(n_range=(-4, 4))
    >>> alpha_from_window_width(hw)
    46.97871376374779
    """
    if window_half_width <= 0:
        raise ValueError("window_half_width must be strictly positive.")
    return float(1.0 / window_half_width**2)


# ---------------------------------------------------------------------------
# AlphaCalibrator — fit α from empirical threshold data
# ---------------------------------------------------------------------------

class AlphaCalibrator:
    """Fit the gate sharpness parameter α from empirical coupling-efficiency data.

    Usage
    -----
    1. Collect pairs of (Δφ, G_measured) from a thrust/coupling experiment.
    2. Call :meth:`fit` to estimate α via nonlinear least squares.
    3. Use :attr:`alpha_fitted` in ``coherence_gate()`` calls.

    The fitting minimises

        Σᵢ  [G_measured(i)  −  exp(−α · Δφ(i)²)]²

    over α > 0.  A geometry-derived prior from :func:`alpha_from_window_width`
    is used as the initial guess if not otherwise specified.

    Parameters
    ----------
    f0 : float, optional
        Base frequency [Hz].  Used to compute the geometry-derived α prior.
        If None, the prior defaults to α = 2.0.
    n_range : tuple, optional
        φ-exponent range for the prior computation.

    Attributes
    ----------
    alpha_prior : float
        Geometry-derived prior estimate of α (before fitting).
    alpha_fitted : float or None
        Fitted α after calling :meth:`fit`.  None until fitted.
    residuals : np.ndarray or None
        Residuals at the fitted α.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> true_alpha = 3.0
    >>> phi_vals = np.linspace(0, 1.5, 30)
    >>> G_obs = np.exp(-true_alpha * phi_vals**2) + rng.normal(0, 0.02, 30)
    >>> cal = AlphaCalibrator(f0=100.0)
    >>> cal.fit(phi_vals, G_obs)
    >>> round(cal.alpha_fitted, 2)  # doctest: +SKIP
    3.0
    """

    def __init__(
        self,
        f0: float | None = None,
        n_range: tuple[int, int] = (-4, 4),
    ) -> None:
        self.f0 = f0
        self.n_range = n_range
        self.alpha_fitted: float | None = None
        self.residuals: np.ndarray | None = None

        if f0 is not None:
            hw = acceptance_window_width(n_range)
            self.alpha_prior = alpha_from_window_width(hw)
        else:
            self.alpha_prior = 2.0  # neutral fallback

    def fit(
        self,
        delta_phi: ArrayLike,
        G_measured: ArrayLike,
        alpha_init: float | None = None,
    ) -> "AlphaCalibrator":
        """Fit α to measured (Δφ, G) pairs.

        Parameters
        ----------
        delta_phi : array-like, shape (N,)
            Measured phase-deviation values Δφ.
        G_measured : array-like, shape (N,)
            Measured coupling-efficiency values G ∈ [0, 1].
        alpha_init : float, optional
            Initial guess for α.  Defaults to :attr:`alpha_prior`.

        Returns
        -------
        self
            Returns the calibrator instance for method chaining.

        Raises
        ------
        ImportError
            If SciPy is not installed.
        RuntimeError
            If the optimiser fails to converge.
        """
        try:
            from scipy.optimize import curve_fit
        except ImportError as exc:
            raise ImportError(
                "SciPy is required for AlphaCalibrator.fit(). "
                "Install it with: pip install scipy"
            ) from exc

        delta_phi = np.asarray(delta_phi, dtype=float)
        G_measured = np.asarray(G_measured, dtype=float)

        if delta_phi.shape != G_measured.shape:
            raise ValueError("delta_phi and G_measured must have the same shape.")

        G_measured = np.clip(G_measured, 0.0, 1.0)

        def gate_model(phi, alpha):
            return np.exp(-alpha * phi**2)

        p0 = [alpha_init if alpha_init is not None else self.alpha_prior]

        try:
            popt, _ = curve_fit(
                gate_model,
                delta_phi,
                G_measured,
                p0=p0,
                bounds=(1e-6, np.inf),
                maxfev=10_000,
            )
        except RuntimeError as exc:
            raise RuntimeError(
                f"AlphaCalibrator failed to converge: {exc}. "
                "Try providing a better alpha_init or more data points."
            ) from exc

        self.alpha_fitted = float(popt[0])
        self.residuals = G_measured - gate_model(delta_phi, self.alpha_fitted)
        return self

    def summary(self) -> str:
        """Return a human-readable summary of the calibration result."""
        lines = [
            "AlphaCalibrator summary",
            "─" * 40,
            f"  f0            : {self.f0} Hz" if self.f0 else "  f0            : (not set)",
            f"  α prior       : {self.alpha_prior:.4f}  (geometry-derived)",
        ]
        if self.alpha_fitted is not None:
            rmse = float(np.sqrt(np.mean(self.residuals**2)))
            lines += [
                f"  α fitted      : {self.alpha_fitted:.4f}",
                f"  RMSE          : {rmse:.6f}",
            ]
        else:
            lines.append("  α fitted      : (not yet fitted — call .fit())")
        return "\n".join(lines)

    def __repr__(self) -> str:
        fitted = f"{self.alpha_fitted:.4f}" if self.alpha_fitted is not None else "None"
        return (
            f"AlphaCalibrator(f0={self.f0}, alpha_prior={self.alpha_prior:.4f}, "
            f"alpha_fitted={fitted})"
        )
