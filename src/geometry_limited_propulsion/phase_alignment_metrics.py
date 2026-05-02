"""
phase_alignment_metrics.py
==========================
Measurable proxies for the phase/geometric deviation Δφ(t) introduced in
the HLV–AILEE formulation (Section 7.1 of the paper).

Three proxy families are provided, each addressing a different observable
channel:

    Temporal   Δφ_temp  ∼  |f_sys − f_lat| / f_lat             (Eq. 7)
    Spatial    Δφ_spat  ∼  ⟨|∇ × v_flow|⟩ / ⟨|v_flow|⟩       (Eq. 8)
    Spectral   Δφ_spec  ∼  ∫ dω |S_thrust(ω) − S_lat(ω)|       (Eq. 9)

A composite Δφ combiner is also provided for use with ``coherence_gate.py``.

All functions accept NumPy arrays and return dimensionless scalars or arrays
suitable for direct input to ``coherence_gate()``.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


# ---------------------------------------------------------------------------
# Temporal proxy  (Eq. 7)
# ---------------------------------------------------------------------------

def delta_phi_temporal(
    f_sys: ArrayLike,
    f_lat: float | ArrayLike,
) -> np.ndarray:
    """Temporal phase-deviation proxy Δφ_temp (Eq. 7).

    Measures fractional mismatch between the system's characteristic
    frequency and the host-lattice reference frequency:

        Δφ_temp(t)  =  |f_sys(t) − f_lat| / f_lat

    In a propulsion context *f_sys* might be the dominant frequency of the
    thrust time-series, the fuel-injection cycle, or a combustion oscillation
    mode.  *f_lat* is the theoretically admissible lattice frequency for the
    current operating condition.

    Parameters
    ----------
    f_sys : array-like, shape (N,) or scalar
        System frequency (or frequencies) [Hz or normalised].
    f_lat : float or array-like, shape (N,)
        Lattice reference frequency [same units as *f_sys*].
        Must be strictly positive everywhere.

    Returns
    -------
    np.ndarray
        Dimensionless Δφ_temp ≥ 0, same shape as *f_sys*.

    Examples
    --------
    >>> delta_phi_temporal(f_sys=100.0, f_lat=95.0)
    array(0.05263158)
    >>> import numpy as np
    >>> delta_phi_temporal(np.array([95, 100, 110]), f_lat=100.0)
    array([0.05, 0.  , 0.1 ])
    """
    f_sys = np.asarray(f_sys, dtype=float)
    f_lat = np.asarray(f_lat, dtype=float)
    if np.any(f_lat <= 0):
        raise ValueError("f_lat must be strictly positive.")
    return np.abs(f_sys - f_lat) / f_lat


# ---------------------------------------------------------------------------
# Spatial / vorticity proxy  (Eq. 8)
# ---------------------------------------------------------------------------

def delta_phi_spatial(
    vflow: np.ndarray,
    dx: float = 1.0,
    dy: float = 1.0,
) -> float:
    """Spatial phase-deviation proxy Δφ_spat (Eq. 8).

    Approximates the ratio of mean vorticity magnitude to mean flow speed:

        Δφ_spat  =  ⟨|∇ × v_flow|⟩ / ⟨|v_flow|⟩

    The 2-D curl (z-component of ∇ × v) is estimated via finite differences.

    Parameters
    ----------
    vflow : np.ndarray, shape (2, M, N)
        Velocity field array where

            vflow[0]  →  v_x  (x-component), shape (M, N)
            vflow[1]  →  v_y  (y-component), shape (M, N)

        The spatial domain is assumed uniform with grid spacings *dx* and
        *dy*.
    dx : float
        Grid spacing in the x-direction [same units as velocity · time].
    dy : float
        Grid spacing in the y-direction [same units as velocity · time].

    Returns
    -------
    float
        Dimensionless Δφ_spat ≥ 0.

    Raises
    ------
    ValueError
        If *vflow* does not have shape (2, M, N) or if ⟨|v_flow|⟩ ≈ 0.

    Notes
    -----
    For 3-D flow fields the full curl |∇ × v| = √(ωx² + ωy² + ωz²) should
    be used; this 2-D helper is provided as a starting point.

    Examples
    --------
    >>> import numpy as np
    >>> # Uniform flow — zero vorticity → Δφ_spat = 0
    >>> vx = np.ones((4, 4))
    >>> vy = np.zeros((4, 4))
    >>> delta_phi_spatial(np.stack([vx, vy]))
    0.0
    """
    if vflow.ndim != 3 or vflow.shape[0] != 2:
        raise ValueError(
            "vflow must have shape (2, M, N) — [v_x, v_y] over a 2-D grid."
        )
    if dx <= 0 or dy <= 0:
        raise ValueError(f"dx and dy must be strictly positive, got dx={dx}, dy={dy}.")
    vx, vy = vflow[0], vflow[1]

    # z-component of curl: ∂vy/∂x − ∂vx/∂y
    dvy_dx = np.gradient(vy, dx, axis=1)
    dvx_dy = np.gradient(vx, dy, axis=0)
    curl_z = dvy_dx - dvx_dy

    mean_vorticity = float(np.mean(np.abs(curl_z)))
    mean_speed = float(np.mean(np.sqrt(vx**2 + vy**2)))

    if mean_speed < 1e-12:
        raise ValueError(
            "Mean flow speed is effectively zero; Δφ_spat is undefined."
        )
    return mean_vorticity / mean_speed


# ---------------------------------------------------------------------------
# Spectral proxy  (Eq. 9)
# ---------------------------------------------------------------------------

def delta_phi_spectral(
    thrust_signal: ArrayLike,
    lat_signal: ArrayLike,
    sample_rate: float,
    *,
    window: str = "hann",
    nperseg: int | None = None,
) -> float:
    """Spectral phase-deviation proxy Δφ_spec (Eq. 9).

    Computes the L¹ distance between the power spectral density of the
    thrust time-series and a reference (lattice) spectrum:

        Δφ_spec  =  ∫ dω  |S_thrust(ω) − S_lat(ω)|

    Both spectra are normalised to unit area before differencing so that
    the result is a pure shape-mismatch metric, independent of absolute
    power levels.

    Parameters
    ----------
    thrust_signal : array-like, shape (N,)
        Measured thrust time-series [any consistent units].
    lat_signal : array-like, shape (N,)
        Reference / lattice signal [same units as *thrust_signal*].
        This might be a synthesised signal derived from the theoretically
        admissible lattice frequency, or a reference run baseline.
    sample_rate : float
        Sampling rate [Hz].
    window : str
        SciPy-compatible window name passed to ``scipy.signal.welch``.
        Default ``"hann"``.
    nperseg : int or None
        Segment length for Welch PSD estimation.  ``None`` → SciPy default
        (signal length // 8, floored at 256).

    Returns
    -------
    float
        Dimensionless Δφ_spec ≥ 0.  Values near 0 indicate well-matched
        spectra; larger values indicate increasing spectral mismatch.

    Raises
    ------
    ImportError
        If SciPy is not installed.
    ValueError
        If the two signals have different lengths, or if either PSD
        integrates to (near) zero.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> t = np.linspace(0, 1, 1024)
    >>> sig_a = np.sin(2 * np.pi * 50 * t)   # 50 Hz tone
    >>> sig_b = np.sin(2 * np.pi * 80 * t)   # 80 Hz tone — mismatch
    >>> delta_phi_spectral(sig_a, sig_b, sample_rate=1024.0)  # doctest: +SKIP
    1.9...
    """
    try:
        from scipy.signal import welch
    except ImportError as exc:
        raise ImportError(
            "SciPy is required for delta_phi_spectral. "
            "Install it with: pip install scipy"
        ) from exc

    thrust_signal = np.asarray(thrust_signal, dtype=float)
    lat_signal = np.asarray(lat_signal, dtype=float)

    if thrust_signal.shape != lat_signal.shape:
        raise ValueError(
            "thrust_signal and lat_signal must have the same length; "
            f"got {thrust_signal.shape} vs {lat_signal.shape}."
        )

    kw = dict(fs=sample_rate, window=window, nperseg=nperseg)
    freqs, S_thrust = welch(thrust_signal, **kw)
    _,     S_lat    = welch(lat_signal,    **kw)

    # Normalise to unit area (probability-density-like comparison)
    dfreq = freqs[1] - freqs[0]
    norm_thrust = np.trapz(S_thrust, dx=dfreq)
    norm_lat    = np.trapz(S_lat,    dx=dfreq)

    if norm_thrust < 1e-30:
        raise ValueError("Thrust PSD integrates to ~0; check the input signal.")
    if norm_lat < 1e-30:
        raise ValueError("Lattice PSD integrates to ~0; check the reference signal.")

    S_thrust_n = S_thrust / norm_thrust
    S_lat_n    = S_lat    / norm_lat

    return float(np.trapz(np.abs(S_thrust_n - S_lat_n), dx=dfreq))


# ---------------------------------------------------------------------------
# Composite combiner
# ---------------------------------------------------------------------------

def composite_delta_phi(
    *proxies: float | np.ndarray,
    weights: ArrayLike | None = None,
) -> np.ndarray:
    """Combine multiple Δφ proxy values into a single composite metric.

    A simple weighted root-mean-square combination is used:

        Δφ_composite = √( Σᵢ wᵢ · Δφᵢ² )   with  Σᵢ wᵢ = 1

    Parameters
    ----------
    *proxies : float or np.ndarray
        Any number of Δφ proxy values.  All array arguments must be
        broadcastable to a common shape.
    weights : array-like of length len(proxies), optional
        Non-negative weights for each proxy.  Will be normalised to sum to 1.
        Default: equal weighting.

    Returns
    -------
    np.ndarray
        Composite Δφ ≥ 0.

    Examples
    --------
    >>> composite_delta_phi(0.1, 0.2, 0.3)           # equal weights
    array(0.21602469)
    >>> composite_delta_phi(0.1, 0.2, 0.3, weights=[0.5, 0.25, 0.25])
    array(0.15811388)
    """
    if len(proxies) == 0:
        raise ValueError("At least one proxy must be supplied.")

    proxy_arrays = [np.asarray(p, dtype=float) for p in proxies]

    if weights is None:
        w = np.ones(len(proxy_arrays)) / len(proxy_arrays)
    else:
        w = np.asarray(weights, dtype=float)
        if w.shape != (len(proxy_arrays),):
            raise ValueError(
                f"weights must have length {len(proxy_arrays)}, got {w.shape}."
            )
        if np.any(w < 0):
            raise ValueError("All weights must be non-negative.")
        total = w.sum()
        if total < 1e-12:
            raise ValueError("Weights must not all be zero.")
        w = w / total

    return np.sqrt(sum(wi * phi_i**2 for wi, phi_i in zip(w, proxy_arrays)))
