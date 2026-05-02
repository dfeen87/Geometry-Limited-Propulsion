import numpy as np
import pytest

from geometry_limited_propulsion.coherence_gate import (
    coherence_gate,
    classical_delta_v,
    hlv_ailee_delta_v,
    gate_vs_phase_sweep,
    hysteresis_sweep,
)
from geometry_limited_propulsion.phase_alignment_metrics import (
    delta_phi_temporal,
    delta_phi_spatial,
    composite_delta_phi,
)
from geometry_limited_propulsion.lattice_reference import (
    phi_harmonic_series,
    nearest_lattice_frequency,
    psi_lat,
    acceptance_window_width,
    alpha_from_window_width,
    AlphaCalibrator,
    PHI,
)


# ---------------------------------------------------------------------------
# coherence_gate
# ---------------------------------------------------------------------------

def test_coherence_gate_identity_at_zero():
    gate = coherence_gate(np.array([0.0, 0.0]), alpha=2.0)
    assert np.allclose(gate, 1.0)


def test_coherence_gate_decays_with_deviation():
    phi = np.array([0.0, 0.5, 1.0])
    gate = coherence_gate(phi, alpha=2.0)
    assert gate[0] > gate[1] > gate[2]
    assert np.allclose(gate, np.exp(-2.0 * phi**2))


def test_coherence_gate_rejects_zero_alpha():
    with pytest.raises(ValueError):
        coherence_gate([0.1], alpha=0.0)


def test_coherence_gate_rejects_negative_alpha():
    with pytest.raises(ValueError):
        coherence_gate([0.1], alpha=-1.0)


# ---------------------------------------------------------------------------
# classical_delta_v
# ---------------------------------------------------------------------------

def test_classical_delta_v_positive():
    dv = classical_delta_v(isp=300.0, mass_initial=1000.0, mass_final=800.0)
    assert dv > 0


def test_classical_delta_v_tsiolkovsky():
    # Δv = Isp * g0 * ln(M0/Mf)
    from geometry_limited_propulsion.coherence_gate import G0
    dv = classical_delta_v(isp=300.0, mass_initial=1000.0, mass_final=500.0)
    expected = 300.0 * G0 * np.log(1000.0 / 500.0)
    assert np.isclose(dv, expected)


def test_classical_delta_v_rejects_nonpositive_isp():
    with pytest.raises(ValueError):
        classical_delta_v(isp=0.0, mass_initial=1000.0, mass_final=800.0)


def test_classical_delta_v_rejects_bad_masses():
    with pytest.raises(ValueError):
        classical_delta_v(isp=300.0, mass_initial=800.0, mass_final=1000.0)
    with pytest.raises(ValueError):
        classical_delta_v(isp=300.0, mass_initial=1000.0, mass_final=0.0)


# ---------------------------------------------------------------------------
# hlv_ailee_delta_v
# ---------------------------------------------------------------------------

def _simple_arrays(n=10):
    t = np.linspace(0, 1, n)
    p_in = np.ones(n) * 1000.0
    mass = np.ones(n) * 500.0
    velocity = np.ones(n) * 10.0
    delta_phi = np.zeros(n)
    return t, p_in, mass, velocity, delta_phi


def test_hlv_ailee_delta_v_returns_float():
    t, p_in, mass, velocity, delta_phi = _simple_arrays()
    result = hlv_ailee_delta_v(
        t, mass, delta_phi,
        isp=300.0, eta=1.0, alpha=2.0,
        mode='power', p_in=p_in, velocity=velocity,
    )
    assert isinstance(result, float)
    assert result > 0


def test_hlv_ailee_delta_v_rejects_nonpositive_isp():
    t, p_in, mass, velocity, delta_phi = _simple_arrays()
    with pytest.raises(ValueError):
        hlv_ailee_delta_v(
            t, mass, delta_phi,
            isp=0.0, eta=1.0, alpha=2.0,
            mode='power', p_in=p_in, velocity=velocity,
        )


def test_hlv_ailee_delta_v_rejects_bad_eta():
    t, p_in, mass, velocity, delta_phi = _simple_arrays()
    with pytest.raises(ValueError):
        hlv_ailee_delta_v(
            t, mass, delta_phi,
            isp=300.0, eta=0.0, alpha=2.0,
            mode='power', p_in=p_in, velocity=velocity,
        )
    with pytest.raises(ValueError):
        hlv_ailee_delta_v(
            t, mass, delta_phi,
            isp=300.0, eta=1.5, alpha=2.0,
            mode='power', p_in=p_in, velocity=velocity,
        )


def test_hlv_ailee_delta_v_rejects_shape_mismatch():
    t, p_in, mass, velocity, delta_phi = _simple_arrays()
    with pytest.raises(ValueError):
        hlv_ailee_delta_v(
            t, mass[:-1], delta_phi,
            isp=300.0, eta=1.0, alpha=2.0,
            mode='power', p_in=p_in, velocity=velocity,
        )


def test_hlv_ailee_delta_v_zero_gate_gives_zero():
    t, p_in, mass, velocity, _ = _simple_arrays()
    # Very large delta_phi collapses gate to ~0
    delta_phi = np.ones(len(t)) * 1000.0
    result = hlv_ailee_delta_v(
        t, mass, delta_phi,
        isp=300.0, eta=1.0, alpha=2.0,
        mode='power', p_in=p_in, velocity=velocity,
    )
    assert result == pytest.approx(0.0, abs=1e-10)


# ---------------------------------------------------------------------------
# gate_vs_phase_sweep / hysteresis_sweep
# ---------------------------------------------------------------------------

def test_gate_vs_phase_sweep_shape():
    phi = np.linspace(0, 2, 50)
    phi_out, gate_out = gate_vs_phase_sweep(phi, alpha=1.0)
    assert phi_out.shape == gate_out.shape == phi.shape
    assert np.allclose(gate_out, np.exp(-1.0 * phi**2))


def test_hysteresis_sweep_keys_and_shapes():
    phi_up = np.linspace(0, 1, 20)
    phi_down = np.linspace(1, 0, 20)
    result = hysteresis_sweep(phi_up, phi_down, alpha_up=3.0, alpha_down=1.5)
    assert set(result.keys()) == {"phi_up", "gate_up", "phi_down", "gate_down"}
    assert result["gate_up"].shape == phi_up.shape
    assert result["gate_down"].shape == phi_down.shape
    # Up-sweep should collapse faster (larger alpha)
    assert result["gate_up"][10] < result["gate_down"][10]


# ---------------------------------------------------------------------------
# delta_phi_temporal
# ---------------------------------------------------------------------------

def test_delta_phi_temporal_zero_at_exact_match():
    result = delta_phi_temporal(100.0, f_lat=100.0)
    assert float(result) == 0.0


def test_delta_phi_temporal_fractional():
    result = delta_phi_temporal(105.0, f_lat=100.0)
    assert np.isclose(float(result), 0.05)


def test_delta_phi_temporal_rejects_zero_f_lat():
    with pytest.raises(ValueError):
        delta_phi_temporal(100.0, f_lat=0.0)


# ---------------------------------------------------------------------------
# delta_phi_spatial
# ---------------------------------------------------------------------------

def test_delta_phi_spatial_uniform_flow():
    vx = np.ones((4, 4))
    vy = np.zeros((4, 4))
    assert delta_phi_spatial(np.stack([vx, vy])) == 0.0


def test_delta_phi_spatial_rejects_wrong_shape():
    with pytest.raises(ValueError):
        delta_phi_spatial(np.ones((3, 4, 4)))


# ---------------------------------------------------------------------------
# composite_delta_phi
# ---------------------------------------------------------------------------

def test_composite_delta_phi_equal_weights():
    result = composite_delta_phi(0.0, 0.0)
    assert float(result) == 0.0


def test_composite_delta_phi_value():
    result = composite_delta_phi(3.0, 4.0)
    # sqrt((3^2 + 4^2) / 2) = sqrt(12.5)
    assert np.isclose(float(result), np.sqrt(12.5))


def test_composite_delta_phi_custom_weights():
    result = composite_delta_phi(1.0, 0.0, weights=[1.0, 0.0])
    assert np.isclose(float(result), 1.0)


def test_composite_delta_phi_rejects_empty():
    with pytest.raises(ValueError):
        composite_delta_phi()


# ---------------------------------------------------------------------------
# phi_harmonic_series
# ---------------------------------------------------------------------------

def test_phi_harmonic_series_contains_f0():
    series = phi_harmonic_series(100.0, n_range=(0, 0), m_range=(1, 1))
    assert np.isclose(series[0], 100.0)


def test_phi_harmonic_series_sorted_ascending():
    series = phi_harmonic_series(100.0)
    assert np.all(np.diff(series) > 0)


def test_phi_harmonic_series_rejects_nonpositive_f0():
    with pytest.raises(ValueError):
        phi_harmonic_series(0.0)


# ---------------------------------------------------------------------------
# nearest_lattice_frequency
# ---------------------------------------------------------------------------

def test_nearest_lattice_frequency_exact_node():
    f_lat, dphi = nearest_lattice_frequency(100.0, f0=100.0, n_range=(0, 0), m_range=(1, 1))
    assert np.isclose(f_lat[0], 100.0)
    assert np.isclose(dphi[0], 0.0)


def test_nearest_lattice_frequency_shape():
    f_sys = np.array([90.0, 100.0, 110.0])
    f_lat, dphi = nearest_lattice_frequency(f_sys, f0=100.0)
    assert f_lat.shape == f_sys.shape
    assert dphi.shape == f_sys.shape


# ---------------------------------------------------------------------------
# psi_lat
# ---------------------------------------------------------------------------

def test_psi_lat_signal_shape():
    t = np.linspace(0, 1, 512)
    ref = psi_lat(t, f0=50.0)
    assert ref["signal"].shape == t.shape


def test_psi_lat_frequencies():
    t = np.linspace(0, 1, 512)
    ref = psi_lat(t, f0=100.0)
    assert np.isclose(ref["f_phi"], 100.0 * PHI)
    assert np.isclose(ref["f_sub"], 100.0 / PHI)


def test_psi_lat_rejects_large_epsilon():
    t = np.linspace(0, 1, 64)
    with pytest.raises(ValueError):
        psi_lat(t, f0=100.0, epsilon=0.5)


def test_psi_lat_warns_medium_epsilon():
    t = np.linspace(0, 1, 64)
    with pytest.warns(UserWarning):
        psi_lat(t, f0=100.0, epsilon=0.2)


def test_psi_lat_negative_epsilon_warns_correctly():
    # Bug fix: abs() check — epsilon=-0.2 should also trigger the warning
    t = np.linspace(0, 1, 64)
    with pytest.warns(UserWarning):
        psi_lat(t, f0=100.0, epsilon=-0.2)


# ---------------------------------------------------------------------------
# acceptance_window_width / alpha_from_window_width
# ---------------------------------------------------------------------------

def test_acceptance_window_width_default():
    hw = acceptance_window_width(n_range=(-4, 4))
    assert np.isclose(hw, 1.0 / PHI**4)


def test_alpha_from_window_width_roundtrip():
    hw = 0.2
    alpha = alpha_from_window_width(hw)
    assert np.isclose(alpha, 1.0 / hw**2)


def test_alpha_from_window_width_rejects_nonpositive():
    with pytest.raises(ValueError):
        alpha_from_window_width(0.0)


# ---------------------------------------------------------------------------
# AlphaCalibrator
# ---------------------------------------------------------------------------

def test_alpha_calibrator_prior_without_f0():
    cal = AlphaCalibrator()
    assert cal.alpha_prior == 2.0
    assert cal.alpha_fitted is None


def test_alpha_calibrator_prior_with_f0():
    cal = AlphaCalibrator(f0=100.0)
    expected_hw = acceptance_window_width(n_range=(-4, 4))
    expected_alpha = alpha_from_window_width(expected_hw)
    assert np.isclose(cal.alpha_prior, expected_alpha)


def test_alpha_calibrator_fit_recovers_true_alpha():
    pytest.importorskip("scipy")
    rng = np.random.default_rng(42)
    true_alpha = 3.0
    phi_vals = np.linspace(0, 1.5, 50)
    G_obs = np.exp(-true_alpha * phi_vals**2) + rng.normal(0, 0.005, 50)
    cal = AlphaCalibrator(f0=100.0)
    cal.fit(phi_vals, G_obs)
    assert np.isclose(cal.alpha_fitted, true_alpha, atol=0.1)
    assert cal.residuals is not None


def test_alpha_calibrator_summary_unfitted():
    cal = AlphaCalibrator()
    summary = cal.summary()
    assert "not yet fitted" in summary


def test_alpha_calibrator_repr():
    cal = AlphaCalibrator(f0=50.0)
    r = repr(cal)
    assert "AlphaCalibrator" in r
    assert "alpha_fitted=None" in r
