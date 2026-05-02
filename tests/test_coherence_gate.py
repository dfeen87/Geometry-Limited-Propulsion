import numpy as np

from src.coherence_gate import coherence_gate, classical_delta_v


def test_coherence_gate_identity_at_zero():
    gate = coherence_gate(np.array([0.0, 0.0]), alpha=2.0)
    assert np.allclose(gate, 1.0)


def test_coherence_gate_rejects_nonpositive_alpha():
    try:
        coherence_gate([0.1], alpha=0.0)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for nonpositive alpha")


def test_classical_delta_v_positive():
    dv = classical_delta_v(isp=300.0, mass_initial=1000.0, mass_final=800.0)
    assert dv > 0
