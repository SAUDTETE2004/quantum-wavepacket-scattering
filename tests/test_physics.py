from dataclasses import replace
import numpy as np
import pytest
from qscatter import Config, simulate
from qscatter.theory import amplitudes, packet_transmission


def test_unitarity_threshold_and_opaque_limit():
    k = np.r_[np.geomspace(.001, 10, 1000), np.sqrt(2)]
    t, r = amplitudes(k, 1, 2)
    np.testing.assert_allclose(abs(t)**2+abs(r)**2, 1, atol=2e-14)
    assert abs(abs(amplitudes(np.sqrt(2), 1, 2)[0])**2-1/3) < 1e-14
    t, r = amplitudes(.1, 100, 100)
    assert np.isfinite(t) and abs(abs(r)-1) < 1e-14


def test_free_and_resonance_limits():
    k = np.linspace(.01, 5, 101)
    t, r = amplitudes(k, 0, 2)
    np.testing.assert_allclose(t, 1, atol=1e-14)
    np.testing.assert_allclose(r, 0, atol=1e-14)
    assert abs(abs(amplitudes(np.sqrt(2+(np.pi/2)**2), 1, 2)[0])**2-1) < 1e-14


def test_free_gaussian_spreading():
    c = Config(height=0, final_time=20, snapshots=2)
    run = simulate(c)
    x, p = run["x"], abs(run["psi"][-1])**2
    mean = np.sum(x*p)*c.dx
    std = np.sqrt(np.sum((x-mean)**2*p)*c.dx)
    assert abs(mean-(c.x0+np.sqrt(2*c.energy)*c.final_time)) < 1e-9
    assert abs(std-np.sqrt(c.sigma**2+(c.final_time/(2*c.sigma))**2)) < 1e-9


@pytest.mark.parametrize("method", ["fft", "cn"])
def test_scattering_norm_and_spectral_transmission(method):
    c = Config(method=method, snapshots=11)
    run = simulate(c)
    assert np.max(abs(run["probabilities"].sum(axis=1)-1)) < 1e-10
    assert max(run["edge_probability"]) < 1e-5
    assert run["probabilities"][-1, 1] < 1e-6
    expected = packet_transmission(np.sqrt(2*c.energy), c.sigma, c.height, c.width)
    assert abs(run["probabilities"][-1, 2]-expected) < 2e-3


def test_grid_validation():
    with pytest.raises(ValueError):
        simulate(Config(width=2.03))


def test_complex_amplitudes_match_both_edges():
    # Propagate right-edge Cauchy data backward through the barrier.
    # This tests phases and the centered-coordinate convention, not just |t|.
    for k in [.5, 1.2, 2., 3.]:
        t, r = amplitudes(k, 1, 2)
        q = np.sqrt(k*k-2+0j)
        right = np.array([t*np.exp(1j*k), 1j*k*t*np.exp(1j*k)])
        transfer_back = np.array([[np.cos(2*q), -np.sin(2*q)/q],
                                  [q*np.sin(2*q), np.cos(2*q)]])
        left = np.array([np.exp(-1j*k)+r*np.exp(1j*k),
                         1j*k*(np.exp(-1j*k)-r*np.exp(1j*k))])
        np.testing.assert_allclose(transfer_back@right, left, atol=2e-14)



