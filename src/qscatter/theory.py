"""Continuum amplitudes for a barrier centered at zero; hbar = m = 1."""
import numpy as np  # type: ignore
from scipy.integrate import simpson # type: ignore
from scipy.special import ndtr # type: ignore


def amplitudes(k, height, width):
    """Return t,r in the global exp(+-ikx) convention, for k > 0.

    A scaled hyperbolic expression avoids overflow for opaque barriers.
    sin(qL)/q is evaluated by its continuous limit at E = height.
    """
    k = np.asarray(k, dtype=float)
    if np.any(k <= 0) or height < 0 or width < 0:
        raise ValueError("Require positive k and nonnegative barrier height/width")
    q = np.sqrt(k*k - 2*height + 0j)
    z = q*width
    scale = np.abs(z.imag)
    # exp(+iz-scale) and exp(-iz-scale) never overflow.
    ep, em = np.exp(1j*z-scale), np.exp(-1j*z-scale)
    cosine = (ep+em)/2
    sine = (ep-em)/(2j)
    sinc = np.divide(sine, q, out=np.full_like(q, width), where=np.abs(q)>1e-12)
    denominator = cosine - 1j*(k*k+q*q)*sinc/(2*k)
    t = np.exp(-1j*k*width-scale)/denominator
    r = -1j*(k*k-q*q)*sinc/(2*k)*np.exp(-1j*k*width)/denominator
    return t, r


def spectrum(k, k0, sigma):
    """Unconditioned Gaussian momentum density; integral over R is one."""
    return np.sqrt(2/np.pi)*sigma*np.exp(-2*sigma*sigma*(k-k0)**2)


def packet_transmission(k0, sigma, height, width):
    """Positive-momentum conditional transmission; discarded weight reported."""
    k = np.linspace(max(1e-9, k0-9/(2*sigma)), k0+9/(2*sigma), 4097)
    w = spectrum(k, k0, sigma)
    t, _ = amplitudes(k, height, width)
    return float(simpson(w*abs(t)**2, x=k)/simpson(w, x=k))


def negative_weight(k0, sigma):
    return float(ndtr(-2*sigma*k0))


def gaussian_prediction(k0, sigma, x0, time, height, width, channel="T", filtering=True):
    """Second-order expansion of log(t) or log(r) near k0.

    Includes real log-amplitude derivatives when filtering=True. Invalid
    near amplitude zeros or when the local expansion is not integrable.
    """
    h = min(1e-3, k0/100)
    ks = k0 + h*np.arange(-3, 4)
    amp = amplitudes(ks, height, width)[0 if channel == "T" else 1]
    if np.min(abs(amp)) < 1e-12:
        raise ValueError("Phase expansion undefined near an amplitude zero")
    real = np.polynomial.polynomial.polyfit(ks-k0, np.log(abs(amp)), 4)
    imag = np.polynomial.polynomial.polyfit(ks-k0, np.unwrap(np.angle(amp)), 4)
    d1 = real[1] if filtering else 0.0
    d2 = 2*real[2] if filtering else 0.0
    phase1, phase2 = imag[1], 2*imag[2]
    a = sigma*sigma-d2/2
    b = (time-phase2)/2
    if a <= 0:
        raise ValueError("Local Gaussian approximation is not normalizable")
    direction = 1 if channel == "T" else -1
    center = direction*(x0+k0*time-phase1+b*d1/a)
    std = np.sqrt((a*a+b*b)/a)
    return {"center": float(center), "sigma": float(std),
            "velocity": float(direction*(k0+d1/(2*a))),
            "phase_prime": float(phase1), "phase_second": float(phase2),
            "relative_delay": float(phase1/k0)}

print(np.sqrt(2/np.pi)*0.5*np.exp(-2*0.5*0.5*(1-1)**2))