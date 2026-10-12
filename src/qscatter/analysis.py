"""Density fits report mass, center, standard deviation and shape residual."""
import numpy as np
from scipy.optimize import curve_fit


def gaussian(x, mass, center, sigma):
    return mass/(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5*((x-center)/sigma)**2)


def fit_packet(x, density, mask):
    x, y = x[mask], density[mask]
    dx = x[1]-x[0]
    mass = y.sum()*dx
    if mass < 1e-10:
        raise ValueError("Insufficient probability to fit")
    mean = (x*y).sum()*dx/mass
    std = np.sqrt(((x-mean)**2*y).sum()*dx/mass)
    p, _ = curve_fit(gaussian, x, y, p0=[mass, mean, std],
                     bounds=([0, x.min(), dx/2], [2, x.max(), np.ptp(x)]), maxfev=10000)
    fitted = gaussian(x, *p)
    return {"mass": float(p[0]), "center": float(p[1]), "sigma": float(p[2]),
            "density_peak": float(p[0]/(np.sqrt(2*np.pi)*p[2])),
            "wavefunction_peak_magnitude": float(np.sqrt(p[0]/(np.sqrt(2*np.pi)*p[2]))),
            "relative_L2_residual": float(np.linalg.norm(y-fitted)/np.linalg.norm(y))}

