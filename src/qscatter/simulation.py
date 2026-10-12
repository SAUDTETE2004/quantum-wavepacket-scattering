"""Unitary Strang FFT propagation and independent finite-difference CN."""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.sparse import diags, eye
from scipy.sparse.linalg import splu


@dataclass(frozen=True)
class Config:
    height: float = 1.0
    width: float = 2.0
    energy: float = 0.7
    sigma: float = 10.0
    x0: float = -60.0
    half_box: float = 180.0
    dx: float = 0.1
    dt: float = 0.01
    final_time: float = 110.0
    snapshots: int = 161
    method: str = "fft"

    def validate(self):
        if min(self.width, self.energy, self.sigma, self.half_box, self.dx, self.dt, self.final_time) <= 0:
            raise ValueError("Lengths, energy and times must be positive")
        if self.height < 0 or self.snapshots < 2 or self.method not in ("fft", "cn"):
            raise ValueError("Invalid height, snapshot count or method")
        if self.width >= 2*self.half_box or not -self.half_box < self.x0 < -self.width/2:
            raise ValueError("Packet must start left of a barrier inside the box")
        for value in (2*self.half_box/self.dx, self.width/self.dx):
            if not np.isclose(value, round(value), rtol=0, atol=1e-8):
                raise ValueError("Box and barrier widths must be integer multiples of dx")
        if (round(2*self.half_box/self.dx)-round(self.width/self.dx)) % 2:
            raise ValueError("Barrier edges must align with cell faces")


def initial_packet(x, x0, sigma, k0):
    return (2*np.pi*sigma*sigma)**(-0.25)*np.exp(-(x-x0)**2/(4*sigma*sigma)+1j*k0*(x-x0))


def simulate(config=Config()):
    """Return sampled complex states and probabilities without renormalization.

    FFT boundaries are periodic; CN uses zero ghost values. Neither absorbs.
    Cell-centered quadrature partitions left/barrier/right without overlap.
    """
    c = config
    c.validate()
    n = round(2*c.half_box/c.dx)
    x = -c.half_box+(np.arange(n)+0.5)*c.dx
    v = c.height*(np.abs(x)<c.width/2)
    psi = initial_packet(x, c.x0, c.sigma, np.sqrt(2*c.energy))
    steps = int(np.ceil(c.final_time/c.dt))
    dt = c.final_time/steps
    if c.method == "fft":
        k = 2*np.pi*np.fft.fftfreq(n, c.dx)
        kick, drift = np.exp(-0.5j*v*dt), np.exp(-0.5j*k*k*dt)
        def advance(p):
            return kick*np.fft.ifft(drift*np.fft.fft(kick*p))
    else:
        h = diags([-np.ones(n-1)/(2*c.dx**2), 1/c.dx**2+v,
                   -np.ones(n-1)/(2*c.dx**2)], [-1, 0, 1], format="csc")
        lhs = splu(eye(n, format="csc")+0.5j*dt*h)
        rhs = eye(n, format="csc")-0.5j*dt*h
        def advance(p):
            return lhs.solve(rhs@p)
    save = set(np.linspace(0, steps, min(c.snapshots, steps+1), dtype=int))
    times, states, probabilities, edge, high_k = [], [], [], [], []
    masks = [x < -c.width/2, abs(x)<c.width/2, x>c.width/2]
    edge_mask = abs(x) > 0.9*c.half_box
    frequencies = np.fft.fftfreq(n)
    for step in range(steps+1):
        if step in save:
            density = abs(psi)**2
            times.append(step*dt)
            states.append(psi.copy())
            probabilities.append([float(c.dx*density[m].sum()) for m in masks])
            edge.append(float(c.dx*density[edge_mask].sum()))
            power = abs(np.fft.fft(psi))**2
            high_k.append(float(power[abs(frequencies)>0.4].sum()/power.sum()))
        if step < steps:
            psi = advance(psi)
    return {"config": asdict(c), "x": x, "time": np.array(times),
            "psi": np.array(states), "probabilities": np.array(probabilities),
            "edge_probability": np.array(edge), "high_k_fraction": np.array(high_k)}
