# Numerical design and validation

## Primary method

In units $\hbar=m=1$, Strang splitting evolves one step by

$$e^{-iH\Delta t}\simeq e^{-iV\Delta t/2}
\mathcal F^{-1}e^{-ik^2\Delta t/2}\mathcal F e^{-iV\Delta t/2}.$$

Each factor is unitary. No step renormalizes the state. This makes norm drift
an informative implementation diagnostic, but norm conservation alone says
nothing about phase accuracy or boundary contamination. For a fixed spatial
discretization, Strang is second order in time; discontinuous potentials make
continuum error convergence less straightforward than for smooth potentials.
Do not infer accuracy from a generic time-step rule alone.

The grid is cell-centered, with both barrier edges on cell faces. An integer
number of cells represents the exact requested width; unsupported grids raise
an error rather than silently changing it. Regional probabilities are cell
sums times dx, partitioning the norm exactly. FFT kinetic propagation uses
the physical wave numbers from `fftfreq` with the factor $2\pi$.

## Independent check

Crank–Nicolson uses a three-point second derivative and sparse LU factorization:

$$(I+iH\Delta t/2)\psi_{n+1}=(I-iH\Delta t/2)\psi_n.$$

The finite-difference kinetic operator has its own dispersion error. The CN
box has zero ghost values, whereas FFT is periodic. Agreement before a packet
reaches the edges tests two genuinely different discretizations.

## Run the validation suite

```sh
python -m pytest -q
python scripts/run_experiments.py --suite validation
```

`figures/convergence.json` records every input and measured error. The studies:

| Study | Values | Held fixed / interpretation |
|---|---|---|
| Spatial refinement | dx = 0.2, 0.1, 0.05 | dt = 0.0025; compare with spectral T |
| Time refinement | dt = 0.04, 0.02, 0.01, 0.005 | dx = 0.1; spatial floor can remain |
| Box size | half-width = 100, 140, 180, 220 | Same dx, dt and observation time |
| Independent solver | CN dx = 0.2, 0.1, 0.05 | dt = 0.01; compare with continuum theory |

Diagnostics include norm drift, initial norm deficit, probability still inside
the barrier, plateau change, probability in the outer 10% of the domain, and
spectral weight in the highest 20% of represented wave numbers. Edge sampling
is a warning diagnostic, not a rigorous error bound; box enlargement and
agreement before return are essential. No absorber is used, so lost probability
cannot be mistaken for transmission. All runs are finite-time experiments.

## Reproduction and scope

`--suite all --animate` writes PNG plots, two GIFs, probability and fit-history
CSVs, JSON summaries and environment versions. The sweeps plot dense analytic
curves and nine TDSE sample points across height, width and carrier energy.
This is a modest reproducible sweep, not an exhaustive phase diagram.

The tests cover continuum unitarity, free and resonant limits, threshold and
opaque barriers, exact free-Gaussian motion/spreading, norm conservation,
two-solver transmission agreement, and grid rejection. The separate convergence
experiment establishes resolution trends; one transmission tolerance is not a
substitute for it. For new parameter sets, rerun convergence, move the initial
packet far enough from the barrier, and adjust final time and box size together.

Useful extensions are a broader-band distortion study, direct comparison of
the complex outgoing wavefunction with a spectral quadrature, probability
current and arrival-time measurements, or double-barrier resonances. Keep
phase-delay measurements separate from claims about tunneling traversal time.
