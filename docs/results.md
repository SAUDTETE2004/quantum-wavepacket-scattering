# Results from the included experiments

Generated with the configurations stored in `figures/*_summary.json` and
the versions in `figures/environment.json`. All values below are measured,
except columns explicitly labeled as analytic predictions.

## Transmission and probability plateaus

| Quantity | Tunneling, E₀ = 0.7 | Above barrier, E₀ = 1.5 |
|---|---:|---:|
| Final time | 110 | 80 |
| Numerical T | 0.14485508 | 0.77880307 |
| Spectral T | 0.14506169 | 0.77987926 |
| Numerical R | 0.85514492 | 0.22119693 |
| Probability still in barrier | 2.7e-10 | 1.3e-12 |
| T change over final 10% of run | 2.3e-8 | 1.9e-9 |
| Maximum norm drift | 2.7e-13 | 4.1e-13 |
| Maximum outer-edge probability | 4.6e-10 | 2.0e-10 |

The plateau stability is much better than the transmission error. Continuing
the run will therefore not remove the discrepancy: resolution must improve.
The spectral average, not the carrier plane-wave value, is the reference.

## Gaussian fits and predictions

Widths below are conditional standard deviations within the outgoing packet.
Residuals compare density profiles using relative L2 norm.

| Case / channel | Fitted center | Filtered prediction | Fitted width | Predicted width | Shape residual |
|---|---:|---:|---:|---:|---:|
| Tunneling / reflected | -69.7403 | -69.7512 | 11.3868 | 11.3864 | 0.0236% |
| Tunneling / transmitted | 70.9379 | 70.9460 | 11.3691 | 11.3690 | 0.0108% |
| Above / reflected | -76.0615 | -76.0739 | 10.9721 | 10.9810 | 0.0920% |
| Above / transmitted | 77.5844 | 77.5608 | 10.8359 | 10.8376 | 0.0242% |

For tunneling, the phase-only transmitted center is 69.9207. Including the
real derivatives of log(t) moves it to 70.9460, near the fitted 70.9379.
The center shift cannot all be interpreted as phase delay: momentum filtering
changes the group velocity. The fitted transmitted velocity is 1.192699,
compared with the filtered prediction 1.192708 and incident velocity 1.183216.
The relative transmission phase delays are 0.19693 (tunneling) and 0.76991
(above barrier) in simulation time units.

The fits also report channel mass, peak density, and peak wavefunction
magnitude. For the tunneling transmitted packet these are 0.14485465,
0.00508297, and 0.07129494. Full values and late-time histories for both
channels are available as JSON and CSV. No statistical confidence interval
is claimed: these deterministic fits have discretization and model errors.

## Spatial convergence and independent solver

Absolute transmission errors for tunneling:

| dx | FFT, dt = 0.0025 | CN, dt = 0.01 |
|---|---:|---:|
| 0.20 | 7.419e-4 | 2.172e-3 |
| 0.10 | 1.805e-4 | 5.474e-4 |
| 0.05 | 4.749e-5 | 1.371e-4 |

Both methods approach the same continuum result. Roughly fourfold reductions
under halving dx are consistent with second-order spatial error over these
resolutions; this is an observed trend, not a proof of an asymptotic order.
At fixed dx = 0.1, T changes by 2.15e-5 from dt = 0.01 to 0.005 and by
4.61e-6 from 0.005 to 0.0025. These smaller steps demonstrate time convergence
toward a spatial-error floor.

Coarse time steps 0.04 and 0.02 produce nonmonotonic continuum errors and
outer-edge probabilities of 2.13e-4 and 4.29e-5. In particular, the apparently
smaller transmission error at 0.04 does **not** make it a better simulation.
Discontinuity-induced high-frequency errors and error cancellation require
checking several observables, not choosing the smallest scalar discrepancy.

## Boundary experiment

Half-box sizes 140, 180 and 220 agree in T within 1.5e-10. The deliberately
small half-box 100 has an outer-edge probability of 0.0388 and T = 0.147814,
visibly contaminated compared with the large-box plateau. It also truncates
the initial Gaussian more than the production boxes. Both problems make it
unsuitable. FFT wraparound is a boundary artifact; norm conservation cannot
detect it. The CN boundary is reflecting rather than periodic and is likewise
valid only before return.

## Sweeps and reproducibility checks

The height, width and energy sweeps contain dense spectral and plane-wave
curves plus nine numerical points, covering tunneling, threshold, and
above-barrier scattering. See `parameter_sweeps.json` for individual residuals
and boundary diagnostics. The nine points are illustrative, not each separately
resolution-extrapolated.

```sh
python -m pytest -q
python scripts/run_experiments.py --suite all --animate
python scripts/check_results.py
```

The final command checks generated spatial/time trends, large-box agreement,
the deliberately bad boundary diagnostic, norm conservation, transmission
plateaus, Gaussian residuals, centers, widths and velocities. Unit tests and
generated-result checks serve different purposes and should both be run.
