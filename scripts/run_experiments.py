"""Generate portfolio figures, machine-readable results and optional animation.

Run from the repository root after installing with pip install -e '.[dev]'.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import platform
import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from qscatter import Config, simulate
from qscatter.theory import packet_transmission, amplitudes, gaussian_prediction, negative_weight
from qscatter.analysis import fit_packet, gaussian

plt.rcParams.update({"figure.dpi": 130, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.2})


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2)+"\n", encoding="utf-8")


def summarize(run):
    c = Config(**run["config"])
    p = run["probabilities"]
    theory = packet_transmission(np.sqrt(2*c.energy), c.sigma, c.height, c.width)
    return {"config": run["config"], "T_numerical": float(p[-1, 2]),
            "R_numerical": float(p[-1, 0]), "barrier_remaining": float(p[-1, 1]),
            "T_theory": theory, "absolute_T_error": float(abs(p[-1, 2]-theory)),
            "norm_drift": float(np.max(abs(p.sum(axis=1)-p[0].sum()))),
            "initial_norm_error": float(abs(p[0].sum()-1)),
            "max_edge_probability": float(max(run["edge_probability"])),
            "max_high_k_fraction": float(max(run["high_k_fraction"])),
            "plateau_change_last_10_percent": float(abs(p[-1, 2]-p[int(.9*len(p)), 2])),
            "initial_negative_momentum_weight": negative_weight(np.sqrt(2*c.energy), c.sigma)}


def single(c, name, out, animate=False):
    run = simulate(c)
    x, times, density = run["x"], run["time"], abs(run["psi"])**2
    summary = summarize(run)
    print(name, summary["T_numerical"], summary["T_theory"], flush=True)
    fig, axs = plt.subplots(2, 2, figsize=(11, 6), sharex=True, sharey=True)
    collision = -c.x0/np.sqrt(2*c.energy)
    for ax, target, label in zip(axs.flat, [0, collision, collision+18, c.final_time],
                                  ["Incident", "Barrier interaction", "Separation", "Outgoing packets"]):
        i = np.argmin(abs(times-target))
        ax.plot(x, density[i], color="#126782")
        ax.axvspan(-c.width/2, c.width/2, color="#f4a261", alpha=.6)
        ax.set_title(f"{label}: t = {times[i]:.1f}")
        ax.set(xlabel="x", ylabel="Probability density", xlim=(-125, 125))
    fig.suptitle(f"{name.replace('_',' ').title()} | E₀/U₀ = {c.energy/c.height:.2f}")
    fig.tight_layout()
    fig.savefig(out/f"{name}_snapshots.png")
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for j, label in enumerate(["Left region (late-time R)", "Inside barrier", "Right region (T)"]):
        axs[0].plot(times, run["probabilities"][:, j], label=label)
    axs[0].axhline(summary["T_theory"], color="black", ls="--", label="Spectral T")
    axs[0].axhline(1-summary["T_theory"], color="gray", ls=":", label="Spectral R")
    axs[0].set(xlabel="Time", ylabel="Probability", ylim=(-.03, 1.03))
    axs[0].legend(fontsize=8)
    axs[1].semilogy(times, np.maximum(abs(run["probabilities"].sum(axis=1)-1), 1e-17), label="Norm error")
    axs[1].semilogy(times, np.maximum(run["edge_probability"], 1e-17), label="Outer 10% probability")
    axs[1].set(xlabel="Time", ylabel="Diagnostic magnitude")
    axs[1].legend()
    fig.tight_layout(); fig.savefig(out/f"{name}_probabilities.png"); plt.close(fig)
    fig, axs = plt.subplots(2, 2, figsize=(11, 7))
    summary["fits"] = {}
    for col, channel in enumerate(["R", "T"]):
        mask = x < -c.width/2 if channel == "R" else x > c.width/2
        fit = fit_packet(x, density[-1], mask)
        simple = gaussian_prediction(np.sqrt(2*c.energy), c.sigma, c.x0, times[-1], c.height, c.width, channel, False)
        filtered = gaussian_prediction(np.sqrt(2*c.energy), c.sigma, c.x0, times[-1], c.height, c.width, channel, True)
        mass = summary["T_theory"] if channel == "T" else 1-summary["T_theory"]
        filtered["density_peak"] = mass/(np.sqrt(2*np.pi)*filtered["sigma"])
        summary["fits"][channel] = {"fit": fit, "phase_only_prediction": simple, "filtered_prediction": filtered}
        ax = axs[0, col]
        ax.plot(x[mask], density[-1, mask], label="TDSE", lw=2)
        ax.plot(x[mask], gaussian(x[mask], fit["mass"], fit["center"], fit["sigma"]), "--", label="Gaussian fit")
        ax.plot(x[mask], gaussian(x[mask], mass, filtered["center"], filtered["sigma"]), ":", label="Analytic, filtered")
        ax.set(title=f"{channel}: relative L2 residual {fit['relative_L2_residual']:.3g}", xlabel="x", ylabel="Density")
        ax.legend(fontsize=8)
        late = np.where(times >= .85*c.final_time)[0]
        records = []
        for i in late:
            f = fit_packet(x, density[i], mask)
            pred = gaussian_prediction(np.sqrt(2*c.energy), c.sigma, c.x0, times[i], c.height, c.width, channel)
            records.append([times[i], f["center"], f["sigma"], pred["center"], pred["sigma"]])
        a = np.asarray(records)
        velocity = np.polyfit(a[:, 0], a[:, 1], 1)[0]
        summary["fits"][channel]["fitted_velocity"] = float(velocity)
        np.savetxt(out/f"{name}_{channel}_fit_history.csv", a, delimiter=",", header="time,fit_center,fit_sigma,pred_center,pred_sigma", comments="")
        axs[1, col].plot(a[:, 0], a[:, 2], label="Fitted width")
        axs[1, col].plot(a[:, 0], a[:, 4], "--", label="Predicted width")
        axs[1, col].set(xlabel="Time", ylabel="Conditional standard deviation")
        axs[1, col].legend()
    fig.tight_layout(); fig.savefig(out/f"{name}_gaussian_fits.png"); plt.close(fig)
    np.savetxt(out/f"{name}_probabilities.csv", np.column_stack([times, run["probabilities"], run["edge_probability"]]), delimiter=",", header="time,P_left,P_barrier,P_right,P_edge", comments="")
    write_json(out/f"{name}_summary.json", summary)
    if animate:
        fig, ax = plt.subplots(figsize=(9, 4))
        line, = ax.plot(x, density[0], color="#126782")
        ax.axvspan(-c.width/2, c.width/2, color="#f4a261", alpha=.6)
        ax.set(xlim=(-125, 125), ylim=(0, density.max()*1.08), xlabel="x", ylabel="Probability density")
        title = ax.set_title("")
        def update(i):
            line.set_ydata(density[i]); title.set_text(f"{name}: t = {times[i]:.1f}")
            return line, title
        animation = FuncAnimation(fig, update, frames=range(0, len(times), 2), blit=True)
        animation.save(out/f"{name}.gif", writer=PillowWriter(fps=15))
        plt.close(fig)
    return summary


def sweeps(out):
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    all_rows = []
    base = Config(dx=.1, dt=.01, snapshots=31)
    for ax, field, values in zip(axs, ["height", "width", "energy"],
                                 [np.linspace(.3, 1.8, 61), np.linspace(.4, 4, 61), np.linspace(.3, 2, 61)]):
        predictions = [packet_transmission(np.sqrt(2*(v if field=="energy" else base.energy)), base.sigma,
                       v if field=="height" else base.height, v if field=="width" else base.width) for v in values]
        ax.plot(values, predictions, label="Spectral average")
        mono = [abs(amplitudes(np.sqrt(2*(v if field=="energy" else base.energy)),
                 v if field=="height" else base.height, v if field=="width" else base.width)[0])**2 for v in values]
        ax.plot(values, mono, "--", alpha=.65, label="Plane wave")
        for v in ([.4, 1., 1.6] if field=="height" else [.6, 2., 3.6] if field=="width" else [.4, 1., 1.7]):
            c = replace(base, **{field:v})
            c = replace(c, final_time=130/np.sqrt(2*c.energy))
            result = summarize(simulate(c))
            all_rows.append({"sweep":field, "value":v, **result})
            ax.plot(v, result["T_numerical"], "o", color="#e76f51")
        ax.set(xlabel=field, ylabel="Transmission", ylim=(0, 1.03))
        ax.legend(fontsize=8)
    fig.suptitle("Barrier sweeps: lines = theory, points = independent TDSE runs")
    fig.tight_layout(); fig.savefig(out/"parameter_sweeps.png"); plt.close(fig)
    write_json(out/"parameter_sweeps.json", all_rows)


def convergence(out):
    base = Config(snapshots=31)
    cases = [("dx", replace(base, dx=d, dt=.0025)) for d in [.2, .1, .05]]
    cases += [("dt", replace(base, dt=t)) for t in [.04, .02, .01, .005]]
    cases += [("box", replace(base, half_box=b)) for b in [100., 140., 180., 220.]]
    cases += [("cn", replace(base, method="cn", dx=d)) for d in [.2, .1, .05]]
    rows = []
    for label, c in cases:
        result = summarize(simulate(c))
        rows.append({"study":label, **result})
        print("validation", label, c.dx, c.dt, c.half_box, result["absolute_T_error"], flush=True)
    write_json(out/"convergence.json", rows)
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    for ax, study, field in zip(axs, ["dx", "dt", "box"], ["dx", "dt", "half_box"]):
        subset = [r for r in rows if r["study"]==study]
        ax.loglog([r["config"][field] for r in subset], [max(r["absolute_T_error"], 1e-16) for r in subset], "o-")
        ax.set(xlabel=field, ylabel="Absolute transmission error")
    fig.tight_layout(); fig.savefig(out/"convergence.png"); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=["demo", "sweeps", "validation", "all"], default="demo")
    parser.add_argument("--animate", action="store_true")
    parser.add_argument("--out", type=Path, default=Path("figures"))
    args = parser.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    write_json(args.out/"environment.json", {"python":platform.python_version(), "numpy":np.__version__, "scipy":scipy.__version__, "matplotlib":matplotlib.__version__})
    if args.suite in ("demo", "all"):
        single(Config(), "tunneling", args.out, args.animate)
        single(Config(energy=1.5, final_time=80), "above_barrier", args.out, args.animate)
    if args.suite in ("sweeps", "all"): sweeps(args.out)
    if args.suite in ("validation", "all"): convergence(args.out)


if __name__ == "__main__":
    main()
