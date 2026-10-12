"""Validate an existing full experiment suite; run after --suite all.

These checks test generated experiments, not just analytic identities.
"""
import json
from pathlib import Path
import argparse


def check(directory):
    rows = json.loads((directory/"convergence.json").read_text())
    spatial = [r for r in rows if r["study"] == "dx"]
    cn = [r for r in rows if r["study"] == "cn"]
    for study in (spatial, cn):
        assert all(b["absolute_T_error"] < .4*a["absolute_T_error"] for a, b in zip(study, study[1:])), "Spatial refinement failed"
    ref = next(r["T_numerical"] for r in spatial if r["config"]["dx"] == .1)
    temporal = {r["config"]["dt"]: abs(r["T_numerical"]-ref) for r in rows if r["study"] == "dt"}
    assert temporal[.005] < .3*temporal[.01], "Small-step time refinement failed"
    boxes = {r["config"]["half_box"]: r for r in rows if r["study"] == "box"}
    assert abs(boxes[140]["T_numerical"]-boxes[220]["T_numerical"]) < 1e-7, "Large boxes disagree"
    assert boxes[100]["max_edge_probability"] > .01, "Undersized box did not trigger the diagnostic"
    for name in ("tunneling", "above_barrier"):
        r = json.loads((directory/f"{name}_summary.json").read_text())
        assert r["norm_drift"] < 1e-9
        assert r["absolute_T_error"] < .002
        assert r["barrier_remaining"] < 1e-7
        assert r["max_edge_probability"] < 1e-7
        assert r["plateau_change_last_10_percent"] < 1e-6
        for f in r["fits"].values():
            assert f["fit"]["relative_L2_residual"] < .002
            assert abs(f["fit"]["center"]-f["filtered_prediction"]["center"]) < .05
            assert abs(f["fit"]["sigma"]-f["filtered_prediction"]["sigma"]) < .02
            assert abs(f["fitted_velocity"]-f["filtered_prediction"]["velocity"]) < .001
    print("PASS: spatial/time refinement, box independence, boundary warning, norm, plateaus, Gaussian shapes, centers, widths and velocities")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("figures"))
    check(parser.parse_args().out)
