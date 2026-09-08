#!/usr/bin/env python3
"""Refined center-width transfer audit for native ``icd_source=4``.

The preceding profile sweep found its strongest result by narrowing ``W_cd``.
This audit refines that one parameter while holding the shoulder width and
separation at their validated reference values.  It preserves the same
equal-time native baseline, same-profile zero null, four horizon samples, and
current/width gates.

The additional per-unit-amplitude metrics are diagnostic only.  They expose
whether a small absolute response is caused by weak profile coupling,
amplitude normalization, or a nonlinear response region; they do not alter
the native M3D-C1 equations or rescale the actuator.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ["TCT_FEEDBACK_DT"] = "0.01"
os.environ["TCT_FEEDBACK_SEGMENT_STEPS"] = "1"
os.environ["TCT_FEEDBACK_MAX_SEGMENTS"] = "10"

import native_feedback_controller_audit as nfc
import pulse_train_audit as pta

REPO = Path("/home/ubuntu/work/openmc/sweep")
BASE = Path("/home/ubuntu/m3dc1_runs/TCT_MECHANISM_BASELINE")
SRC = Path("/home/ubuntu/M3DC1-official")
BUILD = SRC / "build-ubuntu-2d"
EXE = BUILD / "unstructured/m3dc1_2d"
OUT = REPO / "validation_runs/m3dc1_tct_center_width_refinement"
RUN_ROOT = Path("/tmp/m3dc1_tct_center_width_refinement_runs")

DT = 0.01
HORIZON_STEPS = 10
HORIZON_TIMES = (0.01, 0.02, 0.05, 0.10)
CURRENT_SOURCE = 4

# The prior best was W_cd=0.14025 at amp=-0.01.  The grid resolves that
# neighborhood and extends downward far enough to detect a turnover.
CENTER_WIDTHS = (0.080, 0.100, 0.120, 0.140, 0.160, 0.180, 0.200, 0.240, 0.2805)
AMPLITUDES = (-0.020, -0.010, -0.005, 0.000, 0.005)
SHOULDER_WIDTH = 0.2805
SHOULDER_DELTA = 0.561
R0 = 10.0
Z0 = 1.0

WIDTH_GATE_PCT = 0.02
JPK_GATE_PCT = 0.10
ZERO_ABS_TOL = 1e-12


def write_json(path: Path, payload: object) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def amplitude_label(amp: float) -> str:
    if amp > 0:
        return "amp_p" + f"{abs(amp):.3f}".replace(".", "")
    if amp < 0:
        return "amp_m" + f"{abs(amp):.3f}".replace(".", "")
    return "amp_zero"


def width_label(width: float) -> str:
    return "wc_" + f"{width:.4f}".replace(".", "p")


def set_center_width(width: float) -> None:
    nfc.PROFILE_WIDTH = width
    nfc.SHOULDER_WIDTH = SHOULDER_WIDTH
    nfc.SHOULDER_DELTA = SHOULDER_DELTA
    nfc.R0 = R0
    nfc.Z0 = Z0
    nfc.CURRENT_SOURCE = CURRENT_SOURCE


def rows_at(directory: Path) -> list[dict[str, float]]:
    rows = nfc.safe_extract(directory)
    if not rows:
        raise RuntimeError("no extracted rows in " + str(directory))
    for t in HORIZON_TIMES:
        if min(abs(row["time"] - t) for row in rows) > 1e-8:
            raise RuntimeError(f"equal-time sample missing in {directory} at t={t}")
    return rows


def zero_delta(
    baseline_rows: list[dict[str, float]], zero_rows: list[dict[str, float]]
) -> float:
    maximum = 0.0
    for t in HORIZON_TIMES:
        baseline = nfc.nearest(baseline_rows, t)
        zero = nfc.nearest(zero_rows, t)
        for key in (
            "W_sheet", "Jpk", "Jint_high", "center_abs_current",
            "shoulder_abs_current", "Reconnected_Flux", "magnetic_energy",
        ):
            maximum = max(maximum, abs(zero[key] - baseline[key]))
    return maximum


def metric_row(
    width: float,
    amp: float,
    t: float,
    row: dict[str, float],
    baseline: dict[str, float],
    zero: dict[str, float],
) -> dict[str, float | str | bool]:
    width_gain = pta.pct(row["W_sheet"], baseline["W_sheet"])
    jpk_change = pta.pct(row["Jpk"], baseline["Jpk"])
    return {
        "profile": width_label(width),
        "W_cd": width,
        "W_cd_shoulder": SHOULDER_WIDTH,
        "delta_cd": SHOULDER_DELTA,
        "case": f"{width_label(width)}_{amplitude_label(amp)}",
        "source": CURRENT_SOURCE,
        "amp": amp,
        "time": t,
        "width_gain_pct": width_gain,
        "Jpk_change_pct": jpk_change,
        "high_J_change_pct": pta.pct(row["Jint_high"], baseline["Jint_high"]),
        "center_change_pct": pta.pct(
            row["center_abs_current"], baseline["center_abs_current"]
        ),
        "shoulder_change_pct": pta.pct(
            row["shoulder_abs_current"], baseline["shoulder_abs_current"]
        ),
        "mode_width_gain_pct": pta.pct(row["W_sheet"], zero["W_sheet"]),
        "mode_Jpk_change_pct": pta.pct(row["Jpk"], zero["Jpk"]),
        "width_gain_per_amp_pct": width_gain / amp if amp else 0.0,
        "Jpk_change_per_amp_pct": jpk_change / amp if amp else 0.0,
        "delta_Reconnected_Flux": row["Reconnected_Flux"] - baseline["Reconnected_Flux"],
        "delta_magnetic_energy": row["magnetic_energy"] - baseline["magnetic_energy"],
        "width_gate_pass": width_gain > WIDTH_GATE_PCT,
        "current_gate_pass": jpk_change <= JPK_GATE_PCT,
        "safe_authority_candidate": (
            width_gain > WIDTH_GATE_PCT and jpk_change <= JPK_GATE_PCT
        ),
    }


def main() -> int:
    if not BASE.exists():
        raise FileNotFoundError(BASE)
    if not EXE.exists():
        raise FileNotFoundError(EXE)

    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    nfc.RUN_ROOT = RUN_ROOT
    nfc.OUT = OUT
    nfc.DT = DT
    nfc.SEGMENT_STEPS = 1
    nfc.SEGMENT_DURATION = DT

    pta.install_operator()
    nfc.install_current_redistribution_operator()
    pta.build()

    baseline_dir = nfc.write_input(
        "baseline", 0, 0.0, 0, 0.0, DT * HORIZON_STEPS,
        nmax_steps=HORIZON_STEPS,
    )
    print("[center-width] running baseline", flush=True)
    nfc.execute(baseline_dir)
    baseline_rows = rows_at(baseline_dir)

    all_rows: list[dict[str, float | str | bool]] = []
    raw_rows: list[dict[str, float | str]] = []
    width_results: list[dict] = []
    zero_results: dict[str, dict] = {}

    for width in CENTER_WIDTHS:
        set_center_width(width)
        profile = width_label(width)
        zero_dir = nfc.write_input(
            f"{profile}_zero", CURRENT_SOURCE, 0.0, 0, 0.0, DT * HORIZON_STEPS,
            nmax_steps=HORIZON_STEPS,
        )
        print(f"[center-width] running {profile}_zero", flush=True)
        nfc.execute(zero_dir)
        zero_rows = rows_at(zero_dir)
        null_max = zero_delta(baseline_rows, zero_rows)
        zero_results[profile] = {
            "max_abs_metric_delta": null_max,
            "tolerance": ZERO_ABS_TOL,
            "pass": null_max <= ZERO_ABS_TOL,
        }

        profile_rows: list[dict[str, float | str | bool]] = []
        for amp in AMPLITUDES:
            label = amplitude_label(amp)
            case_dir = nfc.write_input(
                f"{profile}_{label}", CURRENT_SOURCE, amp, 0, 0.0, DT * HORIZON_STEPS,
                nmax_steps=HORIZON_STEPS,
            )
            print(f"[center-width] running {profile}_{label}", flush=True)
            nfc.execute(case_dir)
            case_rows = rows_at(case_dir)
            for t in HORIZON_TIMES:
                row = nfc.nearest(case_rows, t)
                baseline = nfc.nearest(baseline_rows, t)
                zero = nfc.nearest(zero_rows, t)
                metrics = metric_row(width, amp, t, row, baseline, zero)
                all_rows.append(metrics)
                profile_rows.append(metrics)
                raw_rows.append({
                    "profile": profile,
                    "W_cd": width,
                    "W_cd_shoulder": SHOULDER_WIDTH,
                    "delta_cd": SHOULDER_DELTA,
                    "amp": amp,
                    "time": t,
                    "W_sheet": row["W_sheet"],
                    "Jpk": row["Jpk"],
                    "Jint_high": row["Jint_high"],
                    "center_abs_current": row["center_abs_current"],
                    "shoulder_abs_current": row["shoulder_abs_current"],
                    "Reconnected_Flux": row["Reconnected_Flux"],
                    "magnetic_energy": row["magnetic_energy"],
                })

        candidates = [r for r in profile_rows if bool(r["safe_authority_candidate"])]
        best = max(candidates, key=lambda r: float(r["width_gain_pct"]), default=None)
        width_results.append({
            "W_cd": width,
            "W_cd_shoulder": SHOULDER_WIDTH,
            "delta_cd": SHOULDER_DELTA,
            "profile": profile,
            "zero_equivalence": zero_results[profile],
            "best_case": best,
            "max_width_gain_pct": max(
                (float(r["width_gain_pct"]) for r in profile_rows), default=float("nan")
            ),
            "max_width_case": max(
                profile_rows, key=lambda r: float(r["width_gain_pct"]), default=None
            ),
            "max_Jpk_change_pct": max(
                (float(r["Jpk_change_pct"]) for r in profile_rows), default=float("nan")
            ),
            "max_abs_Jpk_change_pct": max(
                (abs(float(r["Jpk_change_pct"])) for r in profile_rows), default=float("nan")
            ),
            "width_gate_pass_any": any(bool(r["width_gate_pass"]) for r in profile_rows),
            "current_gate_pass_all": all(bool(r["current_gate_pass"]) for r in profile_rows),
        })

    candidates = [r for r in all_rows if bool(r["safe_authority_candidate"])]
    best = max(candidates, key=lambda r: float(r["width_gain_pct"]), default=None)

    # Local finite-difference slopes use the symmetric +/-0.005 pair.  These
    # are reported per unit amp and are not used to change the pass/fail gate.
    local_gain = []
    for width in CENTER_WIDTHS:
        profile = width_label(width)
        for t in HORIZON_TIMES:
            pair = [
                r for r in all_rows
                if (
                    r["profile"] == profile
                    and abs(float(r["time"]) - t) < 1e-12
                    and abs(abs(float(r["amp"])) - 0.005) < 1e-12
                )
            ]
            plus = next((r for r in pair if float(r["amp"]) > 0), None)
            minus = next((r for r in pair if float(r["amp"]) < 0), None)
            if plus is not None and minus is not None:
                local_gain.append({
                    "profile": profile,
                    "W_cd": width,
                    "time": t,
                    "symmetric_width_gain_slope_pct_per_amp": (
                        float(plus["width_gain_pct"]) - float(minus["width_gain_pct"])
                    ) / 0.010,
                    "symmetric_Jpk_slope_pct_per_amp": (
                        float(plus["Jpk_change_pct"]) - float(minus["Jpk_change_pct"])
                    ) / 0.010,
                })

    report = {
        "classification": (
            "M3DC1_TCT_CENTER_WIDTH_REFINEMENT_SAFE_CANDIDATE"
            if best is not None
            else "M3DC1_TCT_CENTER_WIDTH_REFINEMENT_NO_SAFE_CANDIDATE"
        ),
        "claim_boundary": (
            "Native normalized M3D-C1 center-width refinement audit only; no "
            "reactor stabilization, RF wave physics, lithium dimensional "
            "transfer, or experimental validation is implied."
        ),
        "audit": {
            "type": "equal_time_signed_current_transfer_center_width_refinement",
            "dt": DT,
            "ntimemax": HORIZON_STEPS,
            "ntimepr": 1,
            "source": CURRENT_SOURCE,
            "center_widths": list(CENTER_WIDTHS),
            "amplitudes": list(AMPLITUDES),
            "horizon_times": list(HORIZON_TIMES),
            "fixed_profile": {
                "R_0cd": R0,
                "Z_0cd": Z0,
                "W_cd_shoulder": SHOULDER_WIDTH,
                "delta_cd": SHOULDER_DELTA,
            },
            "gates": {
                "width_threshold_pct": WIDTH_GATE_PCT,
                "Jpk_threshold_pct": JPK_GATE_PCT,
                "zero_equivalence_tolerance": ZERO_ABS_TOL,
            },
            "comparison": (
                "Every center width and amplitude starts from the same frozen "
                "initial condition and is compared at equal physical times "
                "against the native source=0 baseline and its same-profile "
                "source=4, amp=0 null."
            ),
            "prior_model_context": (
                "Interpret alongside prior BOUT++ reduced-MHD and OpenMC-style "
                "screening evidence. Differences may reflect actuator/source "
                "normalization, observable definition, time horizon, or model "
                "fidelity; this audit does not merge those evidence levels."
            ),
        },
        "zero_equivalence": zero_results,
        "best_candidate": best,
        "center_width_results": width_results,
        "local_symmetric_gain": local_gain,
        "cases": all_rows,
    }

    pta.write_csv(OUT / "center_width_refinement.csv", all_rows)
    pta.write_csv(OUT / "center_width_refinement_raw.csv", raw_rows)
    pta.write_csv(OUT / "center_width_local_gain.csv", local_gain)
    write_json(OUT / "center_width_refinement_summary.json", report)
    (OUT / "runtime_provenance.txt").write_text(
        "repo={}\nsource={}\nbaseline={}\nexecutable={}\n"
        "executable_sha256={}\nrun_root={}\n"
        "dt={}\nntimemax={}\nntimepr=1\nsource={}\n"
        "center_widths={}\namplitudes={}\nhorizon_times={}\n"
        "shoulder_width={}\nshoulder_separation={}\n".format(
            REPO, SRC, BASE, EXE, pta.sha256_file(EXE), RUN_ROOT,
            DT, HORIZON_STEPS, CURRENT_SOURCE, list(CENTER_WIDTHS),
            list(AMPLITUDES), list(HORIZON_TIMES), SHOULDER_WIDTH,
            SHOULDER_DELTA,
        )
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

