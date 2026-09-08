#!/usr/bin/env python3
"""Profile-parameter sweep for the native TCT current-redistribution operator.

This audit keeps the validated signed-current transfer matrix fixed while
varying only the native ``icd_source=4`` spatial-profile parameters:
``W_cd``, ``W_cd_shoulder``, and ``delta_cd``.  Every case starts from the
same frozen baseline and is sampled at equal physical times.

The result is a profile-screening audit, not a reactor-scale claim.  The
native no-source trajectory remains the primary gate reference and each
profile's source=4, amp=0 trajectory is retained as a same-mode null.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

# Keep this audit independent of the feedback-controller defaults.
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
OUT = REPO / "validation_runs/m3dc1_tct_profile_parameter_sweep"
RUN_ROOT = Path("/tmp/m3dc1_tct_profile_parameter_sweep_runs")

DT = 0.01
HORIZON_STEPS = 10
HORIZON_TIMES = (0.01, 0.02, 0.05, 0.10)
CURRENT_SOURCE = 4
AMPLITUDES = (-0.010, -0.005, -0.002, 0.000, 0.002, 0.005)
R0 = 10.0
Z0 = 1.0

# Frozen profile matrix.  The first row is the previously tested reference.
# Each other row changes one geometric parameter at a time so a useful
# response can be attributed to a single profile dimension.
PROFILES = (
    {"name": "default", "W_cd": 0.2805, "W_cd_shoulder": 0.2805, "delta_cd": 0.561},
    {"name": "center_narrow", "W_cd": 0.14025, "W_cd_shoulder": 0.2805, "delta_cd": 0.561},
    {"name": "center_wide", "W_cd": 0.561, "W_cd_shoulder": 0.2805, "delta_cd": 0.561},
    {"name": "shoulder_narrow", "W_cd": 0.2805, "W_cd_shoulder": 0.14025, "delta_cd": 0.561},
    {"name": "shoulder_wide", "W_cd": 0.2805, "W_cd_shoulder": 0.561, "delta_cd": 0.561},
    {"name": "shoulders_close", "W_cd": 0.2805, "W_cd_shoulder": 0.2805, "delta_cd": 0.2805},
    {"name": "shoulders_far", "W_cd": 0.2805, "W_cd_shoulder": 0.2805, "delta_cd": 1.122},
)

WIDTH_GATE_PCT = 0.02
JPK_GATE_PCT = 0.10
JPK_NOISE_PCT = 1e-6
ZERO_ABS_TOL = 1e-12


def write_json(path: Path, payload: object) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def rows_at(directory: Path) -> list[dict[str, float]]:
    rows = nfc.safe_extract(directory)
    if not rows:
        raise RuntimeError("no extracted rows in " + str(directory))
    for t in HORIZON_TIMES:
        if min(abs(row["time"] - t) for row in rows) > 1e-8:
            raise RuntimeError(
                "equal-time sample missing in {} at t={}".format(directory, t)
            )
    return rows


def amplitude_label(amp: float) -> str:
    if amp > 0:
        return "amp_p" + f"{abs(amp):.3f}".replace(".", "")
    if amp < 0:
        return "amp_m" + f"{abs(amp):.3f}".replace(".", "")
    return "amp_zero"


def set_profile(profile: dict[str, float | str]) -> None:
    """Set the globals consumed by the validated native input writer."""
    nfc.PROFILE_WIDTH = float(profile["W_cd"])
    nfc.SHOULDER_WIDTH = float(profile["W_cd_shoulder"])
    nfc.SHOULDER_DELTA = float(profile["delta_cd"])
    nfc.R0 = R0
    nfc.Z0 = Z0
    nfc.CURRENT_SOURCE = CURRENT_SOURCE


def metric_row(
    profile: dict[str, float | str],
    amp: float,
    t: float,
    row: dict[str, float],
    baseline: dict[str, float],
    zero: dict[str, float],
) -> dict[str, float | str | bool]:
    width = pta.pct(row["W_sheet"], baseline["W_sheet"])
    jpk = pta.pct(row["Jpk"], baseline["Jpk"])
    return {
        "profile": profile["name"],
        "W_cd": float(profile["W_cd"]),
        "W_cd_shoulder": float(profile["W_cd_shoulder"]),
        "delta_cd": float(profile["delta_cd"]),
        "case": f"{profile['name']}_{amplitude_label(amp)}",
        "source": CURRENT_SOURCE,
        "amp": amp,
        "time": t,
        "width_gain_pct": width,
        "Jpk_change_pct": jpk,
        "high_J_change_pct": pta.pct(row["Jint_high"], baseline["Jint_high"]),
        "center_change_pct": pta.pct(
            row["center_abs_current"], baseline["center_abs_current"]
        ),
        "shoulder_change_pct": pta.pct(
            row["shoulder_abs_current"], baseline["shoulder_abs_current"]
        ),
        "mode_width_gain_pct": pta.pct(row["W_sheet"], zero["W_sheet"]),
        "mode_Jpk_change_pct": pta.pct(row["Jpk"], zero["Jpk"]),
        "delta_Reconnected_Flux": row["Reconnected_Flux"] - baseline["Reconnected_Flux"],
        "delta_magnetic_energy": row["magnetic_energy"] - baseline["magnetic_energy"],
        "width_gate_pass": width > WIDTH_GATE_PCT,
        "current_gate_pass": jpk <= JPK_GATE_PCT,
        "safe_authority_candidate": width > WIDTH_GATE_PCT and jpk <= JPK_GATE_PCT,
    }


def null_delta(
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

    # Native baseline is shared by every profile and is run only once.
    baseline_dir = nfc.write_input(
        "baseline", 0, 0.0, 0, 0.0, DT * HORIZON_STEPS,
        nmax_steps=HORIZON_STEPS,
    )
    print("[profile-sweep] running baseline", flush=True)
    nfc.execute(baseline_dir)
    baseline_rows = rows_at(baseline_dir)

    all_rows: list[dict[str, float | str | bool]] = []
    raw_rows: list[dict[str, float | str]] = []
    profile_reports: list[dict] = []
    zero_results: dict[str, dict] = {}

    for profile in PROFILES:
        set_profile(profile)
        name = str(profile["name"])
        zero_dir = nfc.write_input(
            f"{name}_zero", CURRENT_SOURCE, 0.0, 0, 0.0, DT * HORIZON_STEPS,
            nmax_steps=HORIZON_STEPS,
        )
        print(f"[profile-sweep] running {name}_zero", flush=True)
        nfc.execute(zero_dir)
        zero_rows = rows_at(zero_dir)
        zero_max = null_delta(baseline_rows, zero_rows)
        zero_results[name] = {
            "max_abs_metric_delta": zero_max,
            "tolerance": ZERO_ABS_TOL,
            "pass": zero_max <= ZERO_ABS_TOL,
        }

        case_rows: list[dict[str, float | str | bool]] = []
        for amp in AMPLITUDES:
            label = amplitude_label(amp)
            case_dir = nfc.write_input(
                f"{name}_{label}", CURRENT_SOURCE, amp, 0, 0.0, DT * HORIZON_STEPS,
                nmax_steps=HORIZON_STEPS,
            )
            print(f"[profile-sweep] running {name}_{label}", flush=True)
            nfc.execute(case_dir)
            case_rows_at_horizon = rows_at(case_dir)
            for t in HORIZON_TIMES:
                row = nfc.nearest(case_rows_at_horizon, t)
                baseline = nfc.nearest(baseline_rows, t)
                zero = nfc.nearest(zero_rows, t)
                metrics = metric_row(profile, amp, t, row, baseline, zero)
                all_rows.append(metrics)
                case_rows.append(metrics)
                raw_rows.append({
                    "profile": name,
                    "W_cd": float(profile["W_cd"]),
                    "W_cd_shoulder": float(profile["W_cd_shoulder"]),
                    "delta_cd": float(profile["delta_cd"]),
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

        candidates = [r for r in case_rows if bool(r["safe_authority_candidate"])]
        best = max(candidates, key=lambda r: float(r["width_gain_pct"]), default=None)
        profile_reports.append({
            "profile": name,
            "geometry": {
                "W_cd": float(profile["W_cd"]),
                "W_cd_shoulder": float(profile["W_cd_shoulder"]),
                "delta_cd": float(profile["delta_cd"]),
            },
            "zero_equivalence": zero_results[name],
            "best_case": best,
            "max_width_gain_pct": max(
                (float(r["width_gain_pct"]) for r in case_rows), default=float("nan")
            ),
            "max_width_case": max(
                case_rows, key=lambda r: float(r["width_gain_pct"]), default=None
            ),
            "max_Jpk_change_pct": max(
                (float(r["Jpk_change_pct"]) for r in case_rows), default=float("nan")
            ),
            "width_gate_pass_any": any(bool(r["width_gate_pass"]) for r in case_rows),
            "current_gate_pass_all": all(bool(r["current_gate_pass"]) for r in case_rows),
        })

    candidates = [r for r in all_rows if bool(r["safe_authority_candidate"])]
    best = max(candidates, key=lambda r: float(r["width_gain_pct"]), default=None)
    classification = (
        "M3DC1_TCT_PROFILE_SWEEP_SAFE_CANDIDATE"
        if best is not None
        else "M3DC1_TCT_PROFILE_SWEEP_NO_SAFE_PROFILE_FOUND"
    )

    report = {
        "classification": classification,
        "claim_boundary": (
            "Native normalized M3D-C1 profile-parameter sweep only; no reactor "
            "stabilization, RF wave physics, lithium dimensional transfer, or "
            "experimental validation is implied."
        ),
        "audit": {
            "type": "equal_time_signed_current_transfer_profile_parameter_sweep",
            "dt": DT,
            "ntimemax": HORIZON_STEPS,
            "ntimepr": 1,
            "source": CURRENT_SOURCE,
            "amplitudes": list(AMPLITUDES),
            "horizon_times": list(HORIZON_TIMES),
            "profile_count": len(PROFILES),
            "profiles": list(PROFILES),
            "gates": {
                "width_threshold_pct": WIDTH_GATE_PCT,
                "Jpk_threshold_pct": JPK_GATE_PCT,
                "Jpk_noise_pct": JPK_NOISE_PCT,
                "zero_equivalence_tolerance": ZERO_ABS_TOL,
            },
            "comparison": (
                "Every profile/amplitude starts from the same frozen initial "
                "condition and is compared at equal physical times against the "
                "native source=0 baseline and that profile's source=4, amp=0 null."
            ),
        },
        "zero_equivalence": zero_results,
        "best_candidate": best,
        "profile_results": profile_reports,
        "cases": all_rows,
    }

    pta.write_csv(OUT / "profile_parameter_sweep.csv", all_rows)
    pta.write_csv(OUT / "profile_parameter_sweep_raw.csv", raw_rows)
    write_json(OUT / "profile_parameter_sweep_summary.json", report)
    (OUT / "runtime_provenance.txt").write_text(
        "repo={}\nsource={}\nbaseline={}\nexecutable={}\n"
        "executable_sha256={}\nrun_root={}\n"
        "dt={}\nntimemax={}\nntimepr=1\nsource={}\n"
        "amplitudes={}\nhorizon_times={}\nprofiles={}\n".format(
            REPO, SRC, BASE, EXE, pta.sha256_file(EXE), RUN_ROOT,
            DT, HORIZON_STEPS, CURRENT_SOURCE, list(AMPLITUDES),
            list(HORIZON_TIMES), list(PROFILES),
        )
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

