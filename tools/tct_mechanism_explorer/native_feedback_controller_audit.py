#!/usr/bin/env python3
"""Native M3D-C1 state-feedback controller audit.

This runner attempts a real supervisory loop:
  run a short native segment -> read W_sheet/Jpk -> select the next
  allow-listed current-profile command -> restart from C1.h5.

It never changes eta, nu, GEM epsilon, mesh, equilibrium, or solver equations.
If the installed M3D-C1 executable cannot restart from C1.h5 with changed
control inputs, the result is explicitly marked unresolved rather than called
closed-loop evidence.
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import time
from pathlib import Path

import pulse_train_audit as pta

REPO = Path("/home/ubuntu/work/openmc/sweep")
BASE = Path("/home/ubuntu/m3dc1_runs/TCT_MECHANISM_BASELINE")
SRC = Path("/home/ubuntu/M3DC1-official")
BUILD = SRC / "build-ubuntu-2d"
EXE = BUILD / "unstructured/m3dc1_2d"
OUT = REPO / "validation_runs/m3dc1_tct_native_feedback"
RUN_ROOT = Path("/tmp/m3dc1_tct_native_feedback_runs")

DT = float(os.environ.get("TCT_FEEDBACK_DT", "0.01"))
# One restart per native output step lets the controller see a developing
# current excursion before a coarse 0.05-time-unit segment commits it.
SEGMENT_STEPS = int(os.environ.get("TCT_FEEDBACK_SEGMENT_STEPS", "1"))
SEGMENT_DURATION = DT * SEGMENT_STEPS
MAX_SEGMENTS = int(os.environ.get("TCT_FEEDBACK_MAX_SEGMENTS", "40"))
NTIMEPR = 1
# Native M3D-C1 restart controls.  irestart selects restart mode;
# iwrite_restart causes the preceding segment to emit restart state.
RESTART_KEY = "irestart"
IRESTART_SLICE_KEY = "irestart_slice"
IWRITE_RESTART_KEY = "iwrite_restart"

# Control policy thresholds are control-layer values, not physics changes.
THINNING_RATE_THRESHOLD = -0.05
JPK_GROWTH_THRESHOLD = 0.05
BIAS_AMP = -0.002
AGGRESSIVE_AMP = -0.02
HOLD_AMP = -0.005
# Absolute-current safety is a control-layer constraint. It does not change
# M3D-C1 physics; it prevents an aggressive command from remaining active
# after the controlled trajectory exceeds the native no-actuator current.
JPK_ABSOLUTE_GUARD_FRACTION = 0.10
JPK_PREDICTIVE_TRIGGER_FRACTION = 0.05
JPK_REARM_FRACTION = 0.03
SAFE_RELEASE_AMP = 0.0
# Opposite sign to the previously tested aggressive command.  This is an
# allow-listed actuator command under audit, not a solver-equation change.
ACTIVE_BRAKE_AMP = 0.005
PROFILE_WIDTH = 0.2805
SHOULDER_WIDTH = 0.2805
SHOULDER_DELTA = 0.561
CURRENT_SOURCE = 4
R0 = 10.0
Z0 = 1.0

NATIVE_CONTROL_KEYS = {
    "dt", "ntimemax", "ntimepr", "irestart", "irestart_slice", "iwrite_restart",
    "imag_control", "mag_ctrl_amp",
    "icd_source", "J_0cd", "R_0cd", "Z_0cd", "W_cd", "W_cd_shoulder",
    "delta_cd", "cd_t_on", "cd_t_ramp", "cd_t_off",
}



def install_current_redistribution_operator() -> bool:
    """Install the previously validated icd_source=4 TCT profile.

    The official icd_source=1 branch is a single Gaussian plus a uniform
    subtraction term.  It is not the center/shoulder redistribution used by
    the earlier TCT authority rung.  This installer restores that separate,
    net-current-neutral source without changing solver physics.
    """
    modules = SRC / "unstructured/M3Dmodules.f90"
    inputf = SRC / "unstructured/input.f90"
    transport = SRC / "unstructured/transport.f90"
    for path in (modules, inputf, transport):
        if not path.exists():
            raise FileNotFoundError(path)
    changed = False

    text = modules.read_text()
    missing = []
    declarations = {
        "w_cd_shoulder": "  real :: w_cd_shoulder  ! width of shoulder Gaussians for neutral redistribution",
        "cd_t_on": "  real :: cd_t_on        ! current-drive turn-on time",
        "cd_t_ramp": "  real :: cd_t_ramp      ! current-drive smooth ramp duration",
        "cd_t_off": "  real :: cd_t_off       ! current-drive turn-off time",
    }
    for name, declaration in declarations.items():
        if not re.search(rf"\b{re.escape(name)}\b", text, re.I):
            missing.append(declaration)
    if missing:
        anchor = re.search(r"^\s*real\s*::\s*delta_cd\b[^\n]*$", text, re.I | re.M)
        if not anchor:
            raise RuntimeError("current-drive module declaration anchor not found")
        addition = anchor.group(0) + "\n" + "\n".join(missing)
        text = text[:anchor.start()] + addition + text[anchor.end():]
        modules.write_text(text)
        changed = True

    text = inputf.read_text()
    if '"W_cd_shoulder"' not in text:
        anchor = re.search(
            r'^\s*call\s+add_var_double\("delta_cd"[^\n]*\n(?:[^\n]*\n){0,2}',
            text, re.I | re.M,
        )
        if not anchor:
            raise RuntimeError("delta_cd input registration anchor not found")
        regs = anchor.group(0) + (
            '  call add_var_double("W_cd_shoulder", w_cd_shoulder, 0., &\n'
            '       "shoulder width for neutral center-plus-shoulder cd source", source_grp)\n'
            '  call add_var_double("cd_t_on", cd_t_on, 0., &\n'
            '       "time when current drive turns on", source_grp)\n'
            '  call add_var_double("cd_t_ramp", cd_t_ramp, 0., &\n'
            '       "smooth current-drive ramp duration", source_grp)\n'
            '  call add_var_double("cd_t_off", cd_t_off, 1.e30, &\n'
            '       "time when current drive turns off", source_grp)\n'
        )
        text = text[:anchor.start()] + regs + text[anchor.end():]
        inputf.write_text(text)
        changed = True

    text = transport.read_text()
    if "icd_source.eq.4" not in text:
        start = text.find("function cd_func")
        end = text.find("cd_func = temp", start)
        if start < 0 or end < 0:
            raise RuntimeError("cd_func source block not found")
        insert_at = text.rfind("  endif", start, end)
        if insert_at < 0:
            raise RuntimeError("cd_func closing endif not found")
        block = """  else if(icd_source.eq.4) then
     cd_w_center = max(w_cd, 1.e-30)
     if(w_cd_shoulder.gt.0.) then
        cd_w_sh = w_cd_shoulder
     else
        cd_w_sh = cd_w_center
     end if
     cd_sep = abs(delta_cd)
     temp79a = 0.
     temp79b = 0.
     do j=1,npoints
        call magnetic_region(pst79(j,OP_1),pst79(j,OP_DR),pst79(j,OP_DZ), &
             x_79(j),z_79(j),iregion)
        if(iregion.eq.REGION_PLASMA) then
           temp79a(j) = -exp( -(x_79(j)-R_0cd)**2/cd_w_center**2 &
                - (z_79(j)-Z_0cd)**2/cd_w_center**2 ) &
                + 0.5*exp( -(x_79(j)-R_0cd)**2/cd_w_sh**2 &
                - (z_79(j)-(Z_0cd-cd_sep))**2/cd_w_sh**2 ) &
                + 0.5*exp( -(x_79(j)-R_0cd)**2/cd_w_sh**2 &
                - (z_79(j)-(Z_0cd+cd_sep))**2/cd_w_sh**2 )
           temp79b(j) = 1.
        end if
     enddo
     cd_area = real(int1(temp79b))
     if(cd_area.gt.0.) then
        cd_net = real(int1(temp79a))/cd_area
        do j=1,npoints
           if(real(temp79b(j)).gt.0.) temp79a(j) = temp79a(j) - cd_net
        enddo
     end if
     temp79a = cd_gate * J_0cd * temp79a
     temp = temp + intx2(mu79(:,:,OP_1),temp79a)
"""
        # Add the gate and neutral-profile locals to the native cd_func.
        decl = re.search(
            r"^\s*real,\s*allocatable\s*::\s*xvals\(:\),\s*yvals\(:\)\s*$",
            text[start:end], re.I | re.M,
        )
        if not decl:
            raise RuntimeError("cd_func local declaration anchor not found")
        decl_end = start + decl.end()
        local_add = "\n  real :: cd_gate, cd_tau, cd_w_center, cd_w_sh, cd_sep\n  real :: cd_area, cd_net"
        text = text[:decl_end] + local_add + text[decl_end:]
        end += len(local_add)
        temp_match = re.search(r"^\s*temp\s*=\s*0\.\s*$", text[decl_end:end], re.I | re.M)
        if not temp_match:
            raise RuntimeError("cd_func temp initialization anchor not found")
        gate_pos = decl_end + temp_match.end()
        gate = (
            "  if(time.lt.cd_t_on .or. time.ge.cd_t_off) then\n"
            "     cd_gate = 0.\n"
            "  else if(cd_t_ramp.gt.0. .and. time.lt.cd_t_on+cd_t_ramp) then\n"
            "     cd_tau = (time-cd_t_on)/cd_t_ramp\n"
            "     cd_gate = cd_tau*cd_tau*(3. - 2.*cd_tau)\n"
            "  else\n"
            "     cd_gate = 1.\n"
            "  end if\n"
        )
        text = text[:gate_pos] + "\n" + gate + text[gate_pos:]
        end += len(gate) + 1
        insert_at = text.rfind("  endif", start, end)
        if insert_at < 0:
            raise RuntimeError("cd_func closing endif not found")
        text = text[:insert_at] + block + text[insert_at:]
        transport.write_text(text)
        changed = True
    return changed

def write_input(name: str, source: int, amp: float, restart: int,
                t_on: float, t_off: float, nmax_steps: int | None = None) -> Path:
    d = RUN_ROOT / name
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    for item in pta.COPY_NAMES:
        src = BASE / item
        if src.is_symlink():
            (d / item).symlink_to(src.readlink())
        elif src.exists():
            shutil.copy2(src, d / item)
    text = (BASE / "C1input").read_text()
    updates = {
        "dt": f"{DT:.10g}",
        # ntimemax is an absolute timestep ceiling; restarted segments must
        # extend it beyond the timestep already stored in the restart state.
        "ntimemax": str(nmax_steps if nmax_steps is not None else SEGMENT_STEPS),
        "ntimepr": str(NTIMEPR),
        RESTART_KEY: str(restart),
        IRESTART_SLICE_KEY: "-1",
        IWRITE_RESTART_KEY: "1",
        "imag_control": "0",
        "mag_ctrl_amp": "0.0",
        "icd_source": str(source),
        "J_0cd": f"{amp:.10g}",
        "R_0cd": f"{R0:.10g}",
        "Z_0cd": f"{Z0:.10g}",
        "W_cd": f"{PROFILE_WIDTH:.10g}",
        "W_cd_shoulder": f"{SHOULDER_WIDTH:.10g}",
        "delta_cd": f"{SHOULDER_DELTA:.10g}",
        "cd_t_on": f"{t_on:.10g}",
        "cd_t_ramp": "0.0",
        "cd_t_off": f"{t_off:.10g}",
    }
    for key, value in updates.items():
        if key not in NATIVE_CONTROL_KEYS:
            raise RuntimeError(f"non-allow-listed control key: {key}")
        text = pta.replace_or_add(text, key, value)
    # Always write native restart controls explicitly. irestart=0 starts a
    # fresh run; irestart=1 consumes restart state copied from the prior segment.
    for key, value in {
        RESTART_KEY: str(restart),
        IRESTART_SLICE_KEY: "-1",
        IWRITE_RESTART_KEY: "1",
    }.items():
        text = pta.replace_or_add(text, key, value)
    (d / "C1input").write_text(text)
    launch = f'''#!/usr/bin/env bash
set -euo pipefail
export TMPDIR="${TMPDIR:-/tmp/tct-$USER}"
export OMPI_MCA_orte_tmpdir_base="${OMPI_MCA_orte_tmpdir_base:-$TMPDIR}"
mkdir -p "$TMPDIR" "$OMPI_MCA_orte_tmpdir_base"
source "$HOME/spack/share/spack/setup-env.sh"
spack env activate m3dc1-deps
cd "{d}"
set +e
timeout 1200s mpirun --oversubscribe -n 1 "{EXE}" -pc_factor_mat_solver_type mumps > C1stdout 2> launcher.stderr
rc=$?
set -e
printf 'return_code=%s\\n' "$rc" > run_status.txt
exit "$rc"
'''
    (d / "launch_command.sh").write_text(launch)
    (d / "launch_command.sh").chmod(0o755)
    return d


def copy_restart_state(previous: Path, current: Path) -> None:
    # C1.h5 is the primary restart state. Preserve any auxiliary restart files
    # emitted by a particular M3D-C1 checkout without assuming their names.
    for item in previous.iterdir():
        if item.name in {"C1input", "launch_command.sh", "C1stdout", "launcher.stderr",
                         "run_status.txt", "wrapper_stdout.log", "elapsed_seconds.txt"}:
            continue
        target = current / item.name
        if item.is_dir():
            # Some builds place restart metadata in a directory.
            shutil.copytree(item, target, dirs_exist_ok=True)
        elif item.is_file():
            shutil.copy2(item, target)


def execute(d: Path) -> None:
    t0 = time.time()
    p = pta.sh(["bash", "launch_command.sh"], cwd=d)
    (d / "wrapper_stdout.log").write_text(p.stdout)
    (d / "elapsed_seconds.txt").write_text(f"{time.time()-t0:.6f}\n")
    if p.returncode:
        raise RuntimeError(f"{d.name} failed rc={p.returncode}\n{p.stdout[-5000:]}")
    if not (d / "C1.h5").exists():
        stdout_tail = (d / "C1stdout").read_text(errors="replace")[-6000:] if (d / "C1stdout").exists() else "missing C1stdout"
        stderr_tail = (d / "launcher.stderr").read_text(errors="replace")[-6000:] if (d / "launcher.stderr").exists() else "missing launcher.stderr"
        status = (d / "run_status.txt").read_text(errors="replace") if (d / "run_status.txt").exists() else "missing run_status.txt"
        raise RuntimeError(f"{d.name} produced no C1.h5; status={status}\nC1stdout tail:\n{stdout_tail}\nlauncher.stderr tail:\n{stderr_tail}")


def safe_extract(d: Path) -> list[dict[str, float]]:
    try:
        return pta.extract(d)
    except Exception as exc:
        raise RuntimeError(f"cannot extract {d}: {exc}") from exc


def policy(
    previous: dict[str, float] | None,
    current: dict[str, float],
    baseline_ref: dict[str, float],
    prior_state: str,
) -> dict[str, float | int | str | bool]:
    """Select a bounded command with predictive current protection.

    The previous controller used a coarse segment and a release-only response.
    Here the restart cadence is one native output step by default.  The policy
    predicts the next Jpk from the measured local slope and applies an
    opposite-sign ACTIVE_BRAKE before or immediately after the absolute guard.
    ACTIVE_BRAKE latches until the measured and predicted current are within the
    rearm band.  All comparisons remain against the equal-time native baseline.
    """
    baseline_jpk = max(abs(baseline_ref["Jpk"]), 1e-300)
    jpk_vs_baseline = current["Jpk"] / baseline_jpk - 1.0
    if previous is None:
        return {
            "state": "BIAS", "source": CURRENT_SOURCE, "amp": BIAS_AMP,
            "guard_triggered": False, "predictive_triggered": False,
            "reason": "initial_bias",
            "jpk_vs_baseline_pct": 100.0 * jpk_vs_baseline,
            "predicted_jpk_vs_baseline_pct": 100.0 * jpk_vs_baseline,
        }

    dt = max(current["time"] - previous["time"], DT)
    dw = (current["W_sheet"] - previous["W_sheet"]) / dt
    dj = (current["Jpk"] - previous["Jpk"]) / max(abs(previous["Jpk"]), 1e-300) / dt
    predicted_jpk = current["Jpk"] + (current["Jpk"] - previous["Jpk"])
    predicted_vs_baseline = predicted_jpk / baseline_jpk - 1.0
    guard_triggered = jpk_vs_baseline > JPK_ABSOLUTE_GUARD_FRACTION
    predictive_triggered = (
        jpk_vs_baseline > JPK_PREDICTIVE_TRIGGER_FRACTION
        or predicted_vs_baseline > JPK_ABSOLUTE_GUARD_FRACTION
    )
    common = {
        "dW_dt": dw,
        "dJpk_dt_fraction": dj,
        "guard_triggered": guard_triggered,
        "predictive_triggered": predictive_triggered,
        "jpk_vs_baseline_pct": 100.0 * jpk_vs_baseline,
        "predicted_jpk_vs_baseline_pct": 100.0 * predicted_vs_baseline,
    }

    if prior_state == "ACTIVE_BRAKE" and max(
        jpk_vs_baseline, predicted_vs_baseline
    ) > JPK_REARM_FRACTION:
        return {
            "state": "ACTIVE_BRAKE", "source": CURRENT_SOURCE,
            "amp": ACTIVE_BRAKE_AMP, "reason": "jpk_brake_latched",
            **common,
        }
    if predictive_triggered:
        reason = (
            "jpk_absolute_guard"
            if guard_triggered else "predictive_jpk_guard"
        )
        return {
            "state": "ACTIVE_BRAKE", "source": CURRENT_SOURCE,
            "amp": ACTIVE_BRAKE_AMP, "reason": reason,
            **common,
        }
    if dw <= THINNING_RATE_THRESHOLD or dj >= JPK_GROWTH_THRESHOLD:
        return {
            "state": "AGGRESSIVE", "source": CURRENT_SOURCE,
            "amp": AGGRESSIVE_AMP, "reason": "precursor_thinning_or_current_growth",
            **common,
        }
    if dw >= 0.0:
        return {
            "state": "HOLD", "source": CURRENT_SOURCE,
            "amp": HOLD_AMP, "reason": "width_recovering",
            **common,
        }
    return {
        "state": "BIAS", "source": CURRENT_SOURCE,
        "amp": BIAS_AMP, "reason": "no_precursor_trigger",
        **common,
    }


def nearest(rows: list[dict[str, float]], t: float) -> dict[str, float]:
    return min(rows, key=lambda row: abs(row["time"] - t))


def deltas(control: list[dict[str, float]], baseline: list[dict[str, float]]) -> list[dict[str, float]]:
    out = []
    for row in control:
        ref = nearest(baseline, row["time"])
        out.append({
            "time": row["time"],
            "width_gain_pct": pta.pct(row["W_sheet"], ref["W_sheet"]),
            "Jpk_change_pct": pta.pct(row["Jpk"], ref["Jpk"]),
            "high_J_change_pct": pta.pct(row["Jint_high"], ref["Jint_high"]),
            "center_change_pct": pta.pct(row["center_abs_current"], ref["center_abs_current"]),
            "shoulder_change_pct": pta.pct(row["shoulder_abs_current"], ref["shoulder_abs_current"]),
            "delta_Reconnected_Flux": row["Reconnected_Flux"] - ref["Reconnected_Flux"],
            "delta_magnetic_energy": row["magnetic_energy"] - ref["magnetic_energy"],
        })
    return out


def write_json_artifact(path: Path, payload: object) -> None:
    """Write one complete JSON document atomically.

    Atomic replacement prevents a stale/partial command history from being
    concatenated with a later run if the wrapper is interrupted.
    """
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def main() -> int:
    if not BASE.exists():
        raise FileNotFoundError(BASE)
    if not EXE.exists():
        raise FileNotFoundError(EXE)
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    pta.install_operator()
    install_current_redistribution_operator()
    pta.build()

    report: dict = {
        "controller": {
            "type": "native_segmented_state_feedback",
            "observable": ["W_sheet", "Jpk", "dW_dt", "dJpk_dt"],
            "states": ["BIAS", "AGGRESSIVE", "HOLD", "ACTIVE_BRAKE"],
            "restart_key": RESTART_KEY,
            "actuator": "icd_source=4 net-current-neutral center-plus-shoulder redistribution",
            "profile": {"R_0cd": R0, "Z_0cd": Z0, "W_cd": PROFILE_WIDTH, "W_cd_shoulder": SHOULDER_WIDTH, "delta_cd": SHOULDER_DELTA},
            "segment_duration": SEGMENT_DURATION,
            "max_segments": MAX_SEGMENTS,
            "policy": {
                "thinning_rate_threshold": THINNING_RATE_THRESHOLD,
                "jpk_growth_threshold": JPK_GROWTH_THRESHOLD,
                "bias_amp": BIAS_AMP,
                "aggressive_amp": AGGRESSIVE_AMP,
                "hold_amp": HOLD_AMP,
                "jpk_absolute_guard_fraction": JPK_ABSOLUTE_GUARD_FRACTION,
                "jpk_predictive_trigger_fraction": JPK_PREDICTIVE_TRIGGER_FRACTION,
                "jpk_rearm_fraction": JPK_REARM_FRACTION,
                "safe_release_amp": SAFE_RELEASE_AMP,
                "active_brake_amp": ACTIVE_BRAKE_AMP,
            },
        },
        "claim_boundary": (
            "Native normalized M3D-C1 supervisory-controller audit only. "
            "No RF wave physics, lithium dimensional transfer, or reactor-scale "
            "stabilization is implied."
        ),
    }

    # Continuous uncontrolled comparison.
    horizon_steps = MAX_SEGMENTS * SEGMENT_STEPS
    horizon_time = MAX_SEGMENTS * SEGMENT_DURATION
    baseline_dir = write_input(
        "baseline_continuous", 0, 0.0, 0, 0.0, horizon_time,
        nmax_steps=horizon_steps,
    )
    print("[native-feedback] running baseline_continuous", flush=True)
    execute(baseline_dir)
    baseline_rows = safe_extract(baseline_dir)

    # Zero-actuation restart seed checks whether the current-drive path is
    # equivalent to no source before any feedback claim is made.
    zero_dir = write_input(
        "zero_current_drive", CURRENT_SOURCE, 0.0, 0, 0.0, horizon_time,
        nmax_steps=horizon_steps,
    )
    print("[native-feedback] running zero_current_drive", flush=True)
    execute(zero_dir)
    zero_rows = safe_extract(zero_dir)
    max_zero = 0.0
    for b, z in zip(baseline_rows, zero_rows):
        for key in ("W_sheet", "Jpk", "Jint_high", "Reconnected_Flux", "magnetic_energy"):
            max_zero = max(max_zero, abs(z[key] - b[key]))
    report["zero_equivalence"] = {
        "comparison": "icd_source=4,J_0cd=0 versus native icd_source=0 baseline over full horizon",
        "max_abs_metric_delta": max_zero,
        "tolerance": 1e-12,
        "pass": max_zero <= 1e-12,
        "interpretation": (
            "A nonzero delta is a source-mode side effect, not controller authority. "
            "Current-drive authority must compare controlled segments against this "
            "same-mode zero case and against the native no-actuator baseline."
        ),
    }

    # First controlled segment starts from the same frozen initial condition.
    previous_dir = write_input("segment_000", CURRENT_SOURCE, BIAS_AMP, 0, 0.0, SEGMENT_DURATION)
    print("[native-feedback] running segment_000", flush=True)
    execute(previous_dir)
    segment_rows = safe_extract(previous_dir)
    control_rows = list(segment_rows)
    command_log: list[dict] = [{
        "segment": 0, "directory": str(previous_dir), "state": "BIAS",
        "source": CURRENT_SOURCE, "amp": BIAS_AMP, "t_on": 0.0,
        "t_off": SEGMENT_DURATION,
        "restart": 0,
        "guard_triggered": False,
        "predictive_triggered": False,
        "reason": "initial_bias",
    }]

    restart_ok = True
    restart_error = None
    previous_state = segment_rows[-1]
    last_state = "BIAS"
    for segment in range(1, MAX_SEGMENTS):
        start_time = float(previous_state["time"])
        prior_row = control_rows[-2] if len(control_rows) >= 2 else None
        baseline_ref = nearest(baseline_rows, start_time)
        decision = policy(prior_row, previous_state, baseline_ref, last_state)
        current_dir = write_input(
            f"segment_{segment:03d}", int(decision["source"]), float(decision["amp"]),
            1, start_time, start_time + SEGMENT_DURATION,
            nmax_steps=int(round((start_time + SEGMENT_DURATION) / DT)),
        )
        copy_restart_state(previous_dir, current_dir)
        print(f"[native-feedback] running segment_{segment:03d} state={decision['state']}", flush=True)
        try:
            execute(current_dir)
            rows = safe_extract(current_dir)
            if not rows:
                raise RuntimeError(
                    f"restart produced no extracted rows: start={start_time}"
                )
            # A valid restart output includes the loaded state itself as the
            # first row at exactly start_time. Require a later emitted row,
            # then remove that duplicated seed before appending the segment.
            advanced_rows = [
                row for row in rows if row["time"] > start_time + 1e-8
            ]
            if not advanced_rows:
                raise RuntimeError(
                    f"restart did not advance physical time: start={start_time}, "
                    f"first={rows[0]['time']}, last={rows[-1]['time']}"
                )
            rows = advanced_rows
        except Exception as exc:
            restart_ok = False
            restart_error = str(exc)
            command_log.append({
                "segment": segment, "directory": str(current_dir),
                "state": decision["state"], "source": decision["source"],
                "amp": decision["amp"], "t_on": start_time,
                "t_off": start_time + SEGMENT_DURATION,
            "ntimemax": int(round((start_time + SEGMENT_DURATION) / DT)),
            "restart": 1,
                "error": str(exc),
            })
            break
        control_rows.extend(rows)
        previous_dir = current_dir
        previous_state = rows[-1]
        command_log.append({
            "segment": segment, "directory": str(current_dir),
            "state": decision["state"], "source": decision["source"],
            "amp": decision["amp"], "t_on": start_time,
            "t_off": start_time + SEGMENT_DURATION, "restart": 1,
            "dW_dt": decision.get("dW_dt"),
            "dJpk_dt_fraction": decision.get("dJpk_dt_fraction"),
            "jpk_vs_baseline_pct": decision.get("jpk_vs_baseline_pct"),
            "predicted_jpk_vs_baseline_pct": decision.get(
                "predicted_jpk_vs_baseline_pct"
            ),
            "guard_triggered": decision.get("guard_triggered", False),
            "predictive_triggered": decision.get("predictive_triggered", False),
            "reason": decision.get("reason"),
        })
        last_state = str(decision["state"])

    report["restart"] = {"pass": restart_ok, "error": restart_error, "segments_completed": len(command_log)}
    report["command_history"] = command_log

    if not restart_ok:
        report["classification"] = "M3DC1_NATIVE_FEEDBACK_RESTART_UNRESOLVED"
    else:
        # Compare at equal physical times. The source-mode zero trajectory is
        # the primary actuator comparison; the native no-source baseline is
        # retained as the overall experiment reference.
        out = deltas(control_rows, baseline_rows)
        mode_out = deltas(control_rows, zero_rows)
        pta.write_csv(OUT / "native_feedback_deltas.csv", out)
        pta.write_csv(OUT / "source_mode_zero_deltas.csv", mode_out)
        mode_max_width = max((r["width_gain_pct"] for r in mode_out), default=-math.inf)
        mode_max_jpk = max((r["Jpk_change_pct"] for r in mode_out), default=math.inf)
        report["classification"] = (
            "M3DC1_NATIVE_FEEDBACK_CONTROLLER_AUTHORITY_PASS"
            if mode_max_width > 0.02 and mode_max_jpk <= 0.10
            else "M3DC1_NATIVE_FEEDBACK_CONTROLLER_NO_AUTHORITY_FOUND"
        )
        report["response"] = {
            "comparison": "native icd_source=0 baseline",
            "max_width_gain_pct": max((r["width_gain_pct"] for r in out), default=math.nan),
            "mean_width_gain_pct": sum(r["width_gain_pct"] for r in out) / max(len(out), 1),
            "max_Jpk_change_pct": max((r["Jpk_change_pct"] for r in out), default=math.nan),
            "final_width_gain_pct": out[-1]["width_gain_pct"] if out else math.nan,
            "final_Jpk_change_pct": out[-1]["Jpk_change_pct"] if out else math.nan,
        }
        report["response_vs_source_mode_zero"] = {
            "comparison": "icd_source=4,J_0cd=0",
            "max_width_gain_pct": mode_max_width,
            "mean_width_gain_pct": sum(r["width_gain_pct"] for r in mode_out) / max(len(mode_out), 1),
            "max_Jpk_change_pct": mode_max_jpk,
            "final_width_gain_pct": mode_out[-1]["width_gain_pct"] if mode_out else math.nan,
            "final_Jpk_change_pct": mode_out[-1]["Jpk_change_pct"] if mode_out else math.nan,
        }
    write_json_artifact(OUT / "command_history.json", command_log)
    write_json_artifact(OUT / "native_feedback_summary.json", report)
    (OUT / "runtime_provenance.txt").write_text(
        f"repo={REPO}\nsource={SRC}\nbaseline={BASE}\nexecutable={EXE}\n"
        f"executable_sha256={pta.sha256_file(EXE)}\nrun_root={RUN_ROOT}\n"
        f"dt={DT}\nsegment_steps={SEGMENT_STEPS}\nmax_segments={MAX_SEGMENTS}\n"
        f"horizon_steps={horizon_steps}\nhorizon_time={horizon_time}\n"
        f"icd_source={CURRENT_SOURCE}\nprofile_width={PROFILE_WIDTH}\nshoulder_width={SHOULDER_WIDTH}\nshoulder_separation={SHOULDER_DELTA}\n"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

