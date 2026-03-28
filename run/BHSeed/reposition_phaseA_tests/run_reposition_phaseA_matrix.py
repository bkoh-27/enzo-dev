#!/usr/bin/env python3
import json
import math
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[3]  # enzo-dev
TS3_DIR = ROOT / "run/BHSeed/TS3_wrap"
BASE_RESTART = TS3_DIR / "DD0001"
ENZO_EXE = ROOT / "src/enzo/enzo.exe"
OUT_ROOT = Path(__file__).resolve().parent / f"matrix_{time.strftime('%Y%m%d_%H%M%S')}"

MODULE_INIT = "source /etc/profile.d/modules.sh && module load gcc openmpi hdf5/1.12.3"
ENV_EXPORTS = (
    "export PATH=/mnt/home/boh10/bh_proj/grackle/src/clib:/mnt/home/boh10/bh_proj/local/hdf5-1.10.11/bin:$PATH && "
    "export LD_LIBRARY_PATH=/mnt/home/boh10/ltu_proj/local/lib:/mnt/home/boh10/bh_proj/local/lib:/mnt/home/boh10/bh_proj/local/hdf5-1.10.11/lib:$LD_LIBRARY_PATH"
)

LENGTH_UNITS = 3.0857e23  # cm
BOX_KPC = LENGTH_UNITS / 3.0857e21  # 100 kpc
CELL_WIDTH = 1.0 / 50.0


@dataclass
class TestResult:
    name: str
    passed: bool
    details: str
    metrics: dict


def ensure_clean_dir(path: Path):
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def set_param(text: str, key: str, value: str) -> str:
    pat = re.compile(rf"^{re.escape(key)}\\s*=.*$", re.MULTILINE)
    line = f"{key:<30}= {value}"
    if pat.search(text):
        return pat.sub(line, text)
    return text + "\n" + line + "\n"


def run_bash(cmd: str):
    return subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)


def parse_kv_line(line: str, tag: str):
    out = {}
    for tok in line.strip().split():
        if "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        if k.startswith(f"[{tag}"):
            continue
        try:
            if re.match(r"^-?\d+$", v):
                out[k] = int(v)
            else:
                out[k] = float(v)
        except Exception:
            out[k] = v
    return out


def parse_tag_lines(log_path: Path, tag: str):
    raw = []
    parsed = []
    if not log_path.exists():
        return parsed, raw
    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if f"[{tag}]" in line:
                raw.append(line.strip())
                parsed.append(parse_kv_line(line, tag))
    return parsed, raw


def get_last_repos(log_path: Path, bh_id: int = None):
    parsed, raw = parse_tag_lines(log_path, "BHREPOS")
    if not parsed:
        raise RuntimeError(f"No [BHREPOS] line in {log_path}")
    if bh_id is None:
        return parsed[-1], raw[-1], parsed, raw
    chosen = None
    for p, r in zip(parsed, raw):
        if int(p.get("bh_id", -999999)) == bh_id:
            chosen = (p, r)
    if chosen is None:
        raise RuntimeError(f"No [BHREPOS] line with bh_id={bh_id} in {log_path}")
    return chosen[0], chosen[1], parsed, raw


def idx_from_center(pos, n=50):
    return int(round(pos * n - 0.5))


def center_from_idx(i, n=50):
    return (i + 0.5) / n


def configure_common_restart_param(text: str, search_radius_kpc: float, diagnose_potential: int, accretion_method: int = 0):
    text = set_param(text, "StopTime", "0.012")
    text = set_param(text, "BHSeedingMethod", "0")
    text = set_param(text, "BHSeedVerbose", "0")
    text = set_param(text, "BHAccretionMethod", str(accretion_method))
    text = set_param(text, "BHAccretionRunEveryTimestep", "1")
    text = set_param(text, "BHRepositionMethod", "0")
    text = set_param(text, "BHRepositionSearchRadius", f"{search_radius_kpc:.6g}")
    text = set_param(text, "BHRepositionMaxDisplacement", "0.5")
    text = set_param(text, "BHRepositionDiagnosePotential", str(diagnose_potential))
    text = set_param(text, "BHRepositionVerbose", "1")
    return text


def prepare_restart_case(name: str, search_radius_kpc: float = 6.0, diagnose_potential: int = 0, accretion_method: int = 0):
    case_dir = OUT_ROOT / name
    ensure_clean_dir(case_dir)
    shutil.copytree(BASE_RESTART, case_dir / "DD0001")

    pfile = case_dir / "DD0001" / "data0001"
    text = pfile.read_text()
    text = configure_common_restart_param(text, search_radius_kpc, diagnose_potential, accretion_method=accretion_method)
    pfile.write_text(text)
    return case_dir, case_dir / "DD0001" / "data0001.cpu0000"


def prepare_ts3_ic_case(name: str):
    case_dir = OUT_ROOT / name
    ensure_clean_dir(case_dir)
    for fname in ["bhseed_ts3wrap.enzo", "make_ts3wrap_ic.py", "assert_ts3wrap.sh"]:
        shutil.copy2(TS3_DIR / fname, case_dir / fname)

    pfile = case_dir / "bhseed_ts3wrap.enzo"
    text = pfile.read_text()
    text = set_param(text, "BHRepositionMethod", "0")
    text = set_param(text, "BHRepositionSearchRadius", "6.0")
    text = set_param(text, "BHRepositionMaxDisplacement", "0.5")
    text = set_param(text, "BHRepositionDiagnosePotential", "1")
    text = set_param(text, "BHRepositionVerbose", "1")
    text = set_param(text, "BHAccretionRunEveryTimestep", "1")
    pfile.write_text(text)
    return case_dir


def run_restart(case_dir: Path, np_ranks: int, log_name: str):
    cmd = (
        f"cd {case_dir} && "
        f"{MODULE_INIT} && "
        f"{ENV_EXPORTS} && "
        f"mpirun -np {np_ranks} {ENZO_EXE} -r DD0001/data0001 > {log_name} 2>&1"
    )
    rc = run_bash(cmd).returncode
    return rc, case_dir / log_name


def run_ts3_ic(case_dir: Path, np_ranks: int, log_name: str):
    cmd = (
        f"cd {case_dir} && "
        f"python3 make_ts3wrap_ic.py > ic.log 2>&1 && "
        f"{MODULE_INIT} && "
        f"{ENV_EXPORTS} && "
        f"mpirun -np {np_ranks} {ENZO_EXE} bhseed_ts3wrap.enzo > {log_name} 2>&1"
    )
    rc = run_bash(cmd).returncode
    return rc, case_dir / log_name


def set_particle_state(g, bh_pos):
    ptype = np.array(g["particle_type"][:])
    bh = np.where(ptype == 8)[0]
    if bh.size == 0:
        raise RuntimeError("No MBH particle in restart")
    i = int(bh[0])
    g["particle_position_x"][i] = bh_pos[0]
    g["particle_position_y"][i] = bh_pos[1]
    g["particle_position_z"][i] = bh_pos[2]
    g["particle_velocity_x"][i] = 0.0
    g["particle_velocity_y"][i] = 0.0
    g["particle_velocity_z"][i] = 0.0
    g["particle_index"][i] = 1


def set_uniform_density(g, rho):
    g["Density"][...] = rho
    g["GasEnergy"][...] = 1.0
    g["TotalEnergy"][...] = 1.0
    g["x-velocity"][...] = 0.0
    g["y-velocity"][...] = 0.0
    g["z-velocity"][...] = 0.0


def set_cell(g, i, j, k, rho):
    g["Density"][k, j, i] = rho


def almost(a, b, atol=1e-8):
    return abs(a - b) <= atol


def run_test_case_restart(name, configure_fn, check_fn, search_radius_kpc=6.0, diagnose_potential=0, accretion_method=0, np_ranks=1, log_name="run.log"):
    case_dir, h5path = prepare_restart_case(name, search_radius_kpc, diagnose_potential, accretion_method=accretion_method)
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        configure_fn(g)

    rc, log_path = run_restart(case_dir, np_ranks=np_ranks, log_name=log_name)
    if rc != 0:
        txt = log_path.read_text(errors="ignore")[-1200:] if log_path.exists() else ""
        return TestResult(name, False, f"run failed rc={rc}; tail={txt}", {"rc": rc, "log": str(log_path)})

    try:
        result = check_fn(case_dir, log_path)
    except Exception as e:
        return TestResult(name, False, f"check exception: {e}", {"log": str(log_path)})
    return result


def t1_config(g):
    bh = (0.15, 0.49, 0.49)
    set_particle_state(g, bh)
    set_uniform_density(g, 1e-6)
    i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
    set_cell(g, i0 + 2, j0, k0, 1e-3)


def t1_check(case_dir, log_path):
    p, raw, _, _ = get_last_repos(log_path, bh_id=1)
    expected = 2 * CELL_WIDTH * BOX_KPC
    ok = (
        int(p["newly_seeded_skip"]) == 0
        and abs(p["diag_peak_offset_kpc"] - expected) < 1e-6
        and abs(p["active_peak_offset_kpc"] - expected) < 1e-6
        and int(p["active_target_exists"]) == 1
    )
    return TestResult("test1_offset_measurement", ok, raw, {
        "diag_peak_offset_kpc": p["diag_peak_offset_kpc"],
        "active_peak_offset_kpc": p["active_peak_offset_kpc"],
        "expected_kpc": expected,
    })


def t2_config(g):
    bh = (0.15, 0.49, 0.49)
    set_particle_state(g, bh)
    set_uniform_density(g, 1e-6)
    i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
    set_cell(g, i0, j0, k0, 1e-3)


def t2_check(case_dir, log_path):
    p, raw, _, _ = get_last_repos(log_path, bh_id=1)
    ok = (
        abs(p["diag_peak_offset_kpc"]) < 1e-12
        and abs(p["active_peak_offset_kpc"]) < 1e-12
        and int(p["active_target_exists"]) == 0
    )
    return TestResult("test2_zero_offset", ok, raw, {
        "diag_peak_offset_kpc": p["diag_peak_offset_kpc"],
        "active_peak_offset_kpc": p["active_peak_offset_kpc"],
        "active_target_exists": int(p["active_target_exists"]),
    })


def t3_config(g):
    bh = (0.15, 0.49, 0.49)
    set_particle_state(g, bh)
    set_uniform_density(g, 1e-6)


def t3_check(case_dir, log_path):
    p, raw, _, _ = get_last_repos(log_path, bh_id=1)
    ok = (
        abs(p["diag_peak_offset_kpc"]) < 1e-12
        and abs(p["active_peak_offset_kpc"]) < 1e-12
        and int(p["active_target_exists"]) == 0
    )
    return TestResult("test3_uniform_no_target", ok, raw, {
        "diag_peak_offset_kpc": p["diag_peak_offset_kpc"],
        "active_peak_offset_kpc": p["active_peak_offset_kpc"],
        "active_target_exists": int(p["active_target_exists"]),
    })


def t4_config(g):
    bh = (0.15, 0.49, 0.49)
    set_particle_state(g, bh)
    set_uniform_density(g, 1e-6)
    i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
    # Two equal maxima at equal distance; lexicographic tie should pick lower x then y then z.
    # Candidate A: (i0+2, j0, k0), Candidate B: (i0, j0+2, k0).
    set_cell(g, i0 + 2, j0, k0, 1e-3)
    set_cell(g, i0, j0 + 2, k0, 1e-3)


def t4_check(case_dir, log_path):
    p, raw, _, _ = get_last_repos(log_path, bh_id=1)
    # Expect Candidate B because x is smaller (0.15 < 0.19)
    expx = center_from_idx(idx_from_center(0.15))
    expy = center_from_idx(idx_from_center(0.49) + 2)
    expz = center_from_idx(idx_from_center(0.49))
    ok = (
        abs(p["active_peak_pos_x"] - expx) < 1e-12
        and abs(p["active_peak_pos_y"] - expy) < 1e-12
        and abs(p["active_peak_pos_z"] - expz) < 1e-12
    )
    return TestResult("test4_tiebreak_determinism", ok, raw, {
        "active_peak_pos_x": p["active_peak_pos_x"],
        "active_peak_pos_y": p["active_peak_pos_y"],
        "active_peak_pos_z": p["active_peak_pos_z"],
        "expected_x": expx,
        "expected_y": expy,
        "expected_z": expz,
    })


def t5_config(g):
    bh = (0.15, 0.49, 0.49)
    set_particle_state(g, bh)
    set_uniform_density(g, 1e-6)


def t5_check(case_dir, log_path):
    p, raw, _, _ = get_last_repos(log_path, bh_id=1)
    logtxt = log_path.read_text(errors="ignore")
    has_warn = "Potential field not available for repositioning diagnostics" in logtxt
    ok = (abs(p["offset_from_potential_kpc"] + 1.0) < 1e-12 and has_warn)
    return TestResult("test5_potential_diagnostic", ok, raw, {
        "offset_from_potential_kpc": p["offset_from_potential_kpc"],
        "warning_present": int(has_warn),
    })


def t6_config(g):
    # BH near x-min boundary; plant a strong peak directly in a ghost cell
    # near the kernel to force diagnostic(active+ghost) != active-only peaks.
    bh = (0.01, 0.49, 0.49)
    set_particle_state(g, bh)
    set_uniform_density(g, 1e-6)
    j0, k0 = idx_from_center(0.49), idx_from_center(0.49)
    set_cell(g, 0, j0, k0, 1e-3)


def t6_check(case_dir, log_path):
    p, raw, _, _ = get_last_repos(log_path, bh_id=1)
    ok = (
        int(p["diag_peak_in_ghost"]) == 1
        and p["diag_peak_offset_kpc"] > 0.0
        and p["search_cells"] > p["search_active_cells"]
        and (p["active_peak_offset_kpc"] != p["diag_peak_offset_kpc"])
    )
    return TestResult("test6_ghost_vs_active_peak", ok, raw, {
        "diag_peak_in_ghost": int(p["diag_peak_in_ghost"]),
        "diag_peak_offset_kpc": p["diag_peak_offset_kpc"],
        "active_peak_offset_kpc": p["active_peak_offset_kpc"],
        "search_cells": int(p["search_cells"]),
        "search_active_cells": int(p["search_active_cells"]),
    })


def t7_config(g):
    # Place BH well inside the domain to avoid subgrid-boundary truncation in np=4.
    bh = (0.37, 0.37, 0.37)
    set_particle_state(g, bh)
    set_uniform_density(g, 1e-6)
    i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
    set_cell(g, i0 + 2, j0, k0, 1e-3)


def t7_run_and_check():
    name = "test7_mpi_determinism"
    case_dir, h5path = prepare_restart_case(name, search_radius_kpc=6.0, diagnose_potential=0, accretion_method=0)
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        t7_config(g)

    rc1, log1 = run_restart(case_dir, np_ranks=1, log_name="run_np1.log")
    rc4, log4 = run_restart(case_dir, np_ranks=4, log_name="run_np4.log")
    if rc1 != 0 or rc4 != 0:
        return TestResult(name, False, f"run failed rc1={rc1} rc4={rc4}", {"rc1": rc1, "rc4": rc4})

    p1, raw1, _, _ = get_last_repos(log1, bh_id=1)
    p4, raw4, _, _ = get_last_repos(log4, bh_id=1)

    skip = {"reposition_wall_ms"}
    keys = sorted(set(p1.keys()) & set(p4.keys()))
    mismatches = []
    for k in keys:
        if k in skip:
            continue
        if str(p1[k]) != str(p4[k]):
            mismatches.append((k, p1[k], p4[k]))

    ok = (len(mismatches) == 0)
    details = "bit-identical fields excluding reposition_wall_ms" if ok else f"mismatch sample: {mismatches[:5]}"
    return TestResult(name, ok, details, {
        "mismatch_count": len(mismatches),
        "np1_line": raw1,
        "np4_line": raw4,
    })


def t8_run_and_check():
    name = "test8_newly_seeded_exemption"
    case_dir = prepare_ts3_ic_case(name)
    rc, log = run_ts3_ic(case_dir, np_ranks=1, log_name="run.log")
    if rc != 0:
        tail = log.read_text(errors="ignore")[-1000:] if log.exists() else ""
        return TestResult(name, False, f"run failed rc={rc}; tail={tail}", {"rc": rc})

    parsed, raw = parse_tag_lines(log, "BHREPOS")
    if not parsed:
        return TestResult(name, False, "no [BHREPOS] line", {})

    chosen = None
    for p, r in zip(parsed, raw):
        if int(p.get("newly_seeded_skip", 0)) == 1:
            chosen = (p, r)
            break
    if chosen is None:
        return TestResult(name, False, "no newly_seeded_skip=1 line", {"n_lines": len(parsed)})

    p, rawline = chosen
    ok = (
        int(p["newly_seeded_skip"]) == 1
        and abs(p["diag_peak_offset_kpc"] + 1.0) < 1e-12
        and abs(p["active_peak_offset_kpc"] + 1.0) < 1e-12
        and abs(p["offset_from_potential_kpc"] + 1.0) < 1e-12
        and int(p["reposition_occurred"]) == 0
    )
    return TestResult(name, ok, rawline, {
        "newly_seeded_skip": int(p["newly_seeded_skip"]),
        "diag_peak_offset_kpc": p["diag_peak_offset_kpc"],
        "active_peak_offset_kpc": p["active_peak_offset_kpc"],
        "offset_from_potential_kpc": p["offset_from_potential_kpc"],
    })


def t9_run_and_check():
    name = "test9_ts3_regression"
    case_dir = prepare_ts3_ic_case(name)
    cmd = (
        f"cd {case_dir} && "
        f"python3 make_ts3wrap_ic.py > ic.log 2>&1 && "
        f"{MODULE_INIT} && "
        f"{ENV_EXPORTS} && "
        f"mpirun -np 1 {ENZO_EXE} bhseed_ts3wrap.enzo > ts3wrap_np1.log 2>&1 && "
        f"mpirun -np 4 {ENZO_EXE} bhseed_ts3wrap.enzo > ts3wrap_np4.log 2>&1 && "
        f"bash assert_ts3wrap.sh > assert.log 2>&1"
    )
    rc = run_bash(cmd).returncode
    assert_log = (case_dir / "assert.log").read_text(errors="ignore") if (case_dir / "assert.log").exists() else ""
    np1_log = case_dir / "ts3wrap_np1.log"
    np1_txt = np1_log.read_text(errors="ignore") if np1_log.exists() else ""

    ok = (
        rc == 0
        and "TS3_WRAP: PASS" in assert_log
        and "[BHREPOS]" in np1_txt
        and "[BHSEED]" in np1_txt
        and "[BHACCR]" in np1_txt
    )
    details = "TS3 wrap assert passed and BHREPOS/BHSEED/BHACCR present" if ok else f"rc={rc}, assert_tail={assert_log[-400:]}"
    return TestResult(name, ok, details, {
        "rc": rc,
        "assert_pass": int("TS3_WRAP: PASS" in assert_log),
        "has_BHREPOS": int("[BHREPOS]" in np1_txt),
        "has_BHSEED": int("[BHSEED]" in np1_txt),
        "has_BHACCR": int("[BHACCR]" in np1_txt),
    })


def baseline_test0_reference():
    p = TS3_DIR / "phaseA_test0_assert.log"
    if p.exists() and "TS3_WRAP: PASS" in p.read_text(errors="ignore"):
        return TestResult("test0_baseline_ts3_wrap", True, "phaseA_test0_assert.log contains TS3_WRAP: PASS", {"log": str(p)})
    return TestResult("test0_baseline_ts3_wrap", False, "missing baseline phaseA_test0_assert.log pass evidence", {"log": str(p)})


def run_validation_param_errors():
    name = "param_validation"
    case_dir, h5path = prepare_restart_case(name, search_radius_kpc=6.0, diagnose_potential=0, accretion_method=0)
    pfile = case_dir / "DD0001" / "data0001"
    txt = pfile.read_text()
    txt = set_param(txt, "BHRepositionSearchRadius", "0.0")
    pfile.write_text(txt)

    rc, log = run_restart(case_dir, np_ranks=1, log_name="run_invalid.log")
    logtxt = log.read_text(errors="ignore") if log.exists() else ""
    ok = (rc != 0 and "BHRepositionSearchRadius" in logtxt)
    return TestResult("param_validation_search_radius", ok, f"rc={rc}", {"rc": rc})


def main():
    ensure_clean_dir(OUT_ROOT)
    results = []

    results.append(baseline_test0_reference())

    results.append(run_test_case_restart("test1_offset_measurement", t1_config, t1_check, search_radius_kpc=6.0))
    results.append(run_test_case_restart("test2_zero_offset", t2_config, t2_check, search_radius_kpc=6.0))
    results.append(run_test_case_restart("test3_uniform_no_target", t3_config, t3_check, search_radius_kpc=6.0))
    results.append(run_test_case_restart("test4_tiebreak_determinism", t4_config, t4_check, search_radius_kpc=6.0))
    results.append(run_test_case_restart("test5_potential_diagnostic", t5_config, t5_check, search_radius_kpc=6.0, diagnose_potential=1))
    results.append(run_test_case_restart("test6_ghost_vs_active_peak", t6_config, t6_check, search_radius_kpc=6.0))

    results.append(t7_run_and_check())
    results.append(t8_run_and_check())
    results.append(t9_run_and_check())

    results.append(run_validation_param_errors())

    passed = [r for r in results if r.passed]
    failed = [r for r in results if not r.passed]

    report = {
        "out_root": str(OUT_ROOT),
        "pass_count": len(passed),
        "fail_count": len(failed),
        "results": [
            {
                "name": r.name,
                "passed": r.passed,
                "details": r.details,
                "metrics": r.metrics,
            }
            for r in results
        ],
    }

    with open(OUT_ROOT / "reposition_phaseA_results.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    with open(OUT_ROOT / "reposition_phaseA_summary.txt", "w", encoding="utf-8") as f:
        for r in results:
            f.write(f"{'PASS' if r.passed else 'FAIL'}  {r.name}: {r.details}\n")
        f.write(f"\nTotal: {len(passed)} passed, {len(failed)} failed\n")

    print(f"Results written to {OUT_ROOT}")
    for r in results:
        print(f"{'PASS' if r.passed else 'FAIL'}  {r.name}: {r.details}")

    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
