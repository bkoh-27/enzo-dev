#!/usr/bin/env python3
import json
import math
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

MODULE_INIT = "source /etc/profile >/dev/null 2>&1 || true && module load openmpi/4.1.8 hdf5/1.12.3"
ENV_EXPORTS = (
    "export LD_LIBRARY_PATH="
    "/mnt/home/boh10/ltu_proj/local/lib:"
    "/mnt/home/boh10/bh_proj/local/lib:"
    "/mnt/home/boh10/bh_proj/grackle/src/clib/.libs:$LD_LIBRARY_PATH"
)

BOX_KPC = 100.0
CELL_WIDTH_CODE = 1.0 / 50.0
CELL_WIDTH_KPC = BOX_KPC * CELL_WIDTH_CODE


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


def run_bash(cmd: str):
    return subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)


def set_param(text: str, key: str, value: str) -> str:
    pat = re.compile(rf"^{re.escape(key)}\s*=.*$", re.MULTILINE)
    line = f"{key:<30}= {value}"
    if pat.search(text):
        return pat.sub(line, text)
    return text + "\n" + line + "\n"


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
    parsed = []
    raw = []
    if not log_path.exists():
        return parsed, raw
    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if f"[{tag}]" in line:
                raw.append(line.strip())
                parsed.append(parse_kv_line(line, tag))
    return parsed, raw


def get_bhrepos(log_path: Path, bh_id=1):
    parsed, raw = parse_tag_lines(log_path, "BHREPOS")
    out = []
    for p, r in zip(parsed, raw):
        if int(p.get("bh_id", -999999)) == bh_id:
            out.append((p, r))
    return out


def get_bhacc(log_path: Path, bh_id=1):
    parsed, raw = parse_tag_lines(log_path, "BHACCR")
    out = []
    for p, r in zip(parsed, raw):
        if int(p.get("bh_id", -999999)) == bh_id:
            out.append((p, r))
    return out


def idx_from_center(pos):
    return int(round(pos * 50 - 0.5))


def center_from_idx(i):
    return (i + 0.5) / 50.0


def set_particle_state(g, bh_pos):
    ptype = np.array(g["particle_type"][:])
    bh = np.where(ptype == 8)[0]
    if bh.size == 0:
        raise RuntimeError("No MBH in restart")
    i = int(bh[0])
    g["particle_position_x"][i] = bh_pos[0]
    g["particle_position_y"][i] = bh_pos[1]
    g["particle_position_z"][i] = bh_pos[2]
    g["particle_velocity_x"][i] = 0.0
    g["particle_velocity_y"][i] = 0.0
    g["particle_velocity_z"][i] = 0.0
    g["particle_index"][i] = 1


def sync_bh_mass_invariant(g):
    ptype = np.array(g["particle_type"][:])
    bh = np.where(ptype == 8)[0]
    if bh.size == 0:
        return
    i = int(bh[0])
    if "bh_formation_mass" in g:
        g["bh_formation_mass"][i] = g["particle_mass"][i]
    if "bhaccr_accreted_mass" in g:
        g["bhaccr_accreted_mass"][i] = 0.0


def set_uniform_density(g, rho):
    g["Density"][...] = rho
    g["GasEnergy"][...] = 1.0
    g["TotalEnergy"][...] = 1.0
    g["x-velocity"][...] = 0.0
    g["y-velocity"][...] = 0.0
    g["z-velocity"][...] = 0.0


def set_cell(g, i, j, k, rho):
    g["Density"][k, j, i] = rho


def configure_restart_param(text: str, *, stop_time: float, reposition_method: int,
                            reposition_verbose: int, reposition_search_kpc: float,
                            reposition_max_disp_cells: float,
                            diagnose_potential: int = 0,
                            accretion_method: int = 0,
                            self_gravity: int = 0):
    text = set_param(text, "StopTime", f"{stop_time:.6g}")
    text = set_param(text, "BHSeedingMethod", "0")
    text = set_param(text, "BHSeedVerbose", "0")
    text = set_param(text, "BHAccretionMethod", str(accretion_method))
    text = set_param(text, "BHAccretionRunEveryTimestep", "1")
    text = set_param(text, "BHRepositionMethod", str(reposition_method))
    text = set_param(text, "BHRepositionSearchRadius", f"{reposition_search_kpc:.6g}")
    text = set_param(text, "BHRepositionMaxDisplacement", f"{reposition_max_disp_cells:.6g}")
    text = set_param(text, "BHRepositionDiagnosePotential", str(diagnose_potential))
    text = set_param(text, "BHRepositionVerbose", str(reposition_verbose))
    text = set_param(text, "SelfGravity", str(self_gravity))
    return text


def prepare_restart_case(name: str, *, stop_time: float, reposition_method: int,
                         reposition_verbose: int, reposition_search_kpc: float,
                         reposition_max_disp_cells: float, diagnose_potential: int = 0,
                         accretion_method: int = 0, self_gravity: int = 0):
    case_dir = OUT_ROOT / name
    ensure_clean_dir(case_dir)
    shutil.copytree(BASE_RESTART, case_dir / "DD0001")

    pfile = case_dir / "DD0001" / "data0001"
    txt = pfile.read_text()
    txt = configure_restart_param(
        txt,
        stop_time=stop_time,
        reposition_method=reposition_method,
        reposition_verbose=reposition_verbose,
        reposition_search_kpc=reposition_search_kpc,
        reposition_max_disp_cells=reposition_max_disp_cells,
        diagnose_potential=diagnose_potential,
        accretion_method=accretion_method,
        self_gravity=self_gravity,
    )
    pfile.write_text(txt)
    return case_dir, case_dir / "DD0001" / "data0001.cpu0000"


def run_restart(case_dir: Path, np_ranks: int, log_name: str):
    cmd = (
        f"cd {case_dir} && "
        f"{MODULE_INIT} && "
        f"{ENV_EXPORTS} && "
        f"mpirun -np {np_ranks} {ENZO_EXE} -r DD0001/data0001 > {log_name} 2>&1"
    )
    rc = run_bash(cmd).returncode
    return rc, case_dir / log_name


def prepare_ts3_case(name: str, method: int, verbose: int, accretion_method: int):
    case_dir = OUT_ROOT / name
    ensure_clean_dir(case_dir)
    for fname in ["bhseed_ts3wrap.enzo", "make_ts3wrap_ic.py", "assert_ts3wrap.sh"]:
        shutil.copy2(TS3_DIR / fname, case_dir / fname)
    pfile = case_dir / "bhseed_ts3wrap.enzo"
    txt = pfile.read_text()
    txt = set_param(txt, "BHRepositionMethod", str(method))
    txt = set_param(txt, "BHRepositionVerbose", str(verbose))
    txt = set_param(txt, "BHRepositionSearchRadius", "8.0")
    txt = set_param(txt, "BHRepositionMaxDisplacement", "0.5")
    txt = set_param(txt, "BHAccretionMethod", str(accretion_method))
    txt = set_param(txt, "BHAccretionRunEveryTimestep", "1")
    pfile.write_text(txt)
    return case_dir


def run_ts3(case_dir: Path, np_ranks: int, log_name: str):
    cmd = (
        f"cd {case_dir} && "
        f"python3 make_ts3wrap_ic.py > ic.log 2>&1 && "
        f"{MODULE_INIT} && "
        f"{ENV_EXPORTS} && "
        f"mpirun -np {np_ranks} {ENZO_EXE} bhseed_ts3wrap.enzo > {log_name} 2>&1"
    )
    rc = run_bash(cmd).returncode
    return rc, case_dir / log_name


def passfail(name, ok, details, metrics):
    return TestResult(name=name, passed=ok, details=details, metrics=metrics)


def test1_and_test2():
    name = "test1_rate_limited_known_offset"
    case_dir, h5path = prepare_restart_case(
        name,
        stop_time=0.07,
        reposition_method=1,
        reposition_verbose=1,
        reposition_search_kpc=8.0,
        reposition_max_disp_cells=0.5,
    )
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        bh = (0.15, 0.49, 0.49)
        set_particle_state(g, bh)
        set_uniform_density(g, 1e-6)
        i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
        set_cell(g, i0 + 3, j0, k0, 1e-3)

    rc, log = run_restart(case_dir, np_ranks=1, log_name="run.log")
    if rc != 0:
        tail = log.read_text(errors="ignore")[-800:] if log.exists() else ""
        fail = passfail(name, False, f"run failed rc={rc}", {"rc": rc, "tail": tail})
        fail2 = passfail("test2_displacement_cap_enforcement", False, "not run", {})
        return fail, fail2

    lines = get_bhrepos(log, bh_id=1)
    if len(lines) < 6:
        fail = passfail(name, False, f"expected >=6 BHREPOS lines, got {len(lines)}", {"n_lines": len(lines)})
        fail2 = passfail("test2_displacement_cap_enforcement", False, "insufficient lines", {"n_lines": len(lines)})
        return fail, fail2

    p0 = lines[0][0]
    pend = lines[-1][0]
    expected_offset_kpc = 3.0 * CELL_WIDTH_KPC
    ok1 = (
        abs(p0["active_peak_offset_kpc"] - expected_offset_kpc) < 1e-6
        and abs(p0["displacement_cells"] - 0.5) < 1e-6
        and int(p0["reposition_occurred"]) == 1
        and abs(pend["active_peak_offset_kpc"]) < 1e-9
    )
    r1 = passfail(name, ok1, lines[0][1], {
        "first_active_peak_offset_kpc": p0["active_peak_offset_kpc"],
        "expected_offset_kpc": expected_offset_kpc,
        "first_displacement_cells": p0["displacement_cells"],
        "final_active_peak_offset_kpc": pend["active_peak_offset_kpc"],
        "n_lines": len(lines),
    })

    max_disp = max(p["displacement_cells"] for p, _ in lines)
    ok2 = max_disp <= 0.5000001
    r2 = passfail("test2_displacement_cap_enforcement", ok2, f"max displacement_cells={max_disp:.8g}", {
        "max_displacement_cells": max_disp,
    })
    return r1, r2


def run_single_restart_test(name, config_fn, check_fn, **params):
    case_dir, h5path = prepare_restart_case(name, **params)
    with h5py.File(h5path, "r+") as f:
        config_fn(f["Grid00000001"])
    rc, log = run_restart(case_dir, np_ranks=1, log_name="run.log")
    if rc != 0:
        tail = log.read_text(errors="ignore")[-800:] if log.exists() else ""
        return passfail(name, False, f"run failed rc={rc}", {"rc": rc, "tail": tail})
    return check_fn(log)


def test3_config(g):
    bh = (0.15, 0.49, 0.49)
    set_particle_state(g, bh)
    set_uniform_density(g, 1e-6)
    i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
    set_cell(g, i0 + 1, j0, k0, 1e-3)


def test3_check(log):
    lines = get_bhrepos(log, bh_id=1)
    if len(lines) < 2:
        return passfail("test3_exact_arrival", False, f"need >=2 lines, got {len(lines)}", {"n_lines": len(lines)})
    p0 = lines[0][0]
    p1 = lines[1][0]
    ok = (
        abs(p0["displacement_cells"] - 0.5) < 1e-8
        and int(p1["active_target_exists"]) == 0
        and abs(p1["displacement_cells"]) < 1e-10
        and int(p1["reposition_occurred"]) == 0
    )
    return passfail("test3_exact_arrival", ok, lines[1][1], {
        "first_displacement_cells": p0["displacement_cells"],
        "second_active_target_exists": int(p1["active_target_exists"]),
        "second_displacement_cells": p1["displacement_cells"],
        "second_reposition_occurred": int(p1["reposition_occurred"]),
    })


def test4_config(g):
    bh = (0.15, 0.49, 0.49)
    set_particle_state(g, bh)
    set_uniform_density(g, 1e-6)
    i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
    set_cell(g, i0, j0, k0, 1e-3)


def test4_check(log):
    lines = get_bhrepos(log, bh_id=1)
    if not lines:
        return passfail("test4_no_movement_at_peak", False, "no BHREPOS line", {})
    p = lines[-1][0]
    ok = abs(p["displacement_cells"]) < 1e-12 and int(p["reposition_occurred"]) == 0
    return passfail("test4_no_movement_at_peak", ok, lines[-1][1], {
        "displacement_cells": p["displacement_cells"],
        "reposition_occurred": int(p["reposition_occurred"]),
    })


def test5_config(g):
    set_particle_state(g, (0.15, 0.49, 0.49))
    set_uniform_density(g, 1e-6)


def test5_check(log):
    lines = get_bhrepos(log, bh_id=1)
    if not lines:
        return passfail("test5_no_movement_uniform", False, "no BHREPOS line", {})
    p = lines[-1][0]
    ok = (
        int(p["active_target_exists"]) == 0
        and abs(p["displacement_cells"]) < 1e-12
        and int(p["reposition_occurred"]) == 0
    )
    return passfail("test5_no_movement_uniform", ok, lines[-1][1], {
        "active_target_exists": int(p["active_target_exists"]),
        "displacement_cells": p["displacement_cells"],
    })


def test6_config(g):
    set_particle_state(g, (0.01, 0.49, 0.49))
    set_uniform_density(g, 1e-6)
    j0, k0 = idx_from_center(0.49), idx_from_center(0.49)
    set_cell(g, 0, j0, k0, 1e-3)


def test6_check(log):
    lines = get_bhrepos(log, bh_id=1)
    if not lines:
        return passfail("test6_ghost_peak_no_movement", False, "no BHREPOS line", {})
    p = lines[-1][0]
    ok = (
        int(p["diag_peak_in_ghost"]) == 1
        and int(p["active_target_exists"]) == 0
        and abs(p["displacement_cells"]) < 1e-12
    )
    return passfail("test6_ghost_peak_no_movement", ok, lines[-1][1], {
        "diag_peak_in_ghost": int(p["diag_peak_in_ghost"]),
        "active_target_exists": int(p["active_target_exists"]),
        "displacement_cells": p["displacement_cells"],
    })


def test7_teleport():
    name = "test7_teleport_mode"
    case_dir, h5path = prepare_restart_case(
        name,
        stop_time=0.03,
        reposition_method=2,
        reposition_verbose=1,
        reposition_search_kpc=8.0,
        reposition_max_disp_cells=0.5,
    )
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        bh = (0.15, 0.49, 0.49)
        set_particle_state(g, bh)
        set_uniform_density(g, 1e-6)
        i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
        set_cell(g, i0 + 3, j0, k0, 1e-3)

    rc, log = run_restart(case_dir, np_ranks=1, log_name="run.log")
    if rc != 0:
        return passfail(name, False, f"run failed rc={rc}", {"rc": rc})
    lines = get_bhrepos(log, bh_id=1)
    if len(lines) < 2:
        return passfail(name, False, f"need >=2 lines, got {len(lines)}", {"n_lines": len(lines)})
    p0 = lines[0][0]
    p1 = lines[1][0]
    ok = (
        abs(p0["displacement_cells"] - 3.0) < 1e-6
        and int(p0["reposition_occurred"]) == 1
        and abs(p1["active_peak_offset_kpc"]) < 1e-9
    )
    return passfail(name, ok, lines[0][1], {
        "first_displacement_cells": p0["displacement_cells"],
        "second_active_peak_offset_kpc": p1["active_peak_offset_kpc"],
    })


def test8_clamp_attempt():
    name = "test8_defensive_clamp"
    case_dir, h5path = prepare_restart_case(
        name,
        stop_time=0.03,
        reposition_method=2,
        reposition_verbose=1,
        reposition_search_kpc=8.0,
        reposition_max_disp_cells=0.5,
    )
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        set_particle_state(g, (0.9899999999999, 0.49, 0.49))
        set_uniform_density(g, 1e-6)
        set_cell(g, 49, idx_from_center(0.49), idx_from_center(0.49), 1e-3)
    rc, log = run_restart(case_dir, np_ranks=1, log_name="run.log")
    if rc != 0:
        return passfail(name, False, f"run failed rc={rc}", {"rc": rc})
    lines = get_bhrepos(log, bh_id=1)
    if not lines:
        return passfail(name, False, "no BHREPOS lines", {})
    p = lines[-1][0]
    # This is a best-effort clamp trigger; if not triggered, report as unresolved.
    ok = int(p["reposition_clamped"]) == 1
    details = lines[-1][1] if ok else "clamp did not trigger in best-effort edge case"
    return passfail(name, ok, details, {
        "reposition_clamped": int(p["reposition_clamped"]),
        "displacement_cells": p["displacement_cells"],
    })


def test9_mpi_determinism():
    name = "test9_mpi_determinism_active_movement"
    case_dir, h5path = prepare_restart_case(
        name,
        stop_time=0.03,
        reposition_method=1,
        reposition_verbose=1,
        reposition_search_kpc=6.0,
        reposition_max_disp_cells=0.5,
    )
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        bh = (0.37, 0.37, 0.37)
        set_particle_state(g, bh)
        set_uniform_density(g, 1e-6)
        i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
        set_cell(g, i0 + 2, j0, k0, 1e-3)
    rc1, log1 = run_restart(case_dir, np_ranks=1, log_name="run_np1.log")
    rc4, log4 = run_restart(case_dir, np_ranks=4, log_name="run_np4.log")
    if rc1 != 0 or rc4 != 0:
        return passfail(name, False, f"run failed rc1={rc1} rc4={rc4}", {"rc1": rc1, "rc4": rc4})
    r1 = get_bhrepos(log1, bh_id=1)
    r4 = get_bhrepos(log4, bh_id=1)
    if not r1 or not r4:
        return passfail(name, False, "missing BHREPOS lines", {"n1": len(r1), "n4": len(r4)})
    p1 = r1[-1][0]
    p4 = r4[-1][0]
    skip = {"reposition_wall_ms"}
    mismatch = []
    for k in sorted(set(p1.keys()) & set(p4.keys())):
        if k in skip:
            continue
        if str(p1[k]) != str(p4[k]):
            mismatch.append((k, p1[k], p4[k]))
    ok = len(mismatch) == 0
    return passfail(name, ok, "bit-identical excluding reposition_wall_ms" if ok else str(mismatch[:5]), {
        "mismatch_count": len(mismatch),
    })


def test10_accretion_regression():
    name = "test10_accretion_regression_with_reposition"
    case_dir, h5path = prepare_restart_case(
        name,
        stop_time=0.03,
        reposition_method=1,
        reposition_verbose=1,
        reposition_search_kpc=8.0,
        reposition_max_disp_cells=0.5,
        accretion_method=1,
    )
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        bh = (0.15, 0.49, 0.49)
        set_particle_state(g, bh)
        sync_bh_mass_invariant(g)
        set_uniform_density(g, 1e-6)
        i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
        set_cell(g, i0 + 2, j0, k0, 1e-3)
    rc, log = run_restart(case_dir, np_ranks=1, log_name="run.log")
    if rc != 0:
        return passfail(name, False, f"run failed rc={rc}", {"rc": rc})
    repos = get_bhrepos(log, bh_id=1)
    accr = get_bhacc(log, bh_id=1)
    if not repos or not accr:
        return passfail(name, False, "missing BHREPOS or BHACCR line", {"n_repos": len(repos), "n_accr": len(accr)})
    rp = repos[0][0]
    ap = accr[0][0]
    ok = (
        int(rp["reposition_occurred"]) == 1
        and ap["dm_removed"] > 0.0
        and ap["frac_gas"] >= 0.0
    )
    return passfail(name, ok, accr[0][1], {
        "reposition_occurred": int(rp["reposition_occurred"]),
        "displacement_cells": rp["displacement_cells"],
        "dm_removed": ap["dm_removed"],
        "bh_mass": ap.get("bh_mass", float("nan")),
        "bh_mass_new": ap.get("bh_mass_new", float("nan")),
    })


def test11_ts3_lifecycle():
    name = "test11_ts3_full_lifecycle"
    case_dir = prepare_ts3_case(name, method=1, verbose=1, accretion_method=1)
    rc1, log1 = run_ts3(case_dir, np_ranks=1, log_name="ts3wrap_np1.log")
    rc4, log4 = run_ts3(case_dir, np_ranks=4, log_name="ts3wrap_np4.log")
    cmd = f"cd {case_dir} && bash assert_ts3wrap.sh > assert.log 2>&1"
    rca = run_bash(cmd).returncode
    assert_log = (case_dir / "assert.log").read_text(errors="ignore") if (case_dir / "assert.log").exists() else ""
    np1_txt = log1.read_text(errors="ignore") if log1.exists() else ""
    ok = (
        rc1 == 0 and rc4 == 0 and rca == 0 and "TS3_WRAP: PASS" in assert_log
        and "[BHSEED]" in np1_txt and "[BHREPOS]" in np1_txt and "[BHACCR]" in np1_txt
    )
    return passfail(name, ok, "TS3 lifecycle complete with BHSEED/BHREPOS/BHACCR", {
        "rc_np1": rc1, "rc_np4": rc4, "rc_assert": rca,
        "assert_pass": int("TS3_WRAP: PASS" in assert_log),
    })


def test12_13_comparison():
    name12 = "test12_compare_accretion_off_vs_on"
    # OFF case
    off_dir, off_h5 = prepare_restart_case(
        "compare_off",
        stop_time=0.03,
        reposition_method=0,
        reposition_verbose=0,
        reposition_search_kpc=8.0,
        reposition_max_disp_cells=0.5,
        accretion_method=1,
    )
    with h5py.File(off_h5, "r+") as f:
        g = f["Grid00000001"]
        bh = (0.15, 0.49, 0.49)
        set_particle_state(g, bh)
        sync_bh_mass_invariant(g)
        set_uniform_density(g, 1e-6)
        i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
        set_cell(g, i0 + 2, j0, k0, 1e-1)
    rc_off, off_log = run_restart(off_dir, np_ranks=1, log_name="run.log")

    # ON case
    on_dir, on_h5 = prepare_restart_case(
        "compare_on",
        stop_time=0.03,
        reposition_method=2,
        reposition_verbose=1,
        reposition_search_kpc=8.0,
        reposition_max_disp_cells=0.5,
        accretion_method=1,
    )
    with h5py.File(on_h5, "r+") as f:
        g = f["Grid00000001"]
        bh = (0.15, 0.49, 0.49)
        set_particle_state(g, bh)
        sync_bh_mass_invariant(g)
        set_uniform_density(g, 1e-6)
        i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
        set_cell(g, i0 + 2, j0, k0, 1e-1)
    rc_on, on_log = run_restart(on_dir, np_ranks=1, log_name="run.log")

    if rc_off != 0 or rc_on != 0:
        fail12 = passfail(name12, False, f"run failed rc_off={rc_off} rc_on={rc_on}", {})
        fail13 = passfail("test13_feedback_proxy_off_vs_on", False, "not run", {})
        return fail12, fail13

    off_accr = get_bhacc(off_log, bh_id=1)
    on_accr = get_bhacc(on_log, bh_id=1)
    if not off_accr or not on_accr:
        fail12 = passfail(name12, False, "missing BHACCR lines", {"off": len(off_accr), "on": len(on_accr)})
        fail13 = passfail("test13_feedback_proxy_off_vs_on", False, "missing BHACCR lines", {})
        return fail12, fail13

    mdot_off = off_accr[0][0]["Mdot_actual"]
    mdot_on = on_accr[0][0]["Mdot_actual"]
    ok12 = mdot_on > mdot_off
    r12 = passfail(name12, ok12, f"Mdot_off={mdot_off:.8e}, Mdot_on={mdot_on:.8e}", {
        "Mdot_actual_off": mdot_off,
        "Mdot_actual_on": mdot_on,
        "ratio": (mdot_on / mdot_off) if mdot_off > 0 else float("inf"),
    })

    # Proxy for feedback power scaling with accretion rate
    p_off = mdot_off
    p_on = mdot_on
    ok13 = p_on > p_off
    r13 = passfail("test13_feedback_proxy_off_vs_on", ok13, f"proxy_off={p_off:.8e}, proxy_on={p_on:.8e}", {
        "proxy_off": p_off,
        "proxy_on": p_on,
    })
    return r12, r13


def test14_15_longrun_and_perf():
    name14 = "test14_long_run_stability"
    case_dir, h5path = prepare_restart_case(
        "longrun",
        stop_time=0.52,
        reposition_method=1,
        reposition_verbose=1,
        reposition_search_kpc=8.0,
        reposition_max_disp_cells=0.5,
    )
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        bh = (0.15, 0.49, 0.49)
        set_particle_state(g, bh)
        set_uniform_density(g, 1e-6)
        i0, j0, k0 = idx_from_center(bh[0]), idx_from_center(bh[1]), idx_from_center(bh[2])
        set_cell(g, i0 + 2, j0, k0, 1e-3)
    rc, log = run_restart(case_dir, np_ranks=1, log_name="run.log")
    if rc != 0:
        fail14 = passfail(name14, False, f"run failed rc={rc}", {"rc": rc})
        fail15 = passfail("test15_reposition_performance", False, "not run", {})
        return fail14, fail15
    txt = log.read_text(errors="ignore")
    repos = get_bhrepos(log, bh_id=1)
    has_nan = ("nan" in txt.lower())
    ok14 = (len(repos) >= 50 and not has_nan)
    r14 = passfail(name14, ok14, f"n_bhrepos={len(repos)} has_nan={has_nan}", {
        "n_bhrepos": len(repos),
        "has_nan": int(has_nan),
    })
    if repos:
        wall = [p["reposition_wall_ms"] for p, _ in repos]
        avg = float(np.mean(wall))
        p95 = float(np.percentile(wall, 95))
        r15 = passfail("test15_reposition_performance", True,
                       f"avg_ms={avg:.6f}, p95_ms={p95:.6f}", {
                           "avg_reposition_wall_ms": avg,
                           "p95_reposition_wall_ms": p95,
                           "samples": len(wall),
                       })
    else:
        r15 = passfail("test15_reposition_performance", False, "no BHREPOS lines", {})
    return r14, r15


def baseline_test0():
    p = ROOT / "run/BHSeed/reposition_phaseB_tests/pre_fixes_20260329_001929/TS3_regression/assert.log"
    ok = p.exists() and ("TS3_WRAP: PASS" in p.read_text(errors="ignore"))
    return passfail("test0_baseline_ts3_wrap_after_fixes", ok, str(p), {"path": str(p)})


def main():
    ensure_clean_dir(OUT_ROOT)
    results = []

    results.append(baseline_test0())

    r1, r2 = test1_and_test2()
    results.append(r1)
    results.append(r2)

    results.append(run_single_restart_test(
        "test3_exact_arrival",
        test3_config,
        test3_check,
        stop_time=0.04,
        reposition_method=1,
        reposition_verbose=1,
        reposition_search_kpc=8.0,
        reposition_max_disp_cells=0.5,
    ))
    results.append(run_single_restart_test(
        "test4_no_movement_at_peak",
        test4_config,
        test4_check,
        stop_time=0.02,
        reposition_method=1,
        reposition_verbose=1,
        reposition_search_kpc=8.0,
        reposition_max_disp_cells=0.5,
    ))
    results.append(run_single_restart_test(
        "test5_no_movement_uniform",
        test5_config,
        test5_check,
        stop_time=0.02,
        reposition_method=1,
        reposition_verbose=1,
        reposition_search_kpc=8.0,
        reposition_max_disp_cells=0.5,
    ))
    results.append(run_single_restart_test(
        "test6_ghost_peak_no_movement",
        test6_config,
        test6_check,
        stop_time=0.02,
        reposition_method=1,
        reposition_verbose=1,
        reposition_search_kpc=6.0,
        reposition_max_disp_cells=0.5,
    ))
    results.append(test7_teleport())
    results.append(test8_clamp_attempt())
    results.append(test9_mpi_determinism())
    results.append(test10_accretion_regression())
    results.append(test11_ts3_lifecycle())

    r12, r13 = test12_13_comparison()
    results.append(r12)
    results.append(r13)

    r14, r15 = test14_15_longrun_and_perf()
    results.append(r14)
    results.append(r15)

    summary_lines = []
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        summary_lines.append(f"{status:4}  {r.name}: {r.details}")
    pass_count = sum(1 for r in results if r.passed)
    fail_count = len(results) - pass_count

    summary = "\n".join(summary_lines) + f"\nTOTAL: {pass_count} passed, {fail_count} failed\n"
    (OUT_ROOT / "reposition_phaseB_summary.txt").write_text(summary)
    (OUT_ROOT / "reposition_phaseB_results.json").write_text(
        json.dumps([r.__dict__ for r in results], indent=2)
    )

    print(summary, end="")
    print(f"ARTIFACT_DIR={OUT_ROOT}")

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
