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

# -----------------------------------------------------------------------------
# Paths and constants
# -----------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[3]  # enzo-dev
BASE_RESTART = ROOT / "run/BHSeed/TS3_wrap/DD0001"
TS3_DIR = ROOT / "run/BHSeed/TS3_wrap"
OUT_ROOT = ROOT / "run/BHSeed/phaseA_tests/matrix_20260324"
ENZO_EXE = ROOT / "src/enzo/enzo.exe"
COOL_RATES = ROOT / "input/cool_rates.in"

MODULE_INIT = "source /etc/profile.d/modules.sh && module load openmpi/4.1.8"
LD_PATH = ":".join([
    "/mnt/sw/nix/store/8r0bpn8x0ljxndwzf7jjkp15ir96qpfg-hdf5-1.10.11/lib",
    "/mnt/sw/nix/store/cprcg9nf8xpr7gljz5d3vgg97blyzwfx-hdf5-1.12.1/lib",
    "/mnt/home/boh10/bh_proj/local/lib",
    "/mnt/home/boh10/ltu_proj/local/lib",
])

# Unit / physics constants (must match test setup)
LENGTH_UNITS = 3.0857e23
TIME_UNITS = 3.1557e13
DENSITY_UNITS = 1.67e-24
VEL_UNITS = LENGTH_UNITS / TIME_UNITS
CELL_WIDTH = 1.0 / 50.0
BOX_KPC = LENGTH_UNITS / 3.0857e21
KERNEL_RADIUS_PHYS_KPC = 3.0
KERNEL_RADIUS_CODE = KERNEL_RADIUS_PHYS_KPC / BOX_KPC
GAMMA = 5.0 / 3.0
MU = 0.6
KBOLTZ = 1.3807e-16
MH = 1.6726e-24
XH = 0.76
G_CGS = 6.67259e-8
SIGMA_T = 6.6524587158e-25
CLIGHT = 2.99792458e10
MSUN = 1.989e33
YR_S = 3.1557e7
BHACCR_CVISC = 6.283

REL_TOL = 3e-3
ABS_TOL = 1e-12


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
    pattern = re.compile(rf"^{re.escape(key)}\\s*=.*$", re.MULTILINE)
    new_line = f"{key:<30}= {value}"
    if pattern.search(text):
        return pattern.sub(new_line, text)
    return text + "\n" + new_line + "\n"


def energy_code_from_temp(T_kelvin):
    return (KBOLTZ * T_kelvin / ((GAMMA - 1.0) * MU * MH)) / (VEL_UNITS ** 2)


def parse_bhaccr_line(line: str):
    # Parse key=value tokens from [BHACCR] line
    data = {}
    for tok in line.strip().split():
        if "=" not in tok:
            continue
        key, val = tok.split("=", 1)
        if key.startswith("[BHACCR"):
            continue
        try:
            if val.lower().startswith("nan"):
                data[key] = float("nan")
            elif val.lower().startswith("inf"):
                data[key] = float("inf")
            elif re.match(r"^-?\d+$", val):
                data[key] = int(val)
            else:
                data[key] = float(val)
        except Exception:
            data[key] = val
    return data


def get_last_bhaccr(log_path: Path, bh_id: int = 1):
    lines = []
    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "[BHACCR]" in line:
                lines.append(line.strip())
    if not lines:
        raise RuntimeError(f"No [BHACCR] line found in {log_path}")

    parsed = [parse_bhaccr_line(l) for l in lines]
    chosen = None
    for d, l in zip(parsed, lines):
        if int(d.get("bh_id", -99999)) == bh_id:
            chosen = (d, l)
    if chosen is None:
        chosen = (parsed[-1], lines[-1])
    return chosen[0], lines


def almost_equal(a, b, rel=REL_TOL, abs_tol=ABS_TOL):
    return abs(a - b) <= max(abs_tol, rel * max(abs(a), abs(b), 1.0))


def prepare_case(case_name: str, radiative_cooling: bool):
    case_dir = OUT_ROOT / case_name
    dd = case_dir / "DD0001"
    ensure_clean_dir(case_dir)
    shutil.copytree(BASE_RESTART, dd)

    # Parameter tuning
    pfile = dd / "data0001"
    text = pfile.read_text()
    # Continue for one additional step beyond InitialTime.
    # Keep native checkpoint flag untouched to avoid OldTime read issues
    # in copied multi-rank restart shards.
    text = set_param(text, "StopTime", "0.02")

    # Disable seeding for pure diagnostics tests
    text = set_param(text, "BHSeedingMethod", "0")
    text = set_param(text, "BHSeedVerbose", "0")

    # Accretion diagnostics settings
    text = set_param(text, "BHAccretionMethod", "1")
    text = set_param(text, "BHAccretionKernelRadius", "3.0")
    text = set_param(text, "BHAccretionTSplitFloor", "5.0e5")
    text = set_param(text, "BHAccretionColdModel", "0")
    text = set_param(text, "BHAccretionCVisc", "6.283")
    text = set_param(text, "BHAccretionNHStar", "0.1")
    text = set_param(text, "BHAccretionBeta", "1.0")
    text = set_param(text, "BHAccretionAlphaMax", "10.0")
    text = set_param(text, "BHAccretionRadiativeEfficiency", "0.1")
    text = set_param(text, "BHAccretionVerbose", "1")
    text = set_param(text, "BHAccretionRunEveryTimestep", "1")

    text = set_param(text, "RadiativeCooling", "1" if radiative_cooling else "0")
    pfile.write_text(text)

    if radiative_cooling:
        shutil.copy2(COOL_RATES, case_dir / "cool_rates.in")

    return case_dir, dd / "data0001.cpu0000"


def set_particle_state(g, bh_pos=(0.15, 0.5, 0.5), bh_vel=(0.0, 0.0, 0.0), bh_mass_code=None):
    ptype = np.array(g["particle_type"][:])
    bh_idx = np.where(ptype == 8)[0]
    if bh_idx.size == 0:
        raise RuntimeError("No MBH particle in restart file")
    i = int(bh_idx[0])

    g["particle_position_x"][i] = bh_pos[0]
    g["particle_position_y"][i] = bh_pos[1]
    g["particle_position_z"][i] = bh_pos[2]
    g["particle_velocity_x"][i] = bh_vel[0]
    g["particle_velocity_y"][i] = bh_vel[1]
    g["particle_velocity_z"][i] = bh_vel[2]
    if bh_mass_code is not None:
        g["particle_mass"][i] = bh_mass_code

    # Keep BH id deterministic for parsing
    g["particle_index"][i] = 1


def write_fields_uniform(g, rho_code, temp_k, vel=(0.0, 0.0, 0.0), metallicity_frac=1e-6):
    shape = g["Density"].shape
    vx, vy, vz = vel
    ge = energy_code_from_temp(temp_k)
    te = ge + 0.5 * (vx * vx + vy * vy + vz * vz)

    g["Density"][...] = rho_code
    g["GasEnergy"][...] = ge
    g["TotalEnergy"][...] = te
    g["x-velocity"][...] = vx
    g["y-velocity"][...] = vy
    g["z-velocity"][...] = vz
    g["Metal_Density"][...] = rho_code * metallicity_frac


def write_fields_mixed_temperature(g, rho_code, t_cold, t_hot, split_x=0.5, vel=(0, 0, 0), metallicity_frac=1e-6):
    nz, ny, nx = g["Density"].shape
    x_centers = (np.arange(nx) + 0.5) / nx
    cold_mask_x = x_centers < split_x
    temp = np.where(cold_mask_x[None, None, :], t_cold, t_hot)

    vx, vy, vz = vel
    ge = energy_code_from_temp(temp)
    te = ge + 0.5 * (vx * vx + vy * vy + vz * vz)

    g["Density"][...] = rho_code
    g["GasEnergy"][...] = ge
    g["TotalEnergy"][...] = te
    g["x-velocity"][...] = vx
    g["y-velocity"][...] = vy
    g["z-velocity"][...] = vz
    g["Metal_Density"][...] = rho_code * metallicity_frac


def write_fields_swirl_cold(g, rho_code, temp_k, bh_pos, omega_code, metallicity_frac=1e-6):
    nz, ny, nx = g["Density"].shape
    x = (np.arange(nx) + 0.5) / nx
    y = (np.arange(ny) + 0.5) / ny
    z = (np.arange(nz) + 0.5) / nz

    X = x[None, None, :]
    Y = y[None, :, None]
    Z = z[:, None, None]

    dx = X - bh_pos[0]
    dy = Y - bh_pos[1]

    vx = -omega_code * dy
    vy = omega_code * dx
    vz = np.zeros_like(vx)

    ge = energy_code_from_temp(temp_k)
    te = ge + 0.5 * (vx * vx + vy * vy + vz * vz)

    g["Density"][...] = rho_code
    g["GasEnergy"][...] = ge
    g["TotalEnergy"][...] = te
    g["x-velocity"][...] = vx
    g["y-velocity"][...] = vy
    g["z-velocity"][...] = vz
    g["Metal_Density"][...] = rho_code * metallicity_frac


def compute_expected_vrot_swirl(rho_code, bh_pos, bh_vel, omega_code):
    nx = ny = nz = 50
    dx_cell = 1.0 / nx
    i0 = int(math.floor((bh_pos[0] - 0.0) / dx_cell))
    j0 = int(math.floor((bh_pos[1] - 0.0) / dx_cell))
    k0 = int(math.floor((bh_pos[2] - 0.0) / dx_cell))
    i0 = max(0, min(nx - 1, i0))
    j0 = max(0, min(ny - 1, j0))
    k0 = max(0, min(nz - 1, k0))

    rcell = max(0, int(math.ceil(KERNEL_RADIUS_CODE / dx_cell)))
    r2 = KERNEL_RADIUS_CODE * KERNEL_RADIUS_CODE
    cell_volume = dx_cell ** 3

    lx = ly = lz = 0.0
    mass = 0.0
    for k in range(max(0, k0 - rcell), min(nz - 1, k0 + rcell) + 1):
        dz = (k - k0) * dx_cell
        zc = (k + 0.5) * dx_cell
        for j in range(max(0, j0 - rcell), min(ny - 1, j0 + rcell) + 1):
            dy = (j - j0) * dx_cell
            yc = (j + 0.5) * dx_cell
            for i in range(max(0, i0 - rcell), min(nx - 1, i0 + rcell) + 1):
                dx = (i - i0) * dx_cell
                if dx * dx + dy * dy + dz * dz > r2:
                    continue
                xc = (i + 0.5) * dx_cell
                mcell = rho_code * cell_volume
                if mcell <= 0:
                    continue

                vx = -omega_code * (yc - bh_pos[1])
                vy = omega_code * (xc - bh_pos[0])
                vz = 0.0
                dvx = vx - bh_vel[0]
                dvy = vy - bh_vel[1]
                dvz = vz - bh_vel[2]
                rx = xc - bh_pos[0]
                ry = yc - bh_pos[1]
                rz = zc - bh_pos[2]

                lx += mcell * (ry * dvz - rz * dvy)
                ly += mcell * (rz * dvx - rx * dvz)
                lz += mcell * (rx * dvy - ry * dvx)
                mass += mcell

    if mass <= 0:
        return 0.0
    lspec = math.sqrt((lx / mass) ** 2 + (ly / mass) ** 2 + (lz / mass) ** 2)
    return (lspec / KERNEL_RADIUS_CODE) * VEL_UNITS


def run_enzo(case_dir: Path, np_ranks: int, log_name: str):
    log_path = case_dir / log_name
    cmd = (
        f"cd {case_dir} && "
        f"{MODULE_INIT} && "
        f"export LD_LIBRARY_PATH={LD_PATH}:$LD_LIBRARY_PATH && "
        f"mpirun -np {np_ranks} {ENZO_EXE} -r DD0001/data0001 > {log_name} 2>&1"
    )
    proc = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
    return proc.returncode, log_path, proc.stdout + proc.stderr


def run_diag_case(name, configure_case_fn, check_fn, radiative=False, np_ranks=1):
    case_dir, h5path = prepare_case(name, radiative_cooling=radiative)
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        configure_case_fn(g)

    rc, log_path, shell_out = run_enzo(case_dir, np_ranks=np_ranks, log_name="run.log")
    if rc != 0:
        details = f"run failed rc={rc}; shell_out={shell_out[-2000:]}"
        return TestResult(name, False, details, {"rc": rc})

    try:
        bh, lines = get_last_bhaccr(log_path, bh_id=1)
    except Exception as e:
        return TestResult(name, False, f"parse failure: {e}", {})

    passed, details, metrics = check_fn(bh, lines, case_dir)
    return TestResult(name, passed, details, metrics)


# -----------------------------------------------------------------------------
# Test definitions
# -----------------------------------------------------------------------------

def make_uniform_case(rho, T, gas_vel=(0, 0, 0), bh_vel=(0, 0, 0), bh_pos=(0.15, 0.5, 0.5), bh_mass=None, Z=1e-6):
    def _cfg(g):
        set_particle_state(g, bh_pos=bh_pos, bh_vel=bh_vel, bh_mass_code=bh_mass)
        write_fields_uniform(g, rho, T, vel=gas_vel, metallicity_frac=Z)
    return _cfg


def test1_check(bh, lines, case_dir):
    ok = (
        bh["f_cold"] > 0.999 and bh["n_hot_cells"] == 0 and bh["n_cold_cells"] > 0 and
        almost_equal(bh["rho_cold_avg"], 1e-3 * DENSITY_UNITS, rel=2e-3)
    )
    return ok, f"f_cold={bh['f_cold']:.6g} rho_cold_avg={bh['rho_cold_avg']:.6e}", {
        "f_cold": bh["f_cold"], "rho_cold_avg": bh["rho_cold_avg"]
    }


def test2_config(g):
    # Build a deterministic 50/50 mixed kernel by painting half of the
    # spherical kernel cells hot in k-j-i order and the rest cold.
    bh_pos = (0.15, 0.5, 0.5)
    set_particle_state(g, bh_pos=bh_pos, bh_vel=(0, 0, 0))
    write_fields_uniform(g, rho_code=1e-3, temp_k=1e4, vel=(0, 0, 0), metallicity_frac=1e-6)

    nz, ny, nx = g["Density"].shape
    dx = 1.0 / nx
    i0 = int(math.floor(bh_pos[0] / dx))
    j0 = int(math.floor(bh_pos[1] / dx))
    k0 = int(math.floor(bh_pos[2] / dx))
    i0 = max(0, min(nx - 1, i0))
    j0 = max(0, min(ny - 1, j0))
    k0 = max(0, min(nz - 1, k0))

    rcell = max(0, int(math.ceil(KERNEL_RADIUS_CODE / dx)))
    r2 = KERNEL_RADIUS_CODE * KERNEL_RADIUS_CODE
    kernel_cells = []
    for k in range(max(0, k0 - rcell), min(nz - 1, k0 + rcell) + 1):
        dz = (k - k0) * dx
        for j in range(max(0, j0 - rcell), min(ny - 1, j0 + rcell) + 1):
            dy = (j - j0) * dx
            for i in range(max(0, i0 - rcell), min(nx - 1, i0 + rcell) + 1):
                ddx = (i - i0) * dx
                if ddx * ddx + dy * dy + dz * dz > r2:
                    continue
                kernel_cells.append((k, j, i))

    nhot = len(kernel_cells) // 2
    ge_hot = energy_code_from_temp(1e6)
    for n, (k, j, i) in enumerate(kernel_cells):
        if n < nhot:
            g["GasEnergy"][k, j, i] = ge_hot
            g["TotalEnergy"][k, j, i] = ge_hot


def test2_check(bh, lines, case_dir):
    ok = (
        0.40 < bh["f_hot"] < 0.60 and
        0.40 < bh["f_cold"] < 0.60 and
        bh["n_hot_cells"] > 0 and
        bh["n_cold_cells"] > 0 and
        bh["T_hot_avg"] > 5.0e5 and
        bh["T_cold_avg"] < 5.0e4
    )
    return ok, f"f_hot={bh['f_hot']:.4f} f_cold={bh['f_cold']:.4f} T_hot={bh['T_hot_avg']:.3e} T_cold={bh['T_cold_avg']:.3e}", {
        "f_hot": bh["f_hot"], "f_cold": bh["f_cold"], "T_hot_avg": bh["T_hot_avg"], "T_cold_avg": bh["T_cold_avg"]
    }


def test3_config(g):
    bh_pos = (0.15, 0.5, 0.5)
    set_particle_state(g, bh_pos=bh_pos, bh_vel=(0, 0, 0))
    write_fields_swirl_cold(g, rho_code=1e-3, temp_k=1e4, bh_pos=bh_pos, omega_code=8.0)


def test3_check(bh, lines, case_dir):
    expected = compute_expected_vrot_swirl(1e-3, bh_pos=(0.15, 0.5, 0.5), bh_vel=(0, 0, 0), omega_code=8.0)
    got = bh["V_rot_cold"]
    ok = almost_equal(got, expected, rel=8e-3)
    return ok, f"V_rot_cold={got:.6e} expected={expected:.6e}", {"V_rot_cold": got, "expected": expected}


def test4_check(bh, lines, case_dir):
    ok = (abs(bh["V_rot_cold"]) < 1e-8 and almost_equal(bh["f_AM"], 1.0, rel=1e-6))
    return ok, f"V_rot_cold={bh['V_rot_cold']:.3e} f_AM={bh['f_AM']:.6f}", {"V_rot_cold": bh["V_rot_cold"], "f_AM": bh["f_AM"]}


def test5_check(bh, lines, case_dir):
    ok = (bh["f_hot"] > 0.999 and bh["f_cold"] < 1e-9 and bh["Mdot_cold_raw"] == 0.0 and bh["n_fallback_cells"] == 0)
    return ok, f"f_hot={bh['f_hot']:.6f} n_fallback={bh['n_fallback_cells']}", {"f_hot": bh["f_hot"], "n_fallback_cells": bh["n_fallback_cells"]}


def test6_check(bh, lines, case_dir):
    ok = (bh["f_cold"] > 0.999 and bh["f_hot"] < 1e-9 and bh["Mdot_hot_raw"] == 0.0 and bh["n_fallback_cells"] == 0)
    return ok, f"f_cold={bh['f_cold']:.6f} n_fallback={bh['n_fallback_cells']}", {"f_cold": bh["f_cold"], "n_fallback_cells": bh["n_fallback_cells"]}


def test7_config(g):
    set_particle_state(g, bh_pos=(0.15, 0.5, 0.5), bh_vel=(0, 0, 0))
    write_fields_mixed_temperature(g, rho_code=1e-3, t_cold=1e4, t_hot=1e6, split_x=0.5)


def test7_check(bh, lines, case_dir):
    nk = int(bh["n_hot_cells"] + bh["n_cold_cells"])
    ok = (bh["n_fallback_cells"] == nk and bh["n_fallback_cells"] > 0)
    return ok, f"n_fallback={bh['n_fallback_cells']} n_kernel={nk}", {"n_fallback_cells": bh["n_fallback_cells"], "n_kernel_cells": nk}


def test8_config(g):
    set_particle_state(g, bh_pos=(0.15, 0.5, 0.5), bh_vel=(0, 0, 0))
    write_fields_mixed_temperature(g, rho_code=1e-3, t_cold=4.99e5, t_hot=5.01e5, split_x=0.5)


def test8_check_np(log1, log4):
    bh1, _ = get_last_bhaccr(log1, bh_id=1)
    bh4, _ = get_last_bhaccr(log4, bh_id=1)
    same_counts = (bh1["n_hot_cells"] == bh4["n_hot_cells"] and bh1["n_cold_cells"] == bh4["n_cold_cells"])
    same_fracs = almost_equal(bh1["f_hot"], bh4["f_hot"], rel=0.0, abs_tol=0.0)
    return same_counts and same_fracs, bh1, bh4


def test9_check(bh, lines, case_dir):
    m_bh = bh["bh_mass"] * MSUN
    rho = bh["rho_hot_avg"]
    cs = bh["cs_hot_avg"]
    vrel = bh["v_rel_hot"]
    denom = (cs * cs + vrel * vrel) ** 1.5
    mdot_bondi = 0.0 if denom <= 0 else 4 * math.pi * (G_CGS ** 2) * (m_bh ** 2) * rho / denom
    n_h = rho * XH / MH
    alpha = min(10.0, max(1.0, n_h / 0.1))
    mdot_expected = mdot_bondi * alpha * YR_S / MSUN
    got = bh["Mdot_hot_raw"]
    ok = almost_equal(got, mdot_expected, rel=1.5e-3)
    return ok, f"Mdot_hot_raw={got:.6e} expected={mdot_expected:.6e}", {"got": got, "expected": mdot_expected, "alpha": alpha}


def test10_check_pair(bh0, bhv):
    # expected multiplicative correction between vrel=0 and vrel>0
    cs = bh0["cs_hot_avg"]
    vrel = bhv["v_rel_hot"]
    expected_ratio = (cs * cs / (cs * cs + vrel * vrel)) ** 1.5
    got_ratio = bhv["Mdot_hot_raw"] / bh0["Mdot_hot_raw"] if bh0["Mdot_hot_raw"] > 0 else 0
    ok = almost_equal(got_ratio, expected_ratio, rel=2e-3)
    return ok, expected_ratio, got_ratio


def test11_check(bh, lines, case_dir):
    rho = bh["rho_hot_avg"]
    n_h = rho * XH / MH
    alpha = min(10.0, max(1.0, n_h / 0.1))
    ok = (bh["alpha_boost"] > 1.0 and bh["alpha_boost"] <= 10.0 and almost_equal(bh["alpha_boost"], alpha, rel=2e-3))
    return ok, f"alpha_boost={bh['alpha_boost']:.6g} expected={alpha:.6g}", {"alpha_boost": bh["alpha_boost"], "expected_alpha": alpha}


def test12_check(bh, lines, case_dir):
    cs = bh["cs_cold_avg"]
    vrot = bh["V_rot_cold"]
    expected = 1.0 if vrot <= 0.0 else min(1.0, (cs / vrot) ** 3 / BHACCR_CVISC)
    ok = (
        bh["Mdot_cold_raw"] > 0.0 and
        vrot > 0.0 and
        almost_equal(bh["f_AM"], expected, rel=2e-3)
    )
    return ok, f"f_AM={bh['f_AM']:.6e} expected={expected:.6e} V_rot_cold={vrot:.6e} cs_cold={cs:.6e}", {
        "f_AM": bh["f_AM"],
        "expected_f_AM": expected,
        "V_rot_cold": vrot,
        "cs_cold_avg": cs,
    }


def test13_check(bh, lines, case_dir):
    ok = (bh["Mdot_total_raw"] < bh["Mdot_Edd"] and bh["cap_active"] == 0 and almost_equal(bh["Mdot_actual"], bh["Mdot_total_raw"], rel=1e-8))
    return ok, f"Mdot_total_raw={bh['Mdot_total_raw']:.3e} Mdot_Edd={bh['Mdot_Edd']:.3e} cap={bh['cap_active']}", {
        "Mdot_total_raw": bh["Mdot_total_raw"], "Mdot_Edd": bh["Mdot_Edd"], "cap_active": bh["cap_active"]
    }


def test14_check(bh, lines, case_dir):
    ok = (bh["Mdot_total_raw"] > bh["Mdot_Edd"] and bh["cap_active"] == 1 and almost_equal(bh["Mdot_hot_actual"] + bh["Mdot_cold_actual"], bh["Mdot_Edd"], rel=1e-8))
    return ok, f"Mdot_total_raw={bh['Mdot_total_raw']:.3e} Mdot_Edd={bh['Mdot_Edd']:.3e} cap={bh['cap_active']}", {
        "Mdot_total_raw": bh["Mdot_total_raw"], "Mdot_Edd": bh["Mdot_Edd"], "cap_active": bh["cap_active"]
    }


def run_test0_baseline_reference():
    # Test 0 already executed on baseline worktree; capture artifact status.
    p = ROOT / "run/BHSeed/TS3_wrap/phaseA_test0_assert.log"
    ok = p.exists() and "TS3_WRAP: PASS" in p.read_text()
    details = "phaseA_test0_assert.log contains TS3_WRAP: PASS" if ok else "missing or failed baseline assert log"
    return TestResult("test0_baseline_ts3_wrap", ok, details, {"log": str(p)})


def run_test15_determinism():
    case_dir, h5path = prepare_case("t15_mpi_determinism", radiative_cooling=False)
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        set_particle_state(g, bh_pos=(0.15, 0.5, 0.5), bh_vel=(0.0, 0.0, 0.0))
        write_fields_uniform(g, rho_code=1e-3, temp_k=1e4, vel=(0, 0, 0), metallicity_frac=1e-6)

    rc1, log1, _ = run_enzo(case_dir, 1, "run_np1.log")
    rc4, log4, _ = run_enzo(case_dir, 4, "run_np4.log")
    if rc1 != 0 or rc4 != 0:
        return TestResult("test15_mpi_determinism", False, f"run rc np1={rc1} np4={rc4}", {"rc1": rc1, "rc4": rc4})

    bh1, _ = get_last_bhaccr(log1, bh_id=1)
    bh4, _ = get_last_bhaccr(log4, bh_id=1)

    keys = [
        "f_hot", "f_cold", "n_hot_cells", "n_cold_cells", "n_fallback_cells",
        "rho_hot_avg", "rho_cold_avg", "T_hot_avg", "T_cold_avg",
        "cs_hot_avg", "V_rot_cold", "v_rel_hot", "v_rel_cold",
        "Mdot_hot_raw", "Mdot_cold_raw", "Mdot_total_raw", "Mdot_Edd",
        "f_Edd", "alpha_boost", "f_AM", "Mdot_actual", "Mdot_hot_actual",
        "Mdot_cold_actual", "cap_active"
    ]
    mismatches = []
    for k in keys:
        a = bh1[k]
        b = bh4[k]
        if isinstance(a, int) and isinstance(b, int):
            if a != b:
                mismatches.append((k, a, b))
        else:
            if a != b:
                mismatches.append((k, a, b))

    ok = (len(mismatches) == 0)
    details = "all physics fields bit-identical (excluding wall_ms)" if ok else f"mismatches: {mismatches[:4]}"
    return TestResult("test15_mpi_determinism", ok, details, {"mismatch_count": len(mismatches)})


def run_test16_ts3_regression():
    # Fresh TS3_wrap re-run with current binary, then assert script.
    case_dir = OUT_ROOT / "t16_ts3_regression"
    ensure_clean_dir(case_dir)
    for fname in ["bhseed_ts3wrap.enzo", "make_ts3wrap_ic.py", "assert_ts3wrap.sh"]:
        shutil.copy2(TS3_DIR / fname, case_dir / fname)

    cmd = (
        f"cd {case_dir} && "
        f"python3 make_ts3wrap_ic.py > phaseA_t16_ic.log 2>&1 && "
        f"{MODULE_INIT} && "
        f"export LD_LIBRARY_PATH={LD_PATH}:$LD_LIBRARY_PATH && "
        f"mpirun -np 1 {ENZO_EXE} bhseed_ts3wrap.enzo > ts3wrap_np1.log 2>&1 && "
        f"mpirun -np 4 {ENZO_EXE} bhseed_ts3wrap.enzo > ts3wrap_np4.log 2>&1 && "
        f"bash assert_ts3wrap.sh > assert.log 2>&1"
    )
    rc = subprocess.run(["bash", "-lc", cmd]).returncode
    assert_log = (case_dir / "assert.log").read_text() if (case_dir / "assert.log").exists() else ""
    ok = (rc == 0 and "TS3_WRAP: PASS" in assert_log)
    details = "TS3_wrap assertions pass (8/8)" if ok else f"TS3 regression failed rc={rc}"
    return TestResult("test16_ts3_regression", ok, details, {"rc": rc})


def run_test17_multi_bh():
    # Multi-BH diagnostic run from restart: convert an existing non-BH particle
    # to MBH so we deterministically measure two BHs in one pass.
    case_dir, h5path = prepare_case("t17_multi_bh", radiative_cooling=False)
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        ptype = g["particle_type"][:]
        pmass = g["particle_mass"][:]
        pidx = g["particle_index"][:]
        px = g["particle_position_x"][:]
        py = g["particle_position_y"][:]
        pz = g["particle_position_z"][:]
        pvx = g["particle_velocity_x"][:]
        pvy = g["particle_velocity_y"][:]
        pvz = g["particle_velocity_z"][:]

        bh = np.where(ptype == 8)[0]
        if bh.size == 0 or ptype.size < 2:
            return TestResult("test17_multi_bh", False,
                              "restart does not contain enough particles to build 2-BH probe",
                              {"n_particles": int(ptype.size), "n_bh": int(bh.size)})

        i0 = int(bh[0])
        i1 = 0 if i0 != 0 else 1

        ptype[i0] = 8
        ptype[i1] = 8
        pidx[i0] = 1
        pidx[i1] = 2
        pmass[i1] = pmass[i0]

        px[i0], py[i0], pz[i0] = 0.15, 0.50, 0.50
        px[i1], py[i1], pz[i1] = 0.75, 0.50, 0.50
        pvx[i0] = pvy[i0] = pvz[i0] = 0.0
        pvx[i1] = pvy[i1] = pvz[i1] = 0.0

        g["particle_type"][...] = ptype
        g["particle_mass"][...] = pmass
        g["particle_index"][...] = pidx
        g["particle_position_x"][...] = px
        g["particle_position_y"][...] = py
        g["particle_position_z"][...] = pz
        g["particle_velocity_x"][...] = pvx
        g["particle_velocity_y"][...] = pvy
        g["particle_velocity_z"][...] = pvz

        write_fields_uniform(g, rho_code=1e-3, temp_k=1e4, vel=(0, 0, 0), metallicity_frac=1e-6)

    rc, log_path, _ = run_enzo(case_dir, np_ranks=1, log_name="run.log")
    if rc != 0:
        return TestResult("test17_multi_bh", False, f"run failed rc={rc}", {"rc": rc})

    lines = []
    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "[BHACCR]" in line:
                lines.append(parse_bhaccr_line(line))

    bh_ids = sorted({int(d.get("bh_id", -1)) for d in lines if int(d.get("bh_id", -1)) > 0})
    ok = len(bh_ids) >= 2
    return TestResult("test17_multi_bh", ok,
                      f"diagnostic BH ids: {bh_ids}",
                      {"bh_ids": bh_ids, "line_count": len(lines)})


def run_test18_perf(results):
    vals = []
    for r in results:
        v = r.metrics.get("accretion_diag_wall_ms")
        if isinstance(v, (int, float)):
            vals.append(float(v))
    if not vals:
        return TestResult("test18_performance", False, "no wall-time metrics collected", {})
    mean_ms = float(np.mean(vals))
    med_ms = float(np.median(vals))
    return TestResult("test18_performance", True, f"mean={mean_ms:.4f} ms median={med_ms:.4f} ms per BH diagnostic", {
        "mean_ms": mean_ms,
        "median_ms": med_ms,
        "samples": len(vals),
    })


def main():
    ensure_clean_dir(OUT_ROOT)
    t0 = time.time()

    results = []

    # Test 0 prerequisite baseline evidence
    results.append(run_test0_baseline_reference())

    # 1 Uniform gas property measurement (fallback cold)
    results.append(run_diag_case(
        "t1_uniform_measure",
        make_uniform_case(rho=1e-3, T=1e4),
        test1_check,
        radiative=False,
    ))

    # 2 Mixed hot/cold by temperature split (fallback)
    results.append(run_diag_case(
        "t2_mixed_temperature",
        test2_config,
        test2_check,
        radiative=False,
    ))

    # 3 Angular momentum computation
    results.append(run_diag_case(
        "t3_angular_momentum",
        test3_config,
        test3_check,
        radiative=False,
    ))

    # 4 Zero angular momentum
    results.append(run_diag_case(
        "t4_zero_angular_momentum",
        make_uniform_case(rho=1e-3, T=1e4, gas_vel=(0, 0, 0), bh_vel=(0, 0, 0)),
        test4_check,
        radiative=False,
    ))

    # 5 All-hot kernel via t_cool/t_dyn (radiative cooling on)
    results.append(run_diag_case(
        "t5_all_hot",
        make_uniform_case(rho=1e-3, T=1e7),
        test5_check,
        radiative=True,
    ))

    # 6 All-cold kernel via t_cool/t_dyn (radiative cooling on)
    results.append(run_diag_case(
        "t6_all_cold",
        make_uniform_case(rho=1.0, T=1e5),
        test6_check,
        radiative=True,
    ))

    # 7 Fallback participation and count
    results.append(run_diag_case(
        "t7_fallback",
        test7_config,
        test7_check,
        radiative=False,
    ))

    # 8 Near-threshold stability (temperature floor boundary) np1 vs np4
    case_dir, h5path = prepare_case("t8_threshold_stability", radiative_cooling=False)
    with h5py.File(h5path, "r+") as f:
        g = f["Grid00000001"]
        test8_config(g)
    rc1, log1, _ = run_enzo(case_dir, 1, "run_np1.log")
    rc4, log4, _ = run_enzo(case_dir, 4, "run_np4.log")
    if rc1 == 0 and rc4 == 0:
        ok, bh1, bh4 = test8_check_np(log1, log4)
        details = f"np1 hot/cold={bh1['n_hot_cells']}/{bh1['n_cold_cells']} np4 hot/cold={bh4['n_hot_cells']}/{bh4['n_cold_cells']}"
        results.append(TestResult("t8_threshold_stability", ok, details, {
            "np1_hot": bh1["n_hot_cells"],
            "np1_cold": bh1["n_cold_cells"],
            "np4_hot": bh4["n_hot_cells"],
            "np4_cold": bh4["n_cold_cells"],
        }))
    else:
        results.append(TestResult("t8_threshold_stability", False, f"run rc np1={rc1} np4={rc4}", {"rc1": rc1, "rc4": rc4}))

    # 9 Bondi analytic (hot, v_rel=0)
    r9 = run_diag_case(
        "t9_bondi_analytic",
        make_uniform_case(rho=1e-3, T=1e7, gas_vel=(0, 0, 0), bh_vel=(0, 0, 0)),
        test9_check,
        radiative=False,
    )
    results.append(r9)

    # 10 Bondi-Hoyle correction (compare to v_rel=0 case)
    r10_capture = run_diag_case(
        "t10_bondi_hoyle",
        make_uniform_case(rho=1e-3, T=1e7, gas_vel=(0.02, 0, 0), bh_vel=(0, 0, 0)),
        lambda bh, lines, case: (True, "captured", {"bh": bh}),
        radiative=False,
    )
    if not r10_capture.passed:
        results.append(TestResult("t10_bondi_hoyle", False, r10_capture.details, r10_capture.metrics))
    else:
        try:
            bh9, _ = get_last_bhaccr(OUT_ROOT / "t9_bondi_analytic" / "run.log", bh_id=1)
            bh10, _ = get_last_bhaccr(OUT_ROOT / "t10_bondi_hoyle" / "run.log", bh_id=1)
            ok10, exp_ratio, got_ratio = test10_check_pair(bh9, bh10)
            results.append(TestResult("t10_bondi_hoyle", ok10,
                                      f"ratio={got_ratio:.6e} expected={exp_ratio:.6e}",
                                      {"got_ratio": got_ratio, "expected_ratio": exp_ratio}))
        except Exception as e:
            results.append(TestResult("t10_bondi_hoyle", False, f"parse failure: {e}", {}))

    # 11 Boost factor
    results.append(run_diag_case(
        "t11_boost_factor",
        make_uniform_case(rho=1.0, T=1e7, gas_vel=(0, 0, 0), bh_vel=(0, 0, 0)),
        test11_check,
        radiative=False,
    ))

    # 12 AM suppression
    def t12_cfg(g):
        bh_pos = (0.15, 0.5, 0.5)
        set_particle_state(g, bh_pos=bh_pos, bh_vel=(0, 0, 0))
        write_fields_swirl_cold(g, rho_code=1e-2, temp_k=1e4, bh_pos=bh_pos, omega_code=20.0)
    results.append(run_diag_case(
        "t12_am_suppression",
        t12_cfg,
        test12_check,
        radiative=False,
    ))

    # 13 Sub-Eddington
    results.append(run_diag_case(
        "t13_sub_eddington",
        make_uniform_case(rho=1e-5, T=1e7),
        test13_check,
        radiative=False,
    ))

    # 14 Super-Eddington diagnostic
    results.append(run_diag_case(
        "t14_super_eddington",
        make_uniform_case(rho=20.0, T=1e6, bh_mass=5.06744516e-6),
        test14_check,
        radiative=False,
    ))

    # 15 MPI determinism
    results.append(run_test15_determinism())

    # 16 TS3_wrap regression
    results.append(run_test16_ts3_regression())

    # 17 multi-BH diagnostic
    results.append(run_test17_multi_bh())

    # Collect wall times from all single-case logs for Test 18
    for r in results:
        case_log = OUT_ROOT / r.name / "run.log"
        if case_log.exists():
            try:
                bh, _ = get_last_bhaccr(case_log, bh_id=1)
                if "accretion_diag_wall_ms" in bh:
                    r.metrics["accretion_diag_wall_ms"] = float(bh["accretion_diag_wall_ms"])
            except Exception:
                pass
    results.append(run_test18_perf(results))

    elapsed = time.time() - t0

    out = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_s": elapsed,
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

    with open(OUT_ROOT / "phaseA_matrix_results.json", "w") as f:
        json.dump(out, f, indent=2)

    # Human-readable summary
    lines = []
    lines.append("=== Phase A Matrix Summary ===")
    lines.append(f"Output dir: {OUT_ROOT}")
    pass_count = sum(1 for r in results if r.passed)
    lines.append(f"Passed: {pass_count}/{len(results)}")
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        lines.append(f"{r.name}: {status} :: {r.details}")

    (OUT_ROOT / "phaseA_matrix_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
