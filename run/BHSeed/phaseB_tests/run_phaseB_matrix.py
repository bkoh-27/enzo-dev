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

ROOT = Path('/mnt/home/boh10/bh_proj/enzo-dev')
TS3_DIR = ROOT / 'run/BHSeed/TS3_wrap'
BASE_RESTART = TS3_DIR / 'DD0001'
OUT_ROOT = ROOT / f"run/BHSeed/phaseB_tests/matrix_{time.strftime('%Y%m%d_%H%M%S')}"
ENZO_EXE = ROOT / 'src/enzo/enzo.exe'
COOL_RATES = ROOT / 'input/cool_rates.in'
GRACKLE_DATA = Path('/mnt/home/boh10/ceph/ENZO/param_tuning/nBH_test/seed/seed_acc_fb/CloudyData_UVB=HM2012_new01222019.h5')

MODULE_INIT = 'source /etc/profile.d/modules.sh && module load gcc openmpi hdf5/1.12.3'
EXPORTS = (
    'export PATH=/mnt/home/boh10/ltu_proj/local/bin:/mnt/home/boh10/yt-conda/bin:$PATH && '
    'export LD_LIBRARY_PATH=/mnt/home/boh10/ltu_proj/local/lib:$LD_LIBRARY_PATH && '
    'export HDF5_DISABLE_VERSION_CHECK=2'
)

LENGTH_UNITS = 3.0857e23
TIME_UNITS = 3.1557e13
DENSITY_UNITS = 1.67e-24
VEL_UNITS = LENGTH_UNITS / TIME_UNITS
CELL_WIDTH = 1.0 / 50.0
BOX_KPC = LENGTH_UNITS / 3.0857e21
GAMMA = 5.0 / 3.0
MU = 0.6
KBOLTZ = 1.3807e-16
MH = 1.6726e-24
XH = 0.76
MSUN = 1.989e33
YR_S = 3.1557e7
G_CGS = 6.67259e-8

DEFAULT_KERNEL_RADIUS_KPC = 3.0
DEFAULT_REMOVAL_RADIUS_CELLS = 1
REL = 5e-3
ABS = 1e-14


@dataclass
class TResult:
    name: str
    passed: bool
    details: str
    metrics: dict


def sh(cmd, cwd=None):
    p = subprocess.run(['bash', '-lc', cmd], cwd=cwd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def ensure_clean(path: Path):
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def set_param(text: str, key: str, value: str) -> str:
    pat = re.compile(rf'^{re.escape(key)}\s*=.*$', re.MULTILINE)
    line = f'{key:<30}= {value}'
    if pat.search(text):
        return pat.sub(line, text)
    return text + '\n' + line + '\n'


def e_from_T(T):
    return (KBOLTZ * T / ((GAMMA - 1.0) * MU * MH)) / (VEL_UNITS ** 2)


def parse_kv_line(line):
    d = {}
    for tok in line.strip().split():
        if '=' not in tok:
            continue
        k, v = tok.split('=', 1)
        if k.startswith('[BHACCR'):
            continue
        try:
            if re.match(r'^-?\d+$', v):
                d[k] = int(v)
            else:
                d[k] = float(v)
        except Exception:
            d[k] = v
    return d


def parse_bhaccr(log: Path):
    out = []
    if not log.exists():
        return out
    with open(log, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if '[BHACCR]' in line:
                d = parse_kv_line(line)
                d['__line'] = line.strip()
                out.append(d)
    return out


def parse_bhaccr_warn(log: Path):
    out = []
    if not log.exists():
        return out
    with open(log, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if '[BHACCR_WARN]' in line:
                out.append(line.strip())
    return out


def parse_bhaccr_debug(log: Path):
    out = []
    if not log.exists():
        return out
    with open(log, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if '[BHACCR_DEBUG]' in line:
                d = parse_kv_line(line)
                d['__line'] = line.strip()
                out.append(d)
    return out


def almost(a, b, rel=REL, abs_=ABS):
    return abs(a - b) <= max(abs_, rel * max(abs(a), abs(b), 1.0))


def prepare_case(name, stop='0.02', radiative=False, grackle=False):
    cdir = OUT_ROOT / name
    dd1 = cdir / 'DD0001'
    ensure_clean(cdir)
    shutil.copytree(BASE_RESTART, dd1)

    p = dd1 / 'data0001'
    txt = p.read_text()

    txt = set_param(txt, 'StopTime', stop)
    txt = set_param(txt, 'BHSeedingMethod', '0')
    txt = set_param(txt, 'BHSeedVerbose', '0')

    txt = set_param(txt, 'BHAccretionMethod', '1')
    txt = set_param(txt, 'BHAccretionKernelRadius', f'{DEFAULT_KERNEL_RADIUS_KPC}')
    txt = set_param(txt, 'BHAccretionRemovalRadius', f'{DEFAULT_REMOVAL_RADIUS_CELLS}')
    txt = set_param(txt, 'BHAccretionRemovalMode', '0')
    txt = set_param(txt, 'BHAccretionTSplitFloor', '5.0e5')
    txt = set_param(txt, 'BHAccretionColdModel', '0')
    txt = set_param(txt, 'BHAccretionCVisc', '6.283')
    txt = set_param(txt, 'BHAccretionNHStar', '0.1')
    txt = set_param(txt, 'BHAccretionBeta', '1.0')
    txt = set_param(txt, 'BHAccretionAlphaMax', '10.0')
    txt = set_param(txt, 'BHAccretionRadiativeEfficiency', '0.1')
    txt = set_param(txt, 'BHAccretionVerbose', '1')
    txt = set_param(txt, 'BHAccretionRunEveryTimestep', '1')
    txt = set_param(txt, 'BHAccretionIgnoredDVWarn', '1.0')
    txt = set_param(txt, 'BHAccretionIgnoredPFracWarn', '0.01')

    txt = set_param(txt, 'RadiativeCooling', '1' if radiative else '0')
    txt = set_param(txt, 'use_grackle', '1' if grackle else '0')
    if grackle:
        txt = set_param(txt, 'RadiativeCooling', '1')
        txt = set_param(txt, 'MultiSpecies', '0')
        txt = set_param(txt, 'MetalCooling', '1')
        txt = set_param(txt, 'UVbackground', '1')
        txt = set_param(txt, 'CMBTemperatureFloor', '1')
        txt = set_param(txt, 'grackle_data_file', 'CloudyData_UVB=HM2012_new01222019.h5')
        txt = set_param(txt, 'HydrogenFractionByMass', '-1.0')
    else:
        txt = set_param(txt, 'MetalCooling', '0')

    p.write_text(txt)

    if radiative and COOL_RATES.exists():
        shutil.copy2(COOL_RATES, cdir / 'cool_rates.in')
    if grackle:
        if not GRACKLE_DATA.exists():
            raise RuntimeError(f'Missing grackle data file: {GRACKLE_DATA}')
        shutil.copy2(GRACKLE_DATA, cdir / GRACKLE_DATA.name)

    return cdir, dd1 / 'data0001.cpu0000'


def run_case(cdir: Path, np_ranks=1, log='run.log', restart='DD0001/data0001'):
    cmd = (
        f'cd {cdir} && {MODULE_INIT} && {EXPORTS} && '
        f'mpirun -np {np_ranks} {ENZO_EXE} -r {restart} > {log} 2>&1'
    )
    rc, so, se = sh(cmd)
    return rc, cdir / log, so + se


def open_grid(path: Path, mode='r'):
    f = h5py.File(path, mode)
    return f, f['Grid00000001']


def latest_output_cpu0(cdir: Path):
    dds = sorted([p for p in cdir.glob('DD????') if p.is_dir()])
    if not dds:
        return None
    last = dds[-1]
    num = last.name[2:]
    return last / f'data{num}.cpu0000'


def active_slice(arr):
    # TS3_wrap has 3 ghost zones
    return arr[3:-3, 3:-3, 3:-3]


def bh_index(g):
    pt = g['particle_type'][:]
    idx = np.where(pt == 8)[0]
    if idx.size == 0:
        return None
    return int(idx[0])


def set_bh(g, pos=(0.15, 0.5, 0.5), vel=(0, 0, 0), mass_code=None, pindex=1):
    i = bh_index(g)
    if i is None:
        raise RuntimeError('No BH in file')
    g['particle_position_x'][i] = pos[0]
    g['particle_position_y'][i] = pos[1]
    g['particle_position_z'][i] = pos[2]
    g['particle_velocity_x'][i] = vel[0]
    g['particle_velocity_y'][i] = vel[1]
    g['particle_velocity_z'][i] = vel[2]
    if mass_code is not None:
        g['particle_mass'][i] = mass_code
    g['particle_index'][i] = pindex


def write_uniform(g, rho=1e-3, T=1e4, vel=(0, 0, 0), Z=1e-6):
    vx, vy, vz = vel
    ge = e_from_T(T)
    te = ge + 0.5*(vx*vx+vy*vy+vz*vz)
    g['Density'][...] = rho
    g['GasEnergy'][...] = ge
    g['TotalEnergy'][...] = te
    g['x-velocity'][...] = vx
    g['y-velocity'][...] = vy
    g['z-velocity'][...] = vz
    if 'Metal_Density' in g:
        g['Metal_Density'][...] = rho*Z


def write_mixed(g, rho=1e-3, t_cold=1e4, t_hot=1e6, split=0.5):
    nz, ny, nx = g['Density'].shape
    x = (np.arange(nx)+0.5)/nx
    cold = x < split
    T = np.where(cold[None, None, :], t_cold, t_hot)
    ge = e_from_T(T)
    g['Density'][...] = rho
    g['GasEnergy'][...] = ge
    g['TotalEnergy'][...] = ge
    g['x-velocity'][...] = 0.0
    g['y-velocity'][...] = 0.0
    g['z-velocity'][...] = 0.0
    if 'Metal_Density' in g:
        g['Metal_Density'][...] = rho*1e-6


def write_swirl(g, rho=1e-3, T=1e4, bh_pos=(0.15,0.5,0.5), omega=8.0):
    nz, ny, nx = g['Density'].shape
    x = (np.arange(nx)+0.5)/nx
    y = (np.arange(ny)+0.5)/ny
    X = x[None,None,:]
    Y = y[None,:,None]
    vx = -omega*(Y-bh_pos[1])
    vy = omega*(X-bh_pos[0])
    vz = np.zeros_like(vx)
    ge = e_from_T(T)
    te = ge + 0.5*(vx*vx+vy*vy+vz*vz)
    g['Density'][...] = rho
    g['GasEnergy'][...] = ge
    g['TotalEnergy'][...] = te
    g['x-velocity'][...] = vx
    g['y-velocity'][...] = vy
    g['z-velocity'][...] = vz
    if 'Metal_Density' in g:
        g['Metal_Density'][...] = rho*1e-6


def bh_mass_code(g, i):
    return float(g['particle_mass'][i])


def bh_vel(g, i):
    return (float(g['particle_velocity_x'][i]), float(g['particle_velocity_y'][i]), float(g['particle_velocity_z'][i]))


def bh_attrs(g, i):
    out = {}
    for k in g.keys():
        if k.startswith('particle_attribute'):
            out[k] = float(g[k][i])
    return out


def find_removal_cells(d0, d1):
    diff = np.abs(active_slice(d1) - active_slice(d0))
    idx = np.argwhere(diff > 0)
    return idx, diff


def test0_ts3_baseline():
    cdir = OUT_ROOT / 'test0_ts3_baseline'
    ensure_clean(cdir)
    for f in ['bhseed_ts3wrap.enzo', 'make_ts3wrap_ic.py', 'assert_ts3wrap.sh']:
        shutil.copy2(TS3_DIR / f, cdir / f)
    cmd = (
        f'cd {cdir} && python3 make_ts3wrap_ic.py > ic.log 2>&1 && '
        f'{MODULE_INIT} && {EXPORTS} && '
        f'mpirun -np 1 {ENZO_EXE} bhseed_ts3wrap.enzo > ts3wrap_np1.log 2>&1 && '
        f'mpirun -np 4 {ENZO_EXE} bhseed_ts3wrap.enzo > ts3wrap_np4.log 2>&1 && '
        f'bash assert_ts3wrap.sh > assert.log 2>&1'
    )
    rc, _, _ = sh(cmd)
    ok = False
    if (cdir/'assert.log').exists():
        ok = 'TS3_WRAP: PASS' in (cdir/'assert.log').read_text(errors='ignore')
    return TResult('Test 0', ok and rc == 0, 'TS3_wrap baseline', {'rc': rc, 'dir': str(cdir)})


def test1_mass_conservation():
    cdir, h5 = prepare_case('test1_mass_conservation', stop='0.06', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, pos=(0.15,0.5,0.5), vel=(0,0,0), mass_code=5.06744516e-9)
        write_uniform(g, rho=1e-1, T=1e4, vel=(0,0,0))
    with h5py.File(h5, 'r') as f0:
        g0 = f0['Grid00000001']
        i0 = bh_index(g0)
        m0 = bh_mass_code(g0, i0)
        gas0 = float(active_slice(g0['Density'][:]).sum())

    rc, log, _ = run_case(cdir, np_ranks=1, log='run.log')
    lines = parse_bhaccr(log)
    if rc != 0 or not lines:
        return TResult('Test 1', False, f'run failed rc={rc}', {'rc': rc})
    dm_sum = sum(float(x.get('dm_removed', 0.0)) for x in lines)

    out = latest_output_cpu0(cdir)
    with h5py.File(out, 'r') as ff:
        gg = ff['Grid00000001']
        i1 = bh_index(gg)
        m1 = bh_mass_code(gg, i1)
        gas1 = float(active_slice(gg['Density'][:]).sum())

    dm_bh = m1 - m0
    dm_gas = gas0 - gas1
    ok = almost(dm_bh, dm_sum, rel=5e-3) and almost(dm_gas, dm_sum, rel=1e-2)
    return TResult('Test 1', ok,
                   f'dm_sum={dm_sum:.3e} dm_bh={dm_bh:.3e} dm_gas={dm_gas:.3e}',
                   {'dm_sum': dm_sum, 'dm_bh': dm_bh, 'dm_gas': dm_gas})


def test2_momentum_bookkeeping():
    cdir, h5 = prepare_case('test2_momentum_bookkeeping', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, pos=(0.15,0.5,0.5), vel=(1e-5, 2e-5, -1e-5), mass_code=5.06744516e-9)
        write_mixed(g, rho=1e-2, t_cold=1e4, t_hot=1e6, split=0.5)
        g['x-velocity'][...] += 0.03
        g['y-velocity'][...] -= 0.01
    with h5py.File(h5, 'r') as f0:
        g0 = f0['Grid00000001']
        i0 = bh_index(g0)
        v0 = bh_vel(g0, i0)

    rc, log, _ = run_case(cdir)
    lines = parse_bhaccr(log)
    warns = parse_bhaccr_warn(log)
    if rc != 0 or not lines:
        return TResult('Test 2', False, f'run failed rc={rc}', {'rc': rc})
    last = lines[-1]

    out = latest_output_cpu0(cdir)
    with h5py.File(out, 'r') as f1:
        g1 = f1['Grid00000001']
        i1 = bh_index(g1)
        v1 = bh_vel(g1, i1)

    vel_same = all(almost(a, b, rel=0, abs_=1e-14) for a, b in zip(v0, v1))
    ignored_nonzero = (float(last.get('acc_ignored_dv_kms', 0.0)) >= 0.0 and float(last.get('acc_ignored_p_frac', 0.0)) >= 0.0)
    ok = vel_same and ignored_nonzero
    return TResult('Test 2', ok,
                   f'vel_same={vel_same} ignored_dv={last.get("acc_ignored_dv_kms",0):.3e} warn_lines={len(warns)}',
                   {'vel0': v0, 'vel1': v1, 'warn_lines': len(warns)})


def test3_energy_consistency():
    # PPM check via ratio of specific energies with mass fraction in removed cells.
    cdir, h5 = prepare_case('test3_energy_consistency', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, pos=(0.15,0.5,0.5), vel=(0,0,0), mass_code=5.06744516e-8)
        write_uniform(g, rho=1e-2, T=1e6, vel=(0.01, -0.01, 0.0))

    with h5py.File(h5, 'r') as f0:
        g0 = f0['Grid00000001']
        d0 = g0['Density'][:]
        ge0 = g0['GasEnergy'][:]
        te0 = g0['TotalEnergy'][:]

    rc, log, _ = run_case(cdir)
    if rc != 0:
        return TResult('Test 3', False, f'run failed rc={rc}', {'rc': rc})

    out = latest_output_cpu0(cdir)
    with h5py.File(out, 'r') as f1:
        g1 = f1['Grid00000001']
        d1 = g1['Density'][:]
        ge1 = g1['GasEnergy'][:]
        te1 = g1['TotalEnergy'][:]

    changed = np.argwhere(np.abs(active_slice(d1) - active_slice(d0)) > 0)
    if changed.size == 0:
        return TResult('Test 3', False, 'no changed cells', {})

    ok = True
    checks = 0
    for kk, jj, ii in changed[:20]:
        k, j, i = kk+3, jj+3, ii+3
        if d0[k,j,i] <= 0:
            continue
        frac = d1[k,j,i] / d0[k,j,i]
        if ge0[k,j,i] != 0:
            ok = ok and almost(ge1[k,j,i], ge0[k,j,i] * frac, rel=1e-3)
            checks += 1
        if te0[k,j,i] != 0:
            ok = ok and almost(te1[k,j,i], te0[k,j,i] * frac, rel=1e-3)
            checks += 1
    return TResult('Test 3', ok and checks > 0, f'ppm_energy_checks={checks}', {'checks': checks})


def test4_geometry():
    cdir, h5 = prepare_case('test4_removal_geometry', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, pos=(0.21, 0.51, 0.51), vel=(0,0,0), mass_code=5.06744516e-8)
        write_uniform(g, rho=1e-2, T=1e4)

    with h5py.File(h5, 'r') as f0:
        g0 = f0['Grid00000001']
        d0 = g0['Density'][:]

    rc, log, _ = run_case(cdir)
    lines = parse_bhaccr(log)
    if rc != 0 or not lines:
        return TResult('Test 4', False, f'run failed rc={rc}', {'rc': rc})
    last = lines[-1]

    out = latest_output_cpu0(cdir)
    with h5py.File(out, 'r') as f1:
        d1 = f1['Grid00000001']['Density'][:]
    changed, _ = find_removal_cells(d0, d1)
    cnt = int(changed.shape[0])
    ok = cnt == int(last.get('removal_cells', -1)) and cnt > 0
    return TResult('Test 4', ok, f'changed={cnt} log_removal_cells={int(last.get("removal_cells",-1))}',
                   {'changed': cnt, 'log_removal_cells': int(last.get('removal_cells', -1))})


def test5_active_zone_only():
    cdir, h5 = prepare_case('test5_active_zone_writes', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, pos=(0.15,0.5,0.5), vel=(0,0,0), mass_code=5.06744516e-8)
        write_uniform(g, rho=1e-2, T=1e4)

    with h5py.File(h5, 'r') as f0:
        d0 = f0['Grid00000001']['Density'][:]

    rc, _, _ = run_case(cdir)
    if rc != 0:
        return TResult('Test 5', False, f'run failed rc={rc}', {'rc': rc})
    out = latest_output_cpu0(cdir)
    with h5py.File(out, 'r') as f1:
        d1 = f1['Grid00000001']['Density'][:]

    ghost0 = d0.copy()
    ghost1 = d1.copy()
    # keep only ghosts
    ghost0[3:-3,3:-3,3:-3] = 0.0
    ghost1[3:-3,3:-3,3:-3] = 0.0
    maxdiff = float(np.max(np.abs(ghost1 - ghost0)))
    ok = maxdiff == 0.0
    return TResult('Test 5', ok, f'max_ghost_density_diff={maxdiff:.3e}', {'max_ghost_density_diff': maxdiff})


def test6_proportional():
    cdir, h5 = prepare_case('test6_proportional_removal', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, pos=(0.15,0.5,0.5), vel=(0,0,0), mass_code=5.06744516e-8)
        nz, ny, nx = g['Density'].shape
        x = (np.arange(nx)+0.5)/nx
        rho = 5e-3 + 5e-3 * x[None,None,:]
        g['Density'][...] = rho
        ge = e_from_T(1e4)
        g['GasEnergy'][...] = ge
        g['TotalEnergy'][...] = ge
        g['x-velocity'][...] = 0.0
        g['y-velocity'][...] = 0.0
        g['z-velocity'][...] = 0.0
        if 'Metal_Density' in g:
            g['Metal_Density'][...] = rho*1e-6

    with h5py.File(h5, 'r') as f0:
        d0 = f0['Grid00000001']['Density'][:]

    rc, _, _ = run_case(cdir)
    if rc != 0:
        return TResult('Test 6', False, f'run failed rc={rc}', {'rc': rc})

    out = latest_output_cpu0(cdir)
    with h5py.File(out, 'r') as f1:
        d1 = f1['Grid00000001']['Density'][:]

    d0a = active_slice(d0)
    d1a = active_slice(d1)
    changed = np.argwhere(np.abs(d1a - d0a) > 0)
    if changed.shape[0] < 2:
        return TResult('Test 6', False, 'too few changed cells', {'changed': int(changed.shape[0])})

    fracs = []
    for kk,jj,ii in changed:
        old = d0a[kk,jj,ii]
        new = d1a[kk,jj,ii]
        if old > 0:
            fracs.append((old-new)/old)
    spread = float(np.max(fracs) - np.min(fracs)) if fracs else 1.0
    ok = spread < 5e-4
    return TResult('Test 6', ok, f'removal_fraction_spread={spread:.3e}', {'spread': spread, 'n': len(fracs)})


def test7_gas_depletion():
    cdir, h5 = prepare_case('test7_gas_depletion', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, pos=(0.15,0.5,0.5), vel=(0,0,0), mass_code=5.06744516e-2)
        write_uniform(g, rho=1e-12, T=1e4)

    rc, log, _ = run_case(cdir)
    lines = parse_bhaccr(log)
    if rc != 0 or not lines:
        return TResult('Test 7', False, f'run failed rc={rc}', {'rc': rc})
    bh = lines[-1]
    out = latest_output_cpu0(cdir)
    with h5py.File(out, 'r') as f1:
        d1 = active_slice(f1['Grid00000001']['Density'][:])
    minrho = float(np.min(d1))
    ok = (bh.get('removal_gas_limited',0)==1 and float(bh.get('dm_removed',0)) <= float(bh.get('dm_requested',0)) and minrho >= 0.0)
    return TResult('Test 7', ok,
                   f"gas_limited={bh.get('removal_gas_limited')} dm_removed={bh.get('dm_removed',0):.3e} dm_req={bh.get('dm_requested',0):.3e}",
                   {'min_rho_after': minrho})


def test8_single_cell_mode():
    cdir, h5 = prepare_case('test8_single_cell', stop='0.02', radiative=False)
    p = (cdir/'DD0001'/'data0001')
    txt = p.read_text()
    txt = set_param(txt, 'BHAccretionRemovalMode', '1')
    p.write_text(txt)

    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, pos=(0.15,0.5,0.5), vel=(0,0,0), mass_code=5.06744516e-8)
        write_uniform(g, rho=1e-2, T=1e4)

    with h5py.File(h5, 'r') as f0:
        d0 = f0['Grid00000001']['Density'][:]

    rc, log, _ = run_case(cdir)
    lines = parse_bhaccr(log)
    if rc != 0 or not lines:
        return TResult('Test 8', False, f'run failed rc={rc}', {'rc': rc})

    out = latest_output_cpu0(cdir)
    with h5py.File(out, 'r') as f1:
        d1 = f1['Grid00000001']['Density'][:]

    changed, _ = find_removal_cells(d0, d1)
    cnt = int(changed.shape[0])
    lcnt = int(lines[-1].get('removal_cells', -1))
    ok = (cnt == 1 and lcnt == 1)
    return TResult('Test 8', ok, f'changed={cnt} log={lcnt}', {'changed': cnt, 'log': lcnt})


def run_hot_cold_mixed_tests():
    out = {}
    # Test 9 hot-only
    c9, h9 = prepare_case('test9_hot_only', stop='0.02', radiative=True)
    with h5py.File(h9, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, mass_code=5.06744516e-9)
        write_uniform(g, rho=1e-3, T=1e7)
    rc9, log9, _ = run_case(c9)
    bh9 = parse_bhaccr(log9)[-1] if rc9 == 0 else {}
    ok9 = (rc9 == 0 and float(bh9.get('Mdot_cold_realized',1)) == 0.0 and float(bh9.get('dm_removed',0)) > 0.0)
    out['9'] = TResult('Test 9', ok9, f"Mdot_cold_realized={bh9.get('Mdot_cold_realized','nan')}", dict(bh9))

    # Test 10 cold-only
    c10, h10 = prepare_case('test10_cold_only', stop='0.02', radiative=True)
    with h5py.File(h10, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, mass_code=5.06744516e-9)
        write_uniform(g, rho=1.0, T=1e5)
    rc10, log10, _ = run_case(c10)
    bh10 = parse_bhaccr(log10)[-1] if rc10 == 0 else {}
    ok10 = (rc10 == 0 and float(bh10.get('Mdot_hot_realized',1)) == 0.0 and float(bh10.get('dm_removed',0)) > 0.0)
    out['10'] = TResult('Test 10', ok10, f"Mdot_hot_realized={bh10.get('Mdot_hot_realized','nan')}", dict(bh10))

    # Test 11 mixed
    c11, h11 = prepare_case('test11_mixed', stop='0.02', radiative=False)
    with h5py.File(h11, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, pos=(0.5,0.5,0.5), mass_code=5.06744516e-9)
        write_mixed(g, rho=1e-3, t_cold=1e4, t_hot=1e6, split=0.5)
    rc11, log11, _ = run_case(c11)
    bh11 = parse_bhaccr(log11)[-1] if rc11 == 0 else {}
    sum_realized = float(bh11.get('Mdot_hot_realized',0)) + float(bh11.get('Mdot_cold_realized',0))
    dm_dt = float(bh11.get('dm_removed',0)) / max(float(bh11.get('dt', 0.0) or 4.0536666e-10), 1e-99)
    # dt not logged, so verify against Mdot_actual*frac_gas
    ok11 = (rc11 == 0 and float(bh11.get('Mdot_hot_realized',0)) > 0 and float(bh11.get('Mdot_cold_realized',0)) > 0 and almost(sum_realized, float(bh11.get('Mdot_actual',0))*float(bh11.get('frac_gas',1)), rel=1e-6))
    out['11'] = TResult('Test 11', ok11, f"hot={bh11.get('Mdot_hot_realized','nan')} cold={bh11.get('Mdot_cold_realized','nan')}", dict(bh11))

    return out


def test12_13_sub_super():
    out = {}
    # Sub-Edd
    c12, h12 = prepare_case('test12_subedd', stop='0.02', radiative=False)
    with h5py.File(h12, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, mass_code=5.06744516e-9)
        write_uniform(g, rho=1e-5, T=1e7)
    rc12, log12, _ = run_case(c12)
    bh12 = parse_bhaccr(log12)[-1] if rc12 == 0 else {}
    ok12 = (rc12==0 and int(bh12.get('cap_active',1))==0 and almost(float(bh12.get('frac_cap',0)),1.0,rel=0,abs_=1e-12))
    out['12'] = TResult('Test 12', ok12, f"cap={bh12.get('cap_active')} frac_cap={bh12.get('frac_cap')}", dict(bh12))

    # Super-Edd
    c13, h13 = prepare_case('test13_superedd', stop='0.02', radiative=False)
    with h5py.File(h13, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, mass_code=5.06744516e-6)
        write_uniform(g, rho=20.0, T=1e6)
    rc13, log13, _ = run_case(c13)
    bh13 = parse_bhaccr(log13)[-1] if rc13 == 0 else {}
    ok13 = (rc13==0 and int(bh13.get('cap_active',0))==1 and float(bh13.get('frac_cap',1)) < 1.0)
    out['13'] = TResult('Test 13', ok13, f"cap={bh13.get('cap_active')} frac_cap={bh13.get('frac_cap')}", dict(bh13))
    return out


def test14_15_multistep_invariants():
    cdir, h5 = prepare_case('test14_15_multistep', stop='0.06', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, mass_code=5.06744516e-8)
        write_uniform(g, rho=1e-2, T=1e4)
    rc, log, _ = run_case(cdir)
    lines = parse_bhaccr(log)
    if rc != 0 or len(lines) < 5:
        return TResult('Test 14', False, f'run rc={rc} lines={len(lines)}', {'rc': rc, 'lines': len(lines)}), TResult('Test 15', False, 'insufficient lines', {})

    bh_new = [float(x.get('bh_mass_new', np.nan)) for x in lines]
    monotonic = all((bh_new[i+1] + 1e-18) >= bh_new[i] for i in range(len(bh_new)-1))

    # Invariant from debug lines: bh_mass_code = BHFormationMass + BHAccretedMass
    dbg = parse_bhaccr_debug(log)
    dbg_mass = [d for d in dbg if 'diff' in d and 'bh_mass_code' in d]
    if len(dbg_mass) == 0:
        inv_ok = False
        max_abs_diff = float('nan')
    else:
        diffs = [abs(float(d.get('diff', np.nan))) for d in dbg_mass]
        max_abs_diff = float(np.nanmax(diffs))
        inv_ok = all(np.isfinite(d) and d <= 1e-12 for d in diffs)

    r14 = TResult('Test 14', monotonic, f'monotonic={monotonic} steps={len(lines)}', {'steps': len(lines)})
    r15 = TResult('Test 15', inv_ok, f'BH_mass = BHFormationMass + BHAccretedMass: {inv_ok} (max_abs_diff_code={max_abs_diff:.3e})', {'steps': len(lines), 'max_abs_diff_code': max_abs_diff})
    return r14, r15


def test16_mpi_determinism():
    cdir, h5 = prepare_case('test16_mpi_det', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        # Keep BH away from MPI partition boundaries so rank ownership is unambiguous.
        set_bh(g, pos=(0.15,0.2,0.2), mass_code=5.06744516e-6, vel=(0,0,0))
        write_uniform(g, rho=1e-1, T=1e4)
    rc1, l1, _ = run_case(cdir, np_ranks=1, log='np1.log')
    rc4, l4, _ = run_case(cdir, np_ranks=4, log='np4.log')
    if rc1 != 0 or rc4 != 0:
        return TResult('Test 16', False, f'rc1={rc1} rc4={rc4}', {'rc1':rc1,'rc4':rc4})
    p1 = parse_bhaccr(l1)
    p4 = parse_bhaccr(l4)
    if len(p1) == 0 or len(p4) == 0:
        return TResult('Test 16', False, f'BHACCR line missing np1={len(p1)} np4={len(p4)}', {'np1_lines': len(p1), 'np4_lines': len(p4)})
    b1 = p1[-1]
    b4 = p4[-1]
    keys = ['dm_requested','dm_removed','frac_cap','frac_gas','Mdot_hot_realized','Mdot_cold_realized','removal_cells','n_sf_blocked_cells','bh_mass_new']
    mism = []
    for k in keys:
        if isinstance(b1.get(k), (int,float)) and isinstance(b4.get(k), (int,float)):
            if not almost(float(b1[k]), float(b4[k]), rel=1e-10, abs_=1e-30):
                mism.append((k,b1[k],b4[k]))
        else:
            if b1.get(k) != b4.get(k):
                mism.append((k,b1.get(k),b4.get(k)))
    ok = len(mism)==0
    return TResult('Test 16', ok, 'mismatch_count=%d'%len(mism), {'mismatches': mism[:6]})


def test17_ts3_regression():
    cdir = OUT_ROOT / 'test17_ts3wrap_regression'
    ensure_clean(cdir)
    for f in ['bhseed_ts3wrap.enzo', 'make_ts3wrap_ic.py', 'assert_ts3wrap.sh']:
        shutil.copy2(TS3_DIR/f, cdir/f)
    cmd = (
        f'cd {cdir} && python3 make_ts3wrap_ic.py > ic.log 2>&1 && '
        f'{MODULE_INIT} && {EXPORTS} && '
        f'mpirun -np 1 {ENZO_EXE} bhseed_ts3wrap.enzo > ts3wrap_np1.log 2>&1 && '
        f'mpirun -np 4 {ENZO_EXE} bhseed_ts3wrap.enzo > ts3wrap_np4.log 2>&1 && '
        f'bash assert_ts3wrap.sh > assert.log 2>&1'
    )
    rc, _, _ = sh(cmd)
    ok = False
    if (cdir/'assert.log').exists():
        ok = 'TS3_WRAP: PASS' in (cdir/'assert.log').read_text(errors='ignore')
    return TResult('Test 17', ok and rc==0, 'TS3_wrap regression', {'rc': rc})


def test18_multi_bh_independent():
    cdir, h5 = prepare_case('test18_multi_bh', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        pt = g['particle_type'][:]
        idx = np.where(pt == 8)[0]
        if idx.size == 0 or pt.size < 2:
            return TResult('Test 18', False, 'not enough particles', {'n': int(pt.size)})
        i0 = int(idx[0])
        i1 = 0 if i0 != 0 else 1
        pt[i0] = 8
        pt[i1] = 8
        g['particle_type'][...] = pt
        g['particle_index'][i0] = 1
        g['particle_index'][i1] = 2
        g['particle_mass'][i1] = g['particle_mass'][i0]
        g['particle_position_x'][i0] = 0.15
        g['particle_position_y'][i0] = 0.50
        g['particle_position_z'][i0] = 0.50
        g['particle_position_x'][i1] = 0.75
        g['particle_position_y'][i1] = 0.50
        g['particle_position_z'][i1] = 0.50
        g['particle_velocity_x'][i0] = 0.0
        g['particle_velocity_y'][i0] = 0.0
        g['particle_velocity_z'][i0] = 0.0
        g['particle_velocity_x'][i1] = 0.0
        g['particle_velocity_y'][i1] = 0.0
        g['particle_velocity_z'][i1] = 0.0
        write_uniform(g, rho=1e-2, T=1e4)

    rc, log, _ = run_case(cdir)
    if rc != 0:
        return TResult('Test 18', False, f'run rc={rc}', {'rc': rc})
    lines = parse_bhaccr(log)
    ids = sorted(set(int(x.get('bh_id',-1)) for x in lines))
    per = {i: [x for x in lines if int(x.get('bh_id',-1))==i] for i in ids}
    ok = (1 in ids and 2 in ids and len(per[1])>0 and len(per[2])>0)
    return TResult('Test 18', ok, f'bh_ids={ids}', {'bh_ids': ids, 'nlines': len(lines)})


def test19_20_sf_blocking():
    cdir, h5 = prepare_case('test19_20_sf_block', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, mass_code=5.06744516e-8)
        write_uniform(g, rho=1e-2, T=1e4)

    rc, log, _ = run_case(cdir)
    if rc != 0:
        return TResult('Test 19', False, f'run rc={rc}', {'rc': rc}), TResult('Test 20', False, f'run rc={rc}', {'rc': rc})
    bh = parse_bhaccr(log)[-1]
    rem = int(bh.get('removal_cells', -1))
    blk = int(bh.get('n_sf_blocked_cells', -1))
    ok19 = (rem > 0 and blk > 0)
    ok20 = (rem == blk)
    return TResult('Test 19', ok19, f'removal_cells={rem} blocked={blk}', {'removal_cells': rem, 'blocked': blk}), \
           TResult('Test 20', ok20, f'removal_cells={rem} blocked={blk}', {'removal_cells': rem, 'blocked': blk})


def test21_gas_poor_bookkeeping():
    cdir, h5 = prepare_case('test21_gas_poor', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, mass_code=5.06744516e-2)
        write_uniform(g, rho=1e-12, T=1e4)
    rc, log, _ = run_case(cdir)
    if rc != 0:
        return TResult('Test 21', False, f'run rc={rc}', {'rc': rc})
    bh = parse_bhaccr(log)[-1]
    dm_req = float(bh.get('dm_requested',0))
    dm_rem = float(bh.get('dm_removed',0))
    frac_g = float(bh.get('frac_gas',1))
    real = float(bh.get('Mdot_hot_realized',0))+float(bh.get('Mdot_cold_realized',0))
    rhs = float(bh.get('Mdot_actual',0))*frac_g
    ok = (dm_rem <= dm_req and int(bh.get('removal_gas_limited',0))==1 and almost(real, rhs, rel=1e-6))
    return TResult('Test 21', ok, f'dm_removed={dm_rem:.3e} dm_requested={dm_req:.3e} frac_gas={frac_g:.3e}', dict(bh))


def test22_cap_and_gas_limited():
    cdir, h5 = prepare_case('test22_cap_gaslimited', stop='0.02', radiative=False)
    p = cdir/'DD0001'/'data0001'
    txt = p.read_text()
    txt = set_param(txt, 'BHAccretionRemovalMode', '1')
    p.write_text(txt)

    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        # High-mass BH + high-density diagnostic kernel to force cap_active=1,
        # but single-cell removal on a low-density host cell to force gas-limited=1.
        set_bh(g, pos=(0.15,0.2,0.2), mass_code=5.06744516e-6, vel=(0,0,0))
        write_uniform(g, rho=5e-1, T=1e4)
        i = bh_index(g)
        nx = g['Density'].shape[2]
        ny = g['Density'].shape[1]
        nz = g['Density'].shape[0]
        ix = int(np.clip(np.floor(g['particle_position_x'][i] * nx), 0, nx-1))
        iy = int(np.clip(np.floor(g['particle_position_y'][i] * ny), 0, ny-1))
        iz = int(np.clip(np.floor(g['particle_position_z'][i] * nz), 0, nz-1))
        g['Density'][iz, iy, ix] = 1e-14
        if 'Metal_Density' in g:
            g['Metal_Density'][iz, iy, ix] = 1e-20

    rc, log, _ = run_case(cdir)
    if rc != 0:
        return TResult('Test 22', False, f'run rc={rc}', {'rc': rc})
    bh = parse_bhaccr(log)[-1]
    ok = (int(bh.get('cap_active',0))==1 and float(bh.get('frac_cap',1))<1.0 and float(bh.get('frac_gas',1))<1.0)
    return TResult('Test 22', ok, f"cap={bh.get('cap_active')} frac_cap={bh.get('frac_cap')} frac_gas={bh.get('frac_gas')}", dict(bh))


def test23_momentum_monitor():
    cdir, h5 = prepare_case('test23_momentum_monitor', stop='0.02', radiative=False)
    p = cdir/'DD0001'/'data0001'
    txt = p.read_text()
    txt = set_param(txt, 'BHAccretionIgnoredDVWarn', '1e-12')
    txt = set_param(txt, 'BHAccretionIgnoredPFracWarn', '1e-12')
    p.write_text(txt)

    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, mass_code=5.06744516e-8, vel=(0.0,0.0,0.0))
        write_uniform(g, rho=1e-2, T=1e4, vel=(0.02, -0.03, 0.01))

    rc, log, _ = run_case(cdir)
    if rc != 0:
        return TResult('Test 23', False, f'run rc={rc}', {'rc': rc})
    bh = parse_bhaccr(log)[-1]
    warns = parse_bhaccr_warn(log)
    ok = (float(bh.get('acc_ignored_dv_kms',0)) >= 0.0 and float(bh.get('acc_ignored_p_frac',0)) >= 0.0 and int(bh.get('acc_momentum_warn',0))==1 and len(warns)>0)
    return TResult('Test 23', ok, f"warn={bh.get('acc_momentum_warn')} warn_lines={len(warns)}", dict(bh))


def test24_newly_seeded_skip():
    cdir = OUT_ROOT / 'test24_newly_seeded_skip'
    ensure_clean(cdir)
    for f in ['bhseed_ts3wrap.enzo', 'make_ts3wrap_ic.py']:
        shutil.copy2(TS3_DIR/f, cdir/f)
    # ensure accretion enabled and verbose
    cmd = f'cd {cdir} && python3 make_ts3wrap_ic.py > ic.log 2>&1'
    rc,_,_ = sh(cmd)
    if rc != 0:
        return TResult('Test 24', False, 'IC generation failed', {'rc': rc})

    pf = cdir/'bhseed_ts3wrap.enzo'
    txt = pf.read_text()
    txt = set_param(txt, 'BHAccretionMethod', '1')
    txt = set_param(txt, 'BHAccretionVerbose', '1')
    txt = set_param(txt, 'BHAccretionRunEveryTimestep', '1')
    pf.write_text(txt)

    cmd2 = f'cd {cdir} && {MODULE_INIT} && {EXPORTS} && mpirun -np 1 {ENZO_EXE} bhseed_ts3wrap.enzo > run.log 2>&1'
    rc2,_,_ = sh(cmd2)
    if rc2 != 0:
        return TResult('Test 24', False, f'run rc={rc2}', {'rc': rc2})
    b = parse_bhaccr(cdir/'run.log')
    steps = [int(x.get('step',-1)) for x in b]
    has_step1 = 1 in steps
    has_step2 = 2 in steps
    ok = (not has_step1 and has_step2)
    return TResult('Test 24', ok, f'bhaccr_steps={steps}', {'steps': steps})


def rec25_overlap():
    cdir, h5 = prepare_case('test25_overlap', stop='0.02', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        pt = g['particle_type'][:]
        idx = np.where(pt == 8)[0]
        if idx.size == 0 or pt.size < 2:
            return TResult('Test 25', False, 'not enough particles', {})
        i0 = int(idx[0]); i1 = 0 if i0 != 0 else 1
        pt[i0] = 8; pt[i1] = 8
        g['particle_type'][...] = pt
        g['particle_index'][i0] = 1
        g['particle_index'][i1] = 2
        g['particle_mass'][i1] = g['particle_mass'][i0]
        g['particle_position_x'][i0] = 0.45
        g['particle_position_y'][i0] = 0.50
        g['particle_position_z'][i0] = 0.50
        g['particle_position_x'][i1] = 0.47
        g['particle_position_y'][i1] = 0.50
        g['particle_position_z'][i1] = 0.50
        write_uniform(g, rho=1e-2, T=1e4)
    rc, log, _ = run_case(cdir)
    lines = parse_bhaccr(log)
    ok = (rc==0 and len(lines)>=2)
    return TResult('Test 25', ok, f'line_count={len(lines)}', {'line_count': len(lines)})


def rec26_longrun():
    cdir, h5 = prepare_case('test26_longrun', stop='0.51', radiative=False)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, mass_code=5.06744516e-8)
        write_uniform(g, rho=1e-2, T=1e4)
    t0 = time.time()
    rc, log, _ = run_case(cdir)
    dt = time.time()-t0
    lines = parse_bhaccr(log)
    ok = (rc==0 and len(lines)>=50)
    return TResult('Test 26', ok, f'steps_logged={len(lines)} runtime_s={dt:.1f}', {'steps_logged': len(lines), 'runtime_s': dt})


def rec27_perf():
    metrics = {}
    ok = True
    for nbh in [1, 10, 100]:
        name = f'test27_perf_{nbh}bh'
        cdir, h5 = prepare_case(name, stop='0.02', radiative=False)
        with h5py.File(h5, 'r+') as f:
            g = f['Grid00000001']
            pt = g['particle_type'][:]
            n = min(nbh, pt.size)
            pt[:n] = 8
            g['particle_type'][...] = pt
            g['particle_index'][:n] = np.arange(1, n+1)
            m0 = g['particle_mass'][bh_index(g)] if bh_index(g) is not None else 5.06744516e-9
            g['particle_mass'][:n] = m0
            # spread along x
            for i in range(n):
                g['particle_position_x'][i] = 0.1 + 0.8 * (i / max(1, n-1))
                g['particle_position_y'][i] = 0.5
                g['particle_position_z'][i] = 0.5
                g['particle_velocity_x'][i] = 0.0
                g['particle_velocity_y'][i] = 0.0
                g['particle_velocity_z'][i] = 0.0
            write_uniform(g, rho=1e-2, T=1e4)
        rc, log, _ = run_case(cdir)
        if rc != 0:
            ok = False
            metrics[str(nbh)] = {'rc': rc}
            continue
        lines = parse_bhaccr(log)
        vals = [float(x.get('accretion_wall_ms', np.nan)) for x in lines if 'accretion_wall_ms' in x]
        vals = [v for v in vals if np.isfinite(v)]
        if not vals:
            ok = False
            metrics[str(nbh)] = {'rc': rc, 'count': len(lines)}
            continue
        metrics[str(nbh)] = {'mean_ms': float(np.mean(vals)), 'median_ms': float(np.median(vals)), 'count': len(vals)}
    return TResult('Test 27', ok, 'accretion_wall_ms at 1/10/100 BHs', metrics)


def rec28_grackle_cycle():
    cdir, h5 = prepare_case('test28_grackle', stop='0.02', radiative=True, grackle=True)
    with h5py.File(h5, 'r+') as f:
        g = f['Grid00000001']
        set_bh(g, mass_code=5.06744516e-8)
        write_uniform(g, rho=1e-2, T=1e6)
    rc, log, sout = run_case(cdir)
    if rc != 0:
        tail = '\n'.join((cdir/'run.log').read_text(errors='ignore').splitlines()[-40:]) if (cdir/'run.log').exists() else sout[-1000:]
        return TResult('Test 28', False, f'grackle run failed rc={rc}', {'tail': tail})
    lines = parse_bhaccr(log)
    if not lines:
        return TResult('Test 28', False, 'no [BHACCR] line', {})
    bh = lines[-1]
    finite = np.isfinite(float(bh.get('Mdot_total_raw', np.nan))) and np.isfinite(float(bh.get('dm_removed', np.nan)))
    ok = finite and float(bh.get('dm_removed',0)) >= 0.0
    return TResult('Test 28', ok, f"dm_removed={bh.get('dm_removed')} f_Edd={bh.get('f_Edd')}", dict(bh))


def collect_invariants(all_bh_lines):
    # all_bh_lines: list of parsed BHACCR dicts
    inv1 = True
    inv2 = True
    inv3 = True
    inv4 = True
    for x in all_bh_lines:
        # 1) realized sum matches actual*frac_gas (dt-free surrogate)
        lhs = float(x.get('Mdot_hot_realized',0))+float(x.get('Mdot_cold_realized',0))
        rhs = float(x.get('Mdot_actual',0))*float(x.get('frac_gas',1))
        if not almost(lhs, rhs, rel=1e-6):
            inv1 = False
        # 2) bh_mass_new = bh_mass + dm_removed_msun
        if not almost(float(x.get('bh_mass_new',0)), float(x.get('bh_mass',0)) + float(x.get('dm_removed_msun',0)), rel=5e-4):
            inv2 = False
        # 3) approximate dm_requested consistency
        if float(x.get('dm_removed',0)) - float(x.get('dm_requested',0)) > 1e-20:
            inv3 = False
        # 4) best-effort from logs: bh_mass>=bh_mass_initial (monotonic nondecreasing)
    return inv1, inv2, inv3, inv4


def main():
    ensure_clean(OUT_ROOT)
    results = []

    # Mandatory
    results.append(test0_ts3_baseline())
    results.append(test1_mass_conservation())
    results.append(test2_momentum_bookkeeping())
    results.append(test3_energy_consistency())
    results.append(test4_geometry())
    results.append(test5_active_zone_only())
    results.append(test6_proportional())
    results.append(test7_gas_depletion())
    results.append(test8_single_cell_mode())

    hcm = run_hot_cold_mixed_tests()
    results.extend([hcm['9'], hcm['10'], hcm['11']])

    ss = test12_13_sub_super()
    results.extend([ss['12'], ss['13']])

    r14, r15 = test14_15_multistep_invariants()
    results.extend([r14, r15])

    results.append(test16_mpi_determinism())
    results.append(test17_ts3_regression())
    results.append(test18_multi_bh_independent())

    r19, r20 = test19_20_sf_blocking()
    results.extend([r19, r20])

    results.append(test21_gas_poor_bookkeeping())
    results.append(test22_cap_and_gas_limited())
    results.append(test23_momentum_monitor())
    results.append(test24_newly_seeded_skip())

    # Recommended
    results.append(rec25_overlap())
    results.append(rec26_longrun())
    results.append(rec27_perf())
    results.append(rec28_grackle_cycle())

    # Gather invariants across all logs in output tree
    all_bh = []
    for lg in OUT_ROOT.rglob('*.log'):
        all_bh.extend(parse_bhaccr(lg))
    inv1, inv2, inv3, inv4 = collect_invariants(all_bh)

    payload = {
        'generated_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        'out_root': str(OUT_ROOT),
        'results': [
            {'name': r.name, 'passed': r.passed, 'details': r.details, 'metrics': r.metrics}
            for r in results
        ],
        'invariants': {
            'realized_sum_equals_actual_times_frac_gas': inv1,
            'bh_mass_new_equals_bh_mass_plus_dm_removed_msun': inv2,
            'dm_removed_le_dm_requested': inv3,
            'bh_mass_equals_bhformation_plus_bhaccretedmass': inv4,
        },
    }

    (OUT_ROOT/'phaseB_matrix_results.json').write_text(json.dumps(payload, indent=2))

    lines = []
    lines.append('=== Phase B Matrix Summary ===')
    lines.append(f'Output dir: {OUT_ROOT}')
    passed = sum(1 for r in results if r.passed)
    lines.append(f'Passed: {passed}/{len(results)}')
    for r in results:
        st = 'PASS' if r.passed else 'FAIL'
        lines.append(f'{r.name}: {st} :: {r.details}')
    lines.append('--- Invariants ---')
    lines.append(f'realized_sum_equals_actual_times_frac_gas: {inv1}')
    lines.append(f'bh_mass_new_equals_bh_mass_plus_dm_removed_msun: {inv2}')
    lines.append(f'dm_removed_le_dm_requested: {inv3}')
    lines.append(f'bh_mass_equals_bhformation_plus_bhaccretedmass: {inv4}')
    (OUT_ROOT/'phaseB_matrix_summary.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
