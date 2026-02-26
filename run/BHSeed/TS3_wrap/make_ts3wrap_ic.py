#!/usr/bin/env python3
import numpy as np
import h5py
import os

NACTIVE = 50
IBUFF = 3
NTOTAL = NACTIVE + 2 * IBUFF

LengthUnits  = 3.0857e23
DensityUnits = 1.67e-24
TimeUnits    = 3.1557e13
kpc_cm       = 3.0857e21
BoxKpc       = LengthUnits / kpc_cm
CellKpc      = BoxKpc / NACTIVE

OD_THRESH    = 5.0
METAL_THRESH = 1e-4
TEMP_THRESH  = 1e4
EXCL_RADIUS  = 30.0

kboltz  = 1.3807e-16
mh      = 1.6726e-24
gamma   = 5.0 / 3.0
mu      = 0.6

bg_density = 1.0
bg_temp    = 2.0e4
bg_metal   = 0.002

cellA_density = 10.0
cellB_density = 9.0
qual_temp     = 500.0
qual_metal    = 1e-5

iA, jA, kA = 0, 24, 24
iB, jB, kB = 49, 24, 24

def temp_to_energy(T):
    energy_cgs = kboltz * T / ((gamma - 1.0) * mu * mh)
    vel_unit = LengthUnits / TimeUnits
    return energy_cgs / (vel_unit ** 2)

def create_ic(outdir):
    os.makedirs(outdir, exist_ok=True)
    fname = os.path.join(outdir, "data0000")
    shape = (NTOTAL, NTOTAL, NTOTAL)

    density = np.full(shape, bg_density, dtype=np.float64)
    density[kA+IBUFF, jA+IBUFF, iA+IBUFF] = cellA_density
    density[kB+IBUFF, jB+IBUFF, iB+IBUFF] = cellB_density

    velocity = np.zeros(shape, dtype=np.float64)

    bg_energy = temp_to_energy(bg_temp)
    q_energy = temp_to_energy(qual_temp)
    energy = np.full(shape, bg_energy, dtype=np.float64)
    energy[kA+IBUFF, jA+IBUFF, iA+IBUFF] = q_energy
    energy[kB+IBUFF, jB+IBUFF, iB+IBUFF] = q_energy
    total_energy = energy.copy()

    metal_density = np.full(shape, bg_metal, dtype=np.float64) * density
    metal_density[kA+IBUFF, jA+IBUFF, iA+IBUFF] = qual_metal * cellA_density
    metal_density[kB+IBUFF, jB+IBUFF, iB+IBUFF] = qual_metal * cellB_density

    with h5py.File(fname, "w") as f:
        f.create_dataset("Density", data=density)
        f.create_dataset("x-velocity", data=velocity)
        f.create_dataset("y-velocity", data=velocity)
        f.create_dataset("z-velocity", data=velocity)
        f.create_dataset("TotalEnergy", data=total_energy)
        f.create_dataset("GasEnergy", data=energy)
        f.create_dataset("Metal_Density", data=metal_density)

    print(f"Created {fname}")
    print(f"  Shape: {shape}")
    print(f"  Cell A (ix={iA}, iy={jA}, iz={kA}) -> HDF5 [{kA+IBUFF},{jA+IBUFF},{iA+IBUFF}]")
    print(f"    density={cellA_density}, temp={qual_temp} K, metal={qual_metal}")
    print(f"    code pos = ({(iA+0.5)/NACTIVE:.4f}, {(jA+0.5)/NACTIVE:.4f}, {(kA+0.5)/NACTIVE:.4f})")
    print(f"  Cell B (ix={iB}, iy={jB}, iz={kB}) -> HDF5 [{kB+IBUFF},{jB+IBUFF},{iB+IBUFF}]")
    print(f"    density={cellB_density}, temp={qual_temp} K, metal={qual_metal}")
    print(f"    code pos = ({(iB+0.5)/NACTIVE:.4f}, {(jB+0.5)/NACTIVE:.4f}, {(kB+0.5)/NACTIVE:.4f})")
    print(f"  Background: density={bg_density}, temp={bg_temp} K, metal={bg_metal}")
    print(f"  BoxKpc = {BoxKpc}, CellKpc = {CellKpc}")
    print(f"  Direct distance A->B = {abs(iB - iA) * CellKpc:.1f} kpc")
    print(f"  Wrapped distance A->B = {(NACTIVE - abs(iB - iA)) * CellKpc:.1f} kpc")
    print(f"  ExclusionRadius = {EXCL_RADIUS} kpc")

    verify_ic(fname)

def verify_ic(fname):
    print("\n--- Verification ---")
    with h5py.File(fname, "r") as f:
        d = f["Density"][:]
        e = f["GasEnergy"][:]
        m = f["Metal_Density"][:]
        vel_unit = LengthUnits / TimeUnits

        for label, ix, iy, iz in [("A", iA, jA, kA), ("B", iB, jB, kB)]:
            hz, hy, hx = iz+IBUFF, iy+IBUFF, ix+IBUFF
            rho = d[hz, hy, hx]
            eng = e[hz, hy, hx]
            met = m[hz, hy, hx] / rho if rho > 0 else 0
            T = eng * (vel_unit**2) * (gamma - 1.0) * mu * mh / kboltz
            print(f"  Cell {label}: density={rho:.1f}, T={T:.0f} K, Z={met:.1e}")
            assert rho > OD_THRESH
            assert T < TEMP_THRESH
            assert met < METAL_THRESH

        mid = NACTIVE // 2 + IBUFF
        rho_bg = d[mid, mid, mid]
        eng_bg = e[mid, mid, mid]
        T_bg = eng_bg * (vel_unit**2) * (gamma - 1.0) * mu * mh / kboltz
        met_bg = m[mid, mid, mid] / rho_bg
        print(f"  Background (center): density={rho_bg:.1f}, T={T_bg:.0f} K, Z={met_bg:.4f}")
        assert (rho_bg < OD_THRESH) or (T_bg > TEMP_THRESH) or (met_bg > METAL_THRESH)

    print("  All verifications PASSED.")

if __name__ == "__main__":
    create_ic("DD0000")
