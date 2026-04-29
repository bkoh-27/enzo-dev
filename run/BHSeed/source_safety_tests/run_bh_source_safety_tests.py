#!/usr/bin/env python3
"""Fast source-level regressions for BH seed defaults and BH source safety."""

import re
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def require(pattern, text, message):
    if not re.search(pattern, text, re.MULTILINE):
        raise AssertionError(message)


def test_defaults_are_opt_in():
    defaults = read("src/enzo/SetDefaultGlobalValues.C")
    require(
        r"^\s*BHSeedingMethod\s*=\s*0\s*;",
        defaults,
        "BHSeedingMethod must default to 0 so missing seed parameters do not run seeding.",
    )
    require(
        r"^\s*BHAccretionMethod\s*=\s*0\s*;",
        defaults,
        "BHAccretionMethod must default to 0 so old parameter files do not run the new active source path.",
    )


def test_write_parameter_file_persists_disabled_state():
    writer = read("src/enzo/WriteParameterFile.C")
    require(
        r'fprintf\(fptr,\s*"BHSeedingMethod\s*=',
        writer,
        "WriteParameterFile must emit BHSeedingMethod, including the disabled state.",
    )
    require(
        r'fprintf\(fptr,\s*"BHAccretionMethod\s*=',
        writer,
        "WriteParameterFile must emit BHAccretionMethod, including the disabled state.",
    )


def test_accretion_safety_helper_compiles_and_passes():
    src = ROOT / "run/BHSeed/source_safety_tests/bh_accretion_safety_unit.C"
    with tempfile.TemporaryDirectory(prefix="bh_accretion_safety_") as td:
        exe = Path(td) / "bh_accretion_safety_unit"
        cmd = ["g++", "-std=c++11", "-Wall", "-Wextra", "-Werror", str(src), "-o", str(exe)]
        subprocess.run(cmd, cwd=ROOT, check=True)
        subprocess.run([str(exe)], cwd=ROOT, check=True)


def test_large_kernels_are_rejected_or_diagnostic_only():
    acc = read("src/enzo/Grid_BHAccretionHandler.C")
    require(
        r"accretion_source_rejected\s*=\s*kernel_truncated",
        acc,
        "BHAccretion must reject gas removal when the diagnostic kernel is truncated.",
    )
    require(
        r"kernel_cells_requested=[\s\S]*kernel_cells_valid=[\s\S]*kernel_cells_active=",
        acc,
        "BHAccretion source diagnostics must report requested, valid, and active kernel cells.",
    )

    repos = read("src/enzo/Grid_BHRepositionHandler.C")
    require(
        r"active_reposition_rejected\s*=\s*\(\s*BHRepositionMethod\s*>\s*0\s*&&\s*search_kernel_truncated\s*\)",
        repos,
        "Active BH repositioning must be rejected when the search kernel is truncated.",
    )
    require(
        r"BHRepositionMethod=0[\s\S]*diagnostics only",
        repos,
        "BHRepositionMethod=0 with verbose logging must be documented in the runtime log.",
    )


def test_legacy_mbh_mass_removal_clamps_energy():
    legacy = read("src/enzo/Grid_SubtractAccretedMassFromSphere.C")
    require(
        r'#include\s+"BHAccretionSafety.h"',
        legacy,
        "Legacy MBH mass removal must use the shared accretion safety helpers.",
    )
    require(
        r"energy_clamped_after_mass_removal=1",
        legacy,
        "Legacy MBH mass removal must log when it clamps energy after gas removal.",
    )


def main():
    tests = [
        test_defaults_are_opt_in,
        test_write_parameter_file_persists_disabled_state,
        test_accretion_safety_helper_compiles_and_passes,
        test_large_kernels_are_rejected_or_diagnostic_only,
        test_legacy_mbh_mass_removal_clamps_energy,
    ]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL {exc}", file=sys.stderr)
        sys.exit(1)
