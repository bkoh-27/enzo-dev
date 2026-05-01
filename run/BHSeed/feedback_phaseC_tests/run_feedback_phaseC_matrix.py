#!/usr/bin/env python3
"""Conservative Phase C regression harness for Enzo BH feedback D0a/D0b-lite.

This script is validation-only. It copies TS3_wrap parameter files into an
output directory, applies per-case overrides there, and never edits source
inputs in place.

BH seeding, accretion, and feedback are opt-in in source defaults. The TS3_wrap
base fixture explicitly sets BHSeedingMethod=1, and this harness explicitly
sets BHAccretionMethod=1 and BHFeedbackMethod=1 in cases that require them.

Environment overrides:
  BH_FEEDBACK_PHASEC_OUT_ROOT   Output root. Default:
                                /tmp/bh_feedback_phaseC_tests/<timestamp>
  ENZO_EXE                      Enzo executable. Default: src/enzo/enzo.exe
  MPIRUN                        MPI launcher. Default:
                                /home/bkoh/miniconda3/envs/enzo/bin/mpirun,
                                falling back to mpirun on PATH
  BH_FEEDBACK_PHASEC_TS3_PARAM  TS3 base parameter file. Default:
                                run/BHSeed/TS3_wrap/bhseed_ts3wrap.enzo
  BH_FEEDBACK_PHASEC_TIMEOUT    Per-case timeout in seconds. Default: 300
  BH_FEEDBACK_PHASEC_RUN_PHASEB Set to 1 to run optional T12-lite.
  BH_FEEDBACK_PHASEB_TIMEOUT    Optional Phase B timeout in seconds. Default: 900
  NP                            MPI rank count for T10 np=1 vs np=NP. Default: 4
  H5LS / H5DUMP                 HDF5 CLI fallbacks if h5py is unavailable.
"""

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
TS3_PARAM = Path(os.environ.get(
    "BH_FEEDBACK_PHASEC_TS3_PARAM",
    ROOT / "run/BHSeed/TS3_wrap/bhseed_ts3wrap.enzo",
))
if not TS3_PARAM.is_absolute():
    TS3_PARAM = ROOT / TS3_PARAM

OUT_ROOT = Path(os.environ.get(
    "BH_FEEDBACK_PHASEC_OUT_ROOT",
    Path("/tmp") / "bh_feedback_phaseC_tests" / time.strftime("%Y%m%d_%H%M%S"),
))

ENZO_EXE = Path(os.environ.get("ENZO_EXE", ROOT / "src/enzo/enzo.exe"))
if not ENZO_EXE.is_absolute():
    ENZO_EXE = ROOT / ENZO_EXE

DEFAULT_CONDA_BIN = Path("/home/bkoh/miniconda3/envs/enzo/bin")
DEFAULT_MPIRUN = DEFAULT_CONDA_BIN / "mpirun"
MPIRUN = os.environ.get("MPIRUN")
if MPIRUN is None:
    MPIRUN = str(DEFAULT_MPIRUN if DEFAULT_MPIRUN.exists() else "mpirun")

H5LS = os.environ.get("H5LS")
if H5LS is None:
    default_h5ls = DEFAULT_CONDA_BIN / "h5ls"
    H5LS = str(default_h5ls if default_h5ls.exists() else "h5ls")

H5DUMP = os.environ.get("H5DUMP")
if H5DUMP is None:
    default_h5dump = DEFAULT_CONDA_BIN / "h5dump"
    H5DUMP = str(default_h5dump if default_h5dump.exists() else "h5dump")

CASE_TIMEOUT = int(os.environ.get("BH_FEEDBACK_PHASEC_TIMEOUT", "300"))
PHASEB_TIMEOUT = int(os.environ.get("BH_FEEDBACK_PHASEB_TIMEOUT", "900"))
T10_NP = int(os.environ.get("NP", "4"))

D0B_BASELINE = "52cf4ceafd036f8f74068797092e064d5459511f"
D0A_COMMIT = "685c7a9c0f5b679a8648464f0011076bb78b79df"
PRESERVED_STASH_SUBSTRING = "Preserve gen_uuid.c outside Phase B merge"

EDD_FACTOR_FATAL = "BHFeedbackEddingtonFactor must be >= 0"
EFFICIENCY_PRODUCT_FATAL = (
    "Product BHAccretionRadiativeEfficiency * "
    "BHFeedbackThermalEfficiency exceeds 1.0"
)
DUAL_MBH_BHACCRETION_FATAL = (
    "MBHAccretion=1 and BHAccretionMethod>0 are mutually exclusive"
)

HDF5_EXPECTED_LABELS = (
    "bh_cumul_reservoir_in",
    "bh_cumul_reservoir_out",
    "bhaccr_last_mdot_realized",
)
HDF5_STALE_LABELS = (
    "particle_attribute_20",
    "particle_attribute_21",
    "particle_attribute_22",
)

PROTECTED_PATHS = (
    "src/enzo/",
    "src/enzo/uuid/gen_uuid.c",
    ".codex",
    "handover_architect_role.md",
    "plan_feedback_phaseC_v2.md",
    "plan_feedback_phaseC_v3.md",
    "run/BHSeed/feedback_phaseB_tests/",
    "run/BHSeed/TS3_wrap/",
    "src/enzo/Make.mach.grammar",
)

F1_BURST_THRESHOLD_CGS = 1e40
D0B_F1_FEEDBACK_KERNEL_RADIUS = "8.0"
D0B_CONSERVATION_RTOL = 1e-5
TINY_ABS_TOL = 1e-30
RESTART_FLOAT_RTOL = 1e-6
EVOLVED_BH_FLOAT_RTOL = 2e-5
RESTART_FLOAT_ATOL = 1e-12
T10_SHORT_RTOL = 1e-10
T10_ABS_TOL = 1e-30
D0C_ACC_BONDI_RTOL = 1e-6
D0C_ACC_STRICT_RTOL = 1e-8
D0C_ACC_STRICT_ATOL = 1e-30
D0C_F6_MIN_MISMATCH = 0.9

CGS_GRAV_CONST = 6.67428e-8
CGS_SOLAR_MASS = 1.9891e33
CGS_YR_S = 3.1556952e7

T7_RESTART_FIELDS = (
    "ParticleMass",
    "ParticlePosition",
    "ParticleVelocity",
    "bhseed_channel",
    "bhseed_redshift",
    "bhseed_patch_mass",
    "bhseed_patch_metallicity",
    "bhseed_patch_density_peak",
    "bhseed_kernel_complete",
    "bhseed_host_dm_density",
    "bhseed_accept_rank",
    "bhaccr_accreted_mass",
    "bhaccr_reservoir_mass",
    "bhaccr_last_accretion_redshift",
    "bhaccr_last_eddington_ratio",
    "bh_formation_mass",
    "bhaccr_last_mdot_actual",
    "bhfdbk_energy_reservoir",
    "bhfdbk_last_feedback_redshift",
    "bh_cumul_reservoir_in",
    "bh_cumul_reservoir_out",
)

T7_CATEGORICAL_FIELDS = {
    "bhseed_channel",
    "bhseed_kernel_complete",
    "bhseed_accept_rank",
}

T7_EVOLVED_BH_FLOAT_FIELDS = {
    "ParticleMass",
    "bhaccr_accreted_mass",
    "bhaccr_last_eddington_ratio",
    "bhaccr_last_mdot_actual",
    "bhfdbk_energy_reservoir",
    "bh_cumul_reservoir_in",
}

BH_PARTICLE_TYPES = {6, 8}
D0C_SEED_EXPECTED_PARTICLE_TYPE = 8
D0C_SEED_EXPECTED_MASS = 4.054e-14
D0C_SEED_MASS_RTOL = 1e-3
D0C_SEED_EXPECTED_POSITION = (0.01, 0.49, 0.49)
D0C_SEED_POSITION_ATOL = 1e-6
D0C_SEED_DETERMINISTIC_POSITION_ATOL = 1e-15
D0C_SEED_DETERMINISTIC_MASS_ATOL = 1e-20
D0C_SEED_MPI_POSITION_ATOL = 1e-12
D0C_REPOS_POSITION_ATOL = 1e-10
D0C_REPOS_ZERO_ATOL = 1e-12
D0C_REPOS_COUNTER_MAX_SAFE = 1000000000
D0C_REPOS_R4_INITIAL_POSITION = (0.15, 0.49, 0.49)
D0C_REPOS_R4_TARGET_POSITION = (0.19, 0.49, 0.49)
D0C_REPOS_R4_POSITION_ATOL = 1e-10
D0C_REPOS_R4_DISPLACEMENT_ATOL = 1e-8
D0C_REPOS_R4_SEARCH_RADIUS_KPC = 6.0
D0C_REPOS_R4_BACKGROUND_DENSITY = 1e-6
D0C_REPOS_R4_TARGET_DENSITY = 1e-3

T10_BHACCR_NUMERIC_FIELDS = (
    ("lambda_edd", ("lambda_edd",)),
    ("edd_factor", ("edd_factor",)),
    ("mdot_total_raw", ("mdot_total_raw", "Mdot_total_raw")),
    ("mdot_actual", ("mdot_actual", "Mdot_actual", "mdot_total_capped", "Mdot_total_capped")),
    ("mdot_edd", ("mdot_edd", "Mdot_Edd")),
    ("frac_cap", ("frac_cap",)),
    ("frac_gas", ("frac_gas",)),
    ("f_hot", ("f_hot",)),
    ("f_cold", ("f_cold",)),
    ("alpha_boost", ("alpha_boost",)),
    ("f_AM", ("f_AM", "f_am")),
    ("dm_requested", ("dm_requested",)),
    ("dm_removed", ("dm_removed",)),
)

T10_BHACCR_CATEGORICAL_FIELDS = (
    "cap_active",
    "removal_gas_limited",
)

T10_BHFDBK_NUMERIC_FIELDS = (
    ("f_Edd", ("f_Edd",)),
    ("E_requested", ("E_requested",)),
    ("reservoir_before", ("reservoir_before",)),
    ("reservoir_after", ("reservoir_after", "reservoir_after_accum")),
    ("reservoir_final", ("reservoir_final",)),
    ("E_deposited", ("E_deposited",)),
    ("cumul_reservoir_in_cgs", ("cumul_reservoir_in_cgs",)),
    ("cumul_reservoir_out_cgs", ("cumul_reservoir_out_cgs",)),
    ("conservation_residual_cgs", ("conservation_residual_cgs",)),
)

T10_BHFDBK_CATEGORICAL_FIELDS = (
    "feedback_mode",
    "burst_diag",
    "burst_occurred",
    "kernel_gas_zero",
)

D2_1_IO_LABEL_FILES = (
    "Grid_ReadGrid.C",
    "Grid_WriteGrid.C",
    "Grid_Group_ReadGrid.C",
    "Grid_Group_WriteGrid.C",
    "New_Grid_ReadGrid.C",
    "New_Grid_WriteGrid.C",
    "AMRH5writer.C",
    "Grid_ConvertToNumpy.C",
    "Grid_WriteCube.C",
    "Grid_WriteGridX.C",
)

D2_1_BHACCR_REQUIRED_FIELDS = (
    "Mdot_actual",
    "Mdot_realized",
    "mdot_actual_code",
    "mdot_realized_code",
    "dm_requested",
    "dm_removed",
    "frac_gas",
    "realized_over_requested",
    "realized_attr_written",
)

D2_1_BHFDBK_REQUIRED_FIELDS = (
    "feedback_mdot_basis",
    "mdot_actual_code",
    "mdot_realized_code",
    "mdot_actual_msunyr",
    "mdot_realized_msunyr",
    "L_feedback",
    "L_feedback_realized_basis",
    "E_requested",
    "E_requested_realized_basis",
    "realized_attr_present",
    "realized_attr_anomalous",
)

D0C_ACC_REQUIRED_FIELDS = (
    "bh_mass",
    "f_hot",
    "f_cold",
    "f_AM",
    "rho_cold_avg",
    "cs_cold_avg",
    "v_rel_cold",
    "Mdot_hot_raw",
    "Mdot_cold_raw",
    "Mdot_total_raw",
    "Mdot_Edd",
    "Mdot_actual",
    "cap_active",
    "frac_cap",
    "frac_gas",
    "dm_requested",
    "dm_removed",
    "removal_gas_limited",
    "lambda_edd",
    "edd_factor",
)


@dataclass
class Result:
    name: str
    status: str
    details: str
    log: str = ""
    metrics: dict = field(default_factory=dict)

    @property
    def failed_required(self):
        return self.status == "FAIL"


@dataclass
class CaseRun:
    name: str
    cdir: Path
    rc: int
    log: Path
    timed_out: bool = False


def run_cmd(args, cwd=ROOT, env=None, timeout=None):
    return subprocess.run(
        [str(a) for a in args],
        cwd=cwd,
        env=env,
        timeout=timeout,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def git_cmd(*args):
    proc = run_cmd(["git", *args], cwd=ROOT)
    return proc.returncode, proc.stdout


def git_output(*args):
    rc, out = git_cmd(*args)
    if rc != 0:
        return f"<git {' '.join(args)} failed rc={rc}>\n{out}"
    return out


def git_is_ancestor(commit, ref="HEAD"):
    rc, _ = git_cmd("merge-base", "--is-ancestor", commit, ref)
    return rc == 0


def executable_available(command):
    path = Path(command)
    if path.is_absolute() or "/" in command:
        return path.exists() and os.access(path, os.X_OK)
    return shutil.which(command) is not None


def pass_result(name, details, metrics=None, log=""):
    return Result(name, "PASS", details, log=log, metrics=metrics or {})


def fail_result(name, details, metrics=None, log=""):
    return Result(name, "FAIL", details, log=log, metrics=metrics or {})


def skip_result(name, details, metrics=None, log=""):
    return Result(name, "SKIP", details, log=log, metrics=metrics or {})


def optional_fail_result(name, details, metrics=None, log=""):
    return Result(name, "OPTIONAL_FAIL", details, log=log, metrics=metrics or {})


def set_param(text, key, value):
    pat = re.compile(rf"^{re.escape(key)}\s*=.*$", re.MULTILINE)
    line = f"{key:<36}= {value}"
    if pat.search(text):
        return pat.sub(line, text)
    return text.rstrip() + "\n" + line + "\n"


def remove_param(text, key):
    pat = re.compile(rf"^{re.escape(key)}\s*=.*\n?", re.MULTILINE)
    return pat.sub("", text)


def make_text(base, overrides, remove_keys=()):
    text = Path(base).read_text()
    for key in remove_keys:
        text = remove_param(text, key)
    for key, value in overrides.items():
        text = set_param(text, key, str(value))
    return text


def baseline_overrides(**overrides):
    # TS3_wrap explicitly sets BHSeedingMethod=1. Accretion and feedback
    # defaults are opt-in, so every Phase C case that needs them enables them
    # here rather than relying on source defaults.
    params = {
        "BHAccretionMethod": 1,
        "BHAccretionRunEveryTimestep": 1,
        "BHAccretionVerbose": 1,
        "BHFeedbackMethod": 1,
        "BHFeedbackModeThreshold": "1e-10",
        "BHFeedbackMinEnergyBurst": "1e50",
        "BHFeedbackVerbose": 1,
    }
    params.update(overrides)
    return params


def d0b_burst_overrides(**overrides):
    # D0b-lite fixture: TS3_wrap's default 1 kpc feedback kernel reaches the
    # burst threshold but can contain no active gas. The smallest probed radius
    # that produces a gas-present deposited burst is 8 physical kpc.
    params = baseline_overrides(
        BHFeedbackMinEnergyBurst="1e40",
        BHFeedbackKernelRadius=D0B_F1_FEEDBACK_KERNEL_RADIUS,
    )
    params.update(overrides)
    return params


def d0b_no_burst_overrides(**overrides):
    params = baseline_overrides(BHFeedbackMinEnergyBurst="1e50")
    params.update(overrides)
    return params


def d0c_acc_base_overrides(**overrides):
    params = baseline_overrides(
        BHSeedingMethod=1,
        BHAccretionMethod=1,
        BHAccretionRunEveryTimestep=1,
        BHAccretionVerbose=1,
        BHFeedbackMethod=1,
        BHFeedbackModeThreshold="1e-10",
        BHFeedbackMinEnergyBurst="1e50",
        BHFeedbackVerbose=1,
        SelfGravity=0,
        BHFeedbackEddingtonFactor=100,
    )
    params.update(overrides)
    return params


def d0c_seed_opt_in_overrides(**overrides):
    params = {
        "BHSeedingMethod": 1,
        "BHSeedRunEveryTimestep": 1,
        "BHSeedVerbose": 1,
        "BHAccretionMethod": 0,
        "BHFeedbackMethod": 0,
    }
    params.update(overrides)
    return params


def d0c_repos_base_overrides(**overrides):
    params = d0c_seed_opt_in_overrides(
        BHRepositionMethod=0,
        BHRepositionVerbose=1,
    )
    params.update(overrides)
    return params


def prepare_case(name, overrides=None, remove_keys=()):
    overrides = overrides or {}
    cdir = OUT_ROOT / name
    if cdir.exists():
        shutil.rmtree(cdir)
    cdir.mkdir(parents=True, exist_ok=True)
    text = make_text(TS3_PARAM, overrides, remove_keys=remove_keys)
    (cdir / "case.enzo").write_text(text)
    return cdir


def enzo_command(param_file="case.enzo", np_ranks=1, restart=False, mpirun=None):
    launcher = mpirun or MPIRUN
    cmd = [launcher, "-n", str(np_ranks), str(ENZO_EXE), "-d"]
    if restart:
        cmd.append("-r")
    cmd.append(str(param_file))
    return cmd


def run_enzo(
        cdir, np_ranks=1, log_name="run.log", param_file="case.enzo",
        restart=False, mpirun=None):
    env = os.environ.copy()
    env["HDF5_DISABLE_VERSION_CHECK"] = "2"
    if DEFAULT_CONDA_BIN.exists():
        env["PATH"] = f"{DEFAULT_CONDA_BIN}:{env.get('PATH', '')}"
    cmd = enzo_command(
        param_file=param_file, np_ranks=np_ranks, restart=restart, mpirun=mpirun)
    log_path = cdir / log_name
    with open(log_path, "w") as log:
        try:
            proc = subprocess.run(
                [str(x) for x in cmd],
                cwd=cdir,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                timeout=CASE_TIMEOUT,
            )
            return CaseRun(cdir.name, cdir, proc.returncode, log_path)
        except subprocess.TimeoutExpired:
            log.write(f"\nTIMEOUT after {CASE_TIMEOUT} seconds\n")
            return CaseRun(cdir.name, cdir, 124, log_path, timed_out=True)


def run_named(name, overrides, remove_keys=()):
    cdir = prepare_case(name, overrides=overrides, remove_keys=remove_keys)
    return run_enzo(cdir)


def parse_kv(line):
    out = {}
    for tok in line.strip().split():
        if "=" not in tok:
            continue
        key, val = tok.split("=", 1)
        if key.startswith("["):
            continue
        val = val.rstrip(",")
        try:
            if re.fullmatch(r"-?\d+", val):
                out[key] = int(val)
            else:
                out[key] = float(val)
        except ValueError:
            out[key] = val
    out["__line"] = line.strip()
    return out


def lines(log, tag):
    if not log.exists():
        return []
    return [ln.strip() for ln in log.read_text(errors="ignore").splitlines()
            if tag in ln]


def rows_for(log, tag, warn_tag=None, min_step=2):
    rows = []
    for line in lines(log, tag):
        if warn_tag and warn_tag in line:
            continue
        row = parse_kv(line)
        step = row.get("step")
        if isinstance(step, int) and step >= min_step:
            rows.append(row)
    return rows


def bhaccr_rows(log, min_step=2):
    return rows_for(log, "[BHACCR]", warn_tag="[BHACCR_WARN]", min_step=min_step)


def bhfdbk_rows(log, min_step=2):
    return rows_for(log, "[BHFDBK]", warn_tag="[BHFDBK_WARN]", min_step=min_step)


def is_finite_number(value):
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def check_field(row, field, predicate, reason):
    if field not in row:
        return f"missing {field}: {row.get('__line', '')}"
    value = row[field]
    if not predicate(value):
        return f"{field}={value!r} {reason}: {row.get('__line', '')}"
    return ""


def lambda_key(row):
    return (row.get("step"), row.get("bh_id"))


def lambda_map(rows):
    return {lambda_key(row): float(row["lambda_edd"])
            for row in rows
            if "lambda_edd" in row and is_finite_number(row.get("lambda_edd"))}


def lambda_comparable(current_rows, baseline_rows, rel_tol=1e-4, abs_tol=1e-12):
    current = lambda_map(current_rows)
    baseline = lambda_map(baseline_rows)
    common = sorted(set(current) & set(baseline))
    if not common:
        return False, "no common (step,bh_id) lambda_edd entries with T1"
    failures = []
    for key in common:
        if not math.isclose(current[key], baseline[key],
                            rel_tol=rel_tol, abs_tol=abs_tol):
            failures.append((key, current[key], baseline[key]))
    if failures:
        return False, f"lambda_edd mismatch examples={failures[:3]}"
    return True, f"{len(common)} lambda_edd entries comparable to T1"


def validate_bhaccr_common(rows, expected_edd_factor=None, cap_expected_zero=True):
    if not rows:
        return ["no [BHACCR] rows with step >= 2"]
    errors = []
    for row in rows:
        err = check_field(
            row, "lambda_edd",
            lambda v: is_finite_number(v) and float(v) > 0.0,
            "must be finite and positive",
        )
        if err:
            errors.append(err)
        if expected_edd_factor is not None:
            err = check_field(
                row, "edd_factor",
                lambda v: is_finite_number(v) and math.isclose(
                    float(v), float(expected_edd_factor),
                    rel_tol=1e-12, abs_tol=1e-12),
                f"must equal {expected_edd_factor}",
            )
            if err:
                errors.append(err)
        if cap_expected_zero and "cap_active" in row and row.get("cap_active") != 0:
            errors.append(
                f"cap_active={row.get('cap_active')} expected 0: "
                f"{row.get('__line', '')}"
            )
    return errors


def thermal_bhfdbk_rows(log):
    return [row for row in bhfdbk_rows(log)
            if row.get("feedback_mode") == "THERMAL"]


def all_thermal_bhfdbk_rows(log):
    return [row for row in bhfdbk_rows(log, min_step=0)
            if row.get("feedback_mode") == "THERMAL"]


def row_value(row, names):
    for name in names:
        if name in row:
            return row[name]
    return None


def row_has_any(row, names):
    return any(name in row for name in names)


def conservation_tolerance(cumul_in):
    if cumul_in > 0.0:
        return D0B_CONSERVATION_RTOL * cumul_in
    return TINY_ABS_TOL


def conservation_check(row):
    for field in ("cumul_reservoir_in_cgs", "conservation_residual_cgs"):
        if field not in row:
            return False, f"missing {field}"
        if not is_finite_number(row[field]):
            return False, f"{field}={row[field]!r} is not finite"
    cumul_in = float(row["cumul_reservoir_in_cgs"])
    residual = float(row["conservation_residual_cgs"])
    tolerance = conservation_tolerance(cumul_in)
    if abs(residual) > tolerance:
        return (
            False,
            f"abs(conservation_residual_cgs)={abs(residual):.8e} exceeds "
            f"tolerance={tolerance:.8e}",
        )
    return True, f"abs_residual={abs(residual):.8e} tolerance={tolerance:.8e}"


def compact_row(row, fields):
    out = {}
    for field in fields:
        if field in row:
            out[field] = row[field]
    out["line"] = row.get("__line", "")
    return out


def validate_t1(case):
    name = "T1 default Phase C smoke"
    if case.timed_out:
        return fail_result(name, f"run timed out after {CASE_TIMEOUT}s", log=str(case.log))
    if case.rc != 0:
        return fail_result(name, f"run exited rc={case.rc}", log=str(case.log))

    errors = validate_bhaccr_common(bhaccr_rows(case.log), expected_edd_factor=1.0)
    thermal_rows = thermal_bhfdbk_rows(case.log)
    if not thermal_rows:
        errors.append("no THERMAL [BHFDBK] rows with step >= 2")
    for row in thermal_rows:
        for field in (
            "cumul_reservoir_in_cgs",
            "cumul_reservoir_out_cgs",
            "conservation_residual_cgs",
        ):
            err = check_field(
                row, field,
                lambda v: is_finite_number(v),
                "must be finite",
            )
            if err:
                errors.append(err)
        for field in ("cumul_reservoir_in_cgs", "cumul_reservoir_out_cgs"):
            if field in row and is_finite_number(row[field]) and float(row[field]) < 0.0:
                errors.append(f"{field}={row[field]} expected non-negative")
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8]}, str(case.log))
    return pass_result(
        name,
        f"rc=0; BHACCR rows={len(bhaccr_rows(case.log))}; "
        f"THERMAL BHFDBK rows={len(thermal_rows)}",
        log=str(case.log),
    )


def validate_edd_factor_case(name, case, expected_edd_factor, baseline_rows):
    if case.timed_out:
        return fail_result(name, f"run timed out after {CASE_TIMEOUT}s", log=str(case.log))
    if case.rc != 0:
        return fail_result(name, f"run exited rc={case.rc}", log=str(case.log))
    rows = bhaccr_rows(case.log)
    errors = validate_bhaccr_common(rows, expected_edd_factor=expected_edd_factor)
    ok, detail = lambda_comparable(rows, baseline_rows)
    if not ok:
        errors.append(detail)
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8]}, str(case.log))
    return pass_result(
        name,
        f"edd_factor={expected_edd_factor}; {detail}",
        {"bhaccr_rows": len(rows)},
        str(case.log),
    )


def rows_by_bh(rows):
    grouped = {}
    for row in rows:
        grouped.setdefault(row.get("bh_id"), []).append(row)
    return grouped


def nondecreasing(values, rel_tol=1e-8, abs_tol=1e-30):
    previous = None
    for value in values:
        if previous is not None:
            tolerance = max(abs_tol, abs(previous) * rel_tol)
            if value + tolerance < previous:
                return False, previous, value
        previous = value
    return True, None, None


def validate_t4(case):
    name = "T4 cumulative reservoir counters"
    if case.rc != 0:
        return fail_result(name, f"T1 run unavailable rc={case.rc}", log=str(case.log))
    errors = []
    thermal_rows = thermal_bhfdbk_rows(case.log)
    if not thermal_rows:
        errors.append("no THERMAL [BHFDBK] rows with step >= 2")
    for bh_id, rows in rows_by_bh(thermal_rows).items():
        for field in ("cumul_reservoir_in_cgs", "cumul_reservoir_out_cgs"):
            values = []
            for row in rows:
                err = check_field(
                    row, field,
                    lambda v: is_finite_number(v) and float(v) >= 0.0,
                    "must be finite and non-negative",
                )
                if err:
                    errors.append(err)
                else:
                    values.append(float(row[field]))
            ok, previous, current = nondecreasing(values)
            if not ok:
                errors.append(
                    f"bh_id={bh_id} {field} decreased from {previous} to {current}"
                )
    inactive_rows = [
        row for row in bhfdbk_rows(case.log)
        if row.get("feedback_mode") == "KINETIC_INACTIVE"
    ]
    for row in inactive_rows:
        for field in ("cumul_reservoir_in_cgs", "cumul_reservoir_out_cgs"):
            err = check_field(
                row, field,
                lambda v: is_finite_number(v) and math.isclose(
                    float(v), 0.0, rel_tol=0.0, abs_tol=1e-30),
                "must remain zero for KINETIC_INACTIVE BH",
            )
            if err:
                errors.append(err)
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8]}, str(case.log))
    return pass_result(
        name,
        f"thermal_bh_count={len(rows_by_bh(thermal_rows))}; "
        f"kinetic_inactive_rows={len(inactive_rows)}",
        log=str(case.log),
    )


def validate_t5(case):
    name = "T5 conservation_residual diagnostic"
    if case.rc != 0:
        return fail_result(name, f"T1 run unavailable rc={case.rc}", log=str(case.log))
    errors = []
    thermal_rows = thermal_bhfdbk_rows(case.log)
    if not thermal_rows:
        errors.append("no THERMAL [BHFDBK] rows with step >= 2")
    for row in thermal_rows:
        for field in ("cumul_reservoir_in_cgs", "conservation_residual_cgs"):
            err = check_field(row, field, is_finite_number, "must be finite")
            if err:
                errors.append(err)
        if "cumul_reservoir_in_cgs" not in row or "conservation_residual_cgs" not in row:
            continue
        cumul_in = float(row["cumul_reservoir_in_cgs"])
        residual = float(row["conservation_residual_cgs"])
        if cumul_in > 0.0:
            if abs(residual) > 1e-5 * cumul_in:
                errors.append(
                    f"abs(conservation_residual_cgs)={abs(residual)} exceeds "
                    f"1e-5*cumul_reservoir_in_cgs={1e-5 * cumul_in}: "
                    f"{row.get('__line', '')}"
                )
        elif abs(residual) > 1e-30:
            errors.append(
                f"conservation_residual_cgs={residual} expected near zero "
                "when cumul_reservoir_in_cgs is zero"
            )
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8]}, str(case.log))
    return pass_result(
        name,
        f"checked {len(thermal_rows)} THERMAL rows",
        log=str(case.log),
    )


def find_hdf5_candidates(cdir):
    candidates = []
    for pattern in ("DD*/data*.cpu*", "RD*/RedshiftOutput*.cpu*", "*.h5", "*.hdf5"):
        candidates.extend(cdir.glob(pattern))
    return sorted({p for p in candidates if p.is_file()})


def flatten_numeric(data):
    try:
        iterator = data.reshape(-1)
    except AttributeError:
        iterator = [data]
    values = []
    for value in iterator:
        values.append(float(value))
    return values


def inspect_hdf5_with_h5py(files):
    import h5py

    info = {
        "method": f"h5py via {sys.executable}",
        "files_inspected": [],
        "present": {},
        "stale": {},
        "values": {label: [] for label in HDF5_EXPECTED_LABELS},
        "errors": [],
    }
    for path in files:
        try:
            with h5py.File(path, "r") as h5f:
                info["files_inspected"].append(str(path))

                def collect(name, obj):
                    base = name.rsplit("/", 1)[-1]
                    if base in HDF5_EXPECTED_LABELS:
                        info["present"].setdefault(base, []).append(f"/{name}")
                        try:
                            info["values"][base].extend(flatten_numeric(obj[()]))
                        except Exception as err:
                            info["errors"].append(
                                f"{path}: failed reading /{name}: {err}")
                    if base in HDF5_STALE_LABELS:
                        info["stale"].setdefault(base, []).append(f"/{name}")

                h5f.visititems(collect)
        except Exception as err:
            info["errors"].append(f"{path}: {err}")
    return info


def inspect_hdf5_with_helper_python(files):
    helper_pythons = [
        DEFAULT_CONDA_BIN / "python",
        Path("/home/bkoh/miniconda3/envs/yt-env/bin/python"),
    ]
    code = r"""
import json
import sys
import h5py

expected = set(sys.argv[1].split(","))
stale = set(sys.argv[2].split(","))
paths = sys.argv[3:]
info = {
    "method": "h5py helper",
    "files_inspected": [],
    "present": {},
    "stale": {},
    "values": {label: [] for label in expected},
    "errors": [],
}
for path in paths:
    try:
        with h5py.File(path, "r") as h5f:
            info["files_inspected"].append(path)
            def collect(name, obj):
                base = name.rsplit("/", 1)[-1]
                if base in expected:
                    info["present"].setdefault(base, []).append("/" + name)
                    try:
                        data = obj[()]
                        try:
                            iterator = data.reshape(-1)
                        except AttributeError:
                            iterator = [data]
                        for value in iterator:
                            info["values"][base].append(float(value))
                    except Exception as err:
                        info["errors"].append(f"{path}: failed reading /{name}: {err}")
                if base in stale:
                    info["stale"].setdefault(base, []).append("/" + name)
            h5f.visititems(collect)
    except Exception as err:
        info["errors"].append(f"{path}: {err}")
print(json.dumps(info))
"""
    for helper in helper_pythons:
        if not helper.exists():
            continue
        proc = run_cmd(
            [helper, "-c", code, ",".join(HDF5_EXPECTED_LABELS),
             ",".join(HDF5_STALE_LABELS), *files],
            cwd=ROOT,
        )
        if proc.returncode == 0:
            info = json.loads(proc.stdout)
            info["method"] = f"h5py via {helper}"
            return info
    return None


def h5ls_paths(path):
    proc = run_cmd([H5LS, "-r", path], cwd=ROOT)
    if proc.returncode != 0:
        return [], proc.stdout
    paths = []
    for line in proc.stdout.splitlines():
        stripped = line.strip()
        if not stripped.startswith("/"):
            continue
        paths.append(stripped.split()[0])
    return paths, ""


def h5dump_values(path, dataset_path):
    proc = run_cmd([H5DUMP, "-d", dataset_path, path], cwd=ROOT)
    if proc.returncode != 0:
        return [], proc.stdout
    match = re.search(r"DATA\s*\{(?P<data>.*?)\n\s*\}", proc.stdout, re.S)
    data = match.group("data") if match else proc.stdout
    values = []
    for token in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?", data):
        values.append(float(token))
    return values, ""


def inspect_hdf5_with_cli(files):
    info = {
        "method": f"{H5LS} + {H5DUMP}",
        "files_inspected": [],
        "present": {},
        "stale": {},
        "values": {label: [] for label in HDF5_EXPECTED_LABELS},
        "errors": [],
    }
    if not executable_available(H5LS) or not executable_available(H5DUMP):
        info["errors"].append(f"HDF5 CLI unavailable: H5LS={H5LS} H5DUMP={H5DUMP}")
        return info
    for path in files:
        paths, err = h5ls_paths(path)
        if err:
            info["errors"].append(f"{path}: h5ls failed: {err.strip()}")
            continue
        info["files_inspected"].append(str(path))
        for dataset_path in paths:
            base = dataset_path.rsplit("/", 1)[-1]
            if base in HDF5_EXPECTED_LABELS:
                info["present"].setdefault(base, []).append(dataset_path)
                values, val_err = h5dump_values(path, dataset_path)
                if val_err:
                    info["errors"].append(
                        f"{path}: h5dump failed for {dataset_path}: {val_err.strip()}")
                info["values"][base].extend(values)
            if base in HDF5_STALE_LABELS:
                info["stale"].setdefault(base, []).append(dataset_path)
    return info


def inspect_hdf5_outputs(cdir):
    files = find_hdf5_candidates(cdir)
    if not files:
        return {
            "method": "none",
            "files_inspected": [],
            "present": {},
            "stale": {},
            "values": {label: [] for label in HDF5_EXPECTED_LABELS},
            "errors": [f"no HDF5 candidate output files found under {cdir}"],
        }
    try:
        return inspect_hdf5_with_h5py(files)
    except ImportError:
        helper_info = inspect_hdf5_with_helper_python(files)
        if helper_info is not None:
            return helper_info
        return inspect_hdf5_with_cli(files)


def validate_t6(case):
    name = "T6 HDF5 label regression guard"
    if case.rc != 0:
        return fail_result(name, f"T1 run unavailable rc={case.rc}", log=str(case.log))
    info = inspect_hdf5_outputs(case.cdir)
    errors = []
    if not info["files_inspected"]:
        errors.append("; ".join(info["errors"]) or "no HDF5 files inspected")
    for label in HDF5_EXPECTED_LABELS:
        if not info["present"].get(label):
            errors.append(f"{label} absent from inspected HDF5 outputs")
        values = info["values"].get(label, [])
        if not values:
            errors.append(f"{label} values were not readable")
        elif not all(math.isfinite(float(value)) for value in values):
            errors.append(f"{label} contains non-finite values")
    for label in HDF5_STALE_LABELS:
        if info["stale"].get(label):
            errors.append(f"stale label {label} present at {info['stale'][label][:3]}")
    if errors:
        return fail_result(
            name,
            errors[0],
            {"errors": errors[:8], "hdf5": info},
            str(case.log),
        )
    return pass_result(
        name,
        f"labels present={','.join(HDF5_EXPECTED_LABELS)}; "
        f"stale labels absent; inspected={len(info['files_inspected'])}; "
        f"method={info['method']}",
        {"hdf5": info},
        str(case.log),
    )


def find_gas_present_burst(rows):
    previous_out = {}
    candidates = []
    for row in rows:
        bh_id = row.get("bh_id")
        cumul_out = row.get("cumul_reservoir_out_cgs")
        previous = previous_out.get(bh_id, 0.0)
        if is_finite_number(cumul_out):
            previous_out[bh_id] = float(cumul_out)
        ok_conservation, conservation_detail = conservation_check(row)
        fields_ok = (
            row.get("feedback_mode") == "THERMAL" and
            row.get("kernel_gas_zero") == 0 and
            row.get("burst_diag") == 1 and
            row.get("burst_occurred") == 1 and
            is_finite_number(row.get("E_deposited")) and
            float(row["E_deposited"]) > 0.0 and
            is_finite_number(row.get("reservoir_final")) and
            abs(float(row["reservoir_final"])) <= TINY_ABS_TOL and
            is_finite_number(cumul_out) and
            float(cumul_out) > previous + max(TINY_ABS_TOL, abs(previous) * 1e-15) and
            is_finite_number(row.get("E_requested")) and
            float(row["E_requested"]) > F1_BURST_THRESHOLD_CGS and
            ok_conservation
        )
        candidate = compact_row(row, (
            "step", "level", "bh_id", "feedback_mode", "f_Edd",
            "E_requested", "reservoir_before", "reservoir_after_accum",
            "reservoir_final", "burst_diag", "burst_occurred",
            "E_deposited", "kernel_gas_zero", "cumul_reservoir_in_cgs",
            "cumul_reservoir_out_cgs", "conservation_residual_cgs",
        ))
        candidate["previous_cumul_reservoir_out_cgs"] = previous
        candidate["conservation"] = conservation_detail
        candidate["passes_f1"] = fields_ok
        candidates.append(candidate)
        if fields_ok:
            return row, candidate, candidates
    return None, None, candidates


def validate_d0b_f1(case):
    name = "D0b-lite F1 gas-present burst"
    if case.timed_out:
        return fail_result(name, f"run timed out after {CASE_TIMEOUT}s", log=str(case.log))
    if case.rc != 0:
        return fail_result(name, f"run exited rc={case.rc}", log=str(case.log))

    rows = all_thermal_bhfdbk_rows(case.log)
    if not rows:
        return fail_result(
            name,
            "no THERMAL [BHFDBK] rows found; cannot exercise burst path",
            log=str(case.log),
        )
    burst_row, burst_summary, candidates = find_gas_present_burst(rows)
    if burst_row is None:
        return fail_result(
            name,
            "no gas-present deposited burst found; kernel_gas_zero=1 or "
            "no-deposition rows do not satisfy F1",
            {"parsed_bhfdbk_candidates": candidates[:8],
             "threshold_cgs": F1_BURST_THRESHOLD_CGS},
            str(case.log),
        )
    return pass_result(
        name,
        "gas-present burst deposited energy with reservoir reset and "
        "conservation residual within tolerance",
        {"burst": burst_summary, "threshold_cgs": F1_BURST_THRESHOLD_CGS},
        str(case.log),
    )


def validate_d0b_f3(case, f1_result):
    name = "D0b-lite F3 conservation audit through burst"
    if case.rc != 0:
        return fail_result(name, f"F1 burst run unavailable rc={case.rc}", log=str(case.log))
    rows = all_thermal_bhfdbk_rows(case.log)
    errors = []
    max_abs_residual = 0.0
    max_relative_residual = 0.0
    max_tolerance = 0.0
    recompute_max_abs_delta = 0.0
    burst_rows = []
    for row in rows:
        ok, detail = conservation_check(row)
        if not ok:
            errors.append(f"{detail}: {row.get('__line', '')}")
            continue
        cumul_in = float(row["cumul_reservoir_in_cgs"])
        cumul_out = float(row.get("cumul_reservoir_out_cgs", 0.0))
        reservoir_final = float(row.get("reservoir_final", 0.0))
        residual = float(row["conservation_residual_cgs"])
        tolerance = conservation_tolerance(cumul_in)
        max_abs_residual = max(max_abs_residual, abs(residual))
        if cumul_in > 0.0:
            max_relative_residual = max(max_relative_residual, abs(residual) / cumul_in)
        max_tolerance = max(max_tolerance, tolerance)
        recomputed = cumul_in - cumul_out - reservoir_final
        recompute_delta = abs(recomputed - residual)
        recompute_max_abs_delta = max(recompute_max_abs_delta, recompute_delta)
        if recompute_delta > max(tolerance, TINY_ABS_TOL):
            errors.append(
                "recomputed cumul_in-cumul_out-reservoir_final differs from "
                f"logged residual by {recompute_delta:.8e}: {row.get('__line', '')}"
            )
        if row.get("burst_occurred") == 1 and row.get("kernel_gas_zero") == 0:
            burst_rows.append(row)
    if not rows:
        errors.append("no THERMAL [BHFDBK] rows found")
    if f1_result.status != "PASS":
        errors.append("F1 did not identify a successful gas-present burst")
    for row in burst_rows:
        if not is_finite_number(row.get("cumul_reservoir_out_cgs")):
            errors.append(f"burst row missing finite cumul_reservoir_out_cgs: {row}")
        elif float(row["cumul_reservoir_out_cgs"]) <= 0.0:
            errors.append(f"burst row cumul_reservoir_out_cgs <= 0: {row}")
        if not is_finite_number(row.get("reservoir_final")):
            errors.append(f"burst row missing finite reservoir_final: {row}")
        elif abs(float(row["reservoir_final"])) > TINY_ABS_TOL:
            errors.append(f"burst row reservoir_final not near zero: {row}")
    if not burst_rows:
        errors.append("no gas-present burst row available for burst-specific F3 checks")
    metrics = {
        "thermal_rows_checked": len(rows),
        "burst_rows_checked": len(burst_rows),
        "max_abs_conservation_residual_cgs": max_abs_residual,
        "max_relative_conservation_residual": max_relative_residual,
        "max_conservation_tolerance_cgs": max_tolerance,
        "recomputed_balance_checked": True,
        "max_recomputed_residual_delta_cgs": recompute_max_abs_delta,
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        f"checked {len(rows)} THERMAL rows; max_abs_residual={max_abs_residual:.8e} "
        f"max_tolerance={max_tolerance:.8e}",
        metrics,
        str(case.log),
    )


def bhseed_rows(log, min_step=0):
    return rows_for(log, "[BHSEED]", min_step=min_step)


def bhseed_seed_rows(log, min_step=0):
    seed_rows = []
    for line in lines(log, "[BHSEED_SEED]"):
        row = parse_kv(line)
        step = row.get("step")
        if isinstance(step, int) and step < min_step:
            continue
        seed_rows.append(row)
    return seed_rows


def diagnostic_lines(log, tag, warn_tag=None):
    return [
        line for line in lines(log, tag)
        if not warn_tag or warn_tag not in line
    ]


def bhrepos_rows(log, min_step=0):
    return rows_for(log, "[BHREPOS]", min_step=min_step)


def bhrepos_warn_rows(log, min_step=0):
    rows = []
    for line in lines(log, "[BHREPOS_WARN]"):
        row = parse_kv(line)
        step = row.get("step")
        if isinstance(step, int) and step < min_step:
            continue
        rows.append(row)
    return rows


def seed_created_total(log, min_step=0):
    seed_rows = bhseed_rows(log, min_step=min_step)
    if not seed_rows:
        return None
    total = 0
    for row in seed_rows:
        if is_finite_number(row.get("created")):
            total += int(row["created"])
    return total


def validate_case_completed(name, case):
    if case.timed_out:
        return f"{name} run timed out after {CASE_TIMEOUT}s"
    if case.rc != 0:
        return f"{name} run exited rc={case.rc}"
    return ""


def required_field_errors(row, fields):
    missing = [field for field in fields if field not in row]
    non_finite = [
        field for field in fields
        if field in row and not is_finite_number(row[field])
    ]
    errors = []
    if missing:
        errors.append(
            f"missing required [BHACCR] fields {missing}: {row.get('__line', '')}"
        )
    if non_finite:
        errors.append(
            f"non-finite required [BHACCR] fields {non_finite}: {row.get('__line', '')}"
        )
    return errors


def select_d0c_acc_row(rows, name, predicate):
    candidates = [row for row in rows if predicate(row)]
    if not candidates:
        return None, f"{name}: no [BHACCR] row matched requested regime"
    return candidates[0], ""


def assert_close(errors, label, actual, expected, rel_tol, abs_tol=D0C_ACC_STRICT_ATOL):
    if not math.isclose(float(actual), float(expected), rel_tol=rel_tol, abs_tol=abs_tol):
        rel = relative_diff(float(actual), float(expected))
        errors.append(
            f"{label}: actual={actual:.12e} expected={expected:.12e} "
            f"rel_diff={rel:.8e} rel_tol={rel_tol:.8e}"
        )


def d0c_bondi_mdot_msunyr(row):
    bh_mass_g = float(row["bh_mass"]) * CGS_SOLAR_MASS
    rho = float(row["rho_cold_avg"])
    cs = float(row["cs_cold_avg"])
    v_rel = float(row["v_rel_cold"])
    denom = (cs * cs + v_rel * v_rel) ** 1.5
    if denom <= 0.0:
        return math.nan
    mdot_g_s = (
        float(row["f_AM"]) * 4.0 * math.pi * CGS_GRAV_CONST * CGS_GRAV_CONST *
        bh_mass_g * bh_mass_g * rho / denom
    )
    return mdot_g_s * CGS_YR_S / CGS_SOLAR_MASS


def validate_d0c_a1(case):
    name = "D0c-acc A1 Bondi formula/unit"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    rows = bhaccr_rows(case.log)
    row, row_error = select_d0c_acc_row(
        rows,
        name,
        lambda item: is_finite_number(item.get("Mdot_cold_raw")) and
        float(item["Mdot_cold_raw"]) > 0.0,
    )
    if row is None:
        return fail_result(name, row_error, {"bhaccr_rows": len(rows)}, str(case.log))

    errors = required_field_errors(row, D0C_ACC_REQUIRED_FIELDS)
    seed_total = seed_created_total(case.log, min_step=0)
    if seed_total is not None and seed_total != 1:
        errors.append(f"seed_created_total={seed_total} expected 1")
    if not errors:
        if float(row["f_hot"]) != 0.0:
            errors.append(f"f_hot={row['f_hot']} expected 0")
        assert_close(errors, "f_cold", float(row["f_cold"]), 1.0, 0.0, 0.0)
        assert_close(errors, "f_AM", float(row["f_AM"]), 1.0, 0.0, 0.0)
        assert_close(errors, "Mdot_hot_raw", float(row["Mdot_hot_raw"]), 0.0, 0.0, 0.0)
        if float(row["Mdot_cold_raw"]) <= 0.0:
            errors.append(f"Mdot_cold_raw={row['Mdot_cold_raw']} expected > 0")
        assert_close(
            errors, "Mdot_total_raw", float(row["Mdot_total_raw"]),
            float(row["Mdot_cold_raw"]), D0C_ACC_STRICT_RTOL,
        )
        if int(row["cap_active"]) != 0:
            errors.append(f"cap_active={row['cap_active']} expected 0")
        assert_close(errors, "frac_gas", float(row["frac_gas"]), 1.0, 1e-12, 1e-12)
        if int(row["removal_gas_limited"]) != 0:
            errors.append(
                f"removal_gas_limited={row['removal_gas_limited']} expected 0"
            )
        if float(row["rho_cold_avg"]) <= 0.0:
            errors.append(f"rho_cold_avg={row['rho_cold_avg']} expected > 0")
        if float(row["cs_cold_avg"]) <= 0.0:
            errors.append(f"cs_cold_avg={row['cs_cold_avg']} expected > 0")

    analytic = math.nan
    rel_error = math.nan
    if not errors:
        analytic = d0c_bondi_mdot_msunyr(row)
        logged = float(row["Mdot_cold_raw"])
        rel_error = relative_diff(analytic, logged)
        if not math.isfinite(analytic):
            errors.append("analytic Mdot_cold_raw is non-finite")
        elif rel_error > D0C_ACC_BONDI_RTOL:
            errors.append(
                f"analytic Mdot_cold_raw={analytic:.12e} logged={logged:.12e} "
                f"relative_error={rel_error:.8e} exceeds {D0C_ACC_BONDI_RTOL:.8e}"
            )

    metrics = {
        "seed_created_total": seed_total,
        "row": compact_row(row, D0C_ACC_REQUIRED_FIELDS),
        "analytic_mdot_cold_raw": analytic,
        "logged_mdot_cold_raw": row.get("Mdot_cold_raw"),
        "relative_error": rel_error,
        "formula": (
            "f_AM*4*pi*G^2*M_bh^2*rho_cold_avg/"
            "(cs_cold_avg^2+v_rel_cold^2)^(3/2)"
        ),
        "unit_convention": {
            "rho_cold_avg": "g/cm^3",
            "cs_cold_avg": "cm/s",
            "v_rel_cold": "cm/s",
            "bh_mass": "Msun",
            "Mdot_cold_raw": "Msun/yr",
        },
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        f"analytic/logged relative_error={rel_error:.8e}",
        metrics,
        str(case.log),
    )


def validate_d0c_a2(case):
    name = "D0c-acc A2 Eddington cap activation"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    rows = bhaccr_rows(case.log)
    row, row_error = select_d0c_acc_row(
        rows,
        name,
        lambda item: item.get("cap_active") == 1,
    )
    if row is None:
        return fail_result(name, row_error, {"bhaccr_rows": len(rows)}, str(case.log))

    errors = required_field_errors(row, D0C_ACC_REQUIRED_FIELDS)
    seed_total = seed_created_total(case.log, min_step=0)
    if seed_total is not None and seed_total != 1:
        errors.append(f"seed_created_total={seed_total} expected 1")
    if not errors:
        if int(row["cap_active"]) != 1:
            errors.append(f"cap_active={row['cap_active']} expected 1")
        if float(row["lambda_edd"]) <= 1.0:
            errors.append(f"lambda_edd={row['lambda_edd']} expected > 1")
        assert_close(errors, "edd_factor", float(row["edd_factor"]), 1.0, 0.0, 0.0)
        if not (float(row["Mdot_total_raw"]) > float(row["Mdot_Edd"]) > 0.0):
            errors.append(
                f"expected Mdot_total_raw > Mdot_Edd > 0, got "
                f"{row['Mdot_total_raw']} and {row['Mdot_Edd']}"
            )
        assert_close(
            errors, "Mdot_actual", float(row["Mdot_actual"]),
            float(row["Mdot_Edd"]), D0C_ACC_STRICT_RTOL,
        )
        if not (0.0 < float(row["frac_cap"]) < 1.0):
            errors.append(f"frac_cap={row['frac_cap']} expected 0 < frac_cap < 1")
        assert_close(errors, "frac_gas", float(row["frac_gas"]), 1.0, 1e-12, 1e-12)
        if int(row["removal_gas_limited"]) != 0:
            errors.append(
                f"removal_gas_limited={row['removal_gas_limited']} expected 0"
            )
        assert_close(
            errors, "dm_removed", float(row["dm_removed"]),
            float(row["dm_requested"]), D0C_ACC_STRICT_RTOL,
        )
        if float(row["dm_removed"]) <= 0.0:
            errors.append(f"dm_removed={row['dm_removed']} expected > 0")

    metrics = {
        "seed_created_total": seed_total,
        "row": compact_row(row, D0C_ACC_REQUIRED_FIELDS),
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        f"cap_active=1 lambda_edd={float(row['lambda_edd']):.8e} "
        f"frac_cap={float(row['frac_cap']):.8e}",
        metrics,
        str(case.log),
    )


def validate_d0c_a5(case):
    name = "D0c-acc A5 gas-limited mismatch"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    rows = bhaccr_rows(case.log)
    row, row_error = select_d0c_acc_row(
        rows,
        name,
        lambda item: item.get("removal_gas_limited") == 1,
    )
    if row is None:
        return fail_result(name, row_error, {"bhaccr_rows": len(rows)}, str(case.log))

    errors = required_field_errors(row, D0C_ACC_REQUIRED_FIELDS)
    seed_total = seed_created_total(case.log, min_step=0)
    if seed_total is not None and seed_total != 1:
        errors.append(f"seed_created_total={seed_total} expected 1")
    if not errors:
        if float(row["Mdot_actual"]) <= 0.0:
            errors.append(f"Mdot_actual={row['Mdot_actual']} expected > 0")
        if float(row["dm_requested"]) <= 0.0:
            errors.append(f"dm_requested={row['dm_requested']} expected > 0")
        if float(row["dm_removed"]) <= 0.0:
            errors.append(f"dm_removed={row['dm_removed']} expected > 0")
        if not (float(row["dm_removed"]) < float(row["dm_requested"])):
            errors.append(
                f"dm_removed={row['dm_removed']} expected < "
                f"dm_requested={row['dm_requested']}"
            )
        if int(row["removal_gas_limited"]) != 1:
            errors.append(
                f"removal_gas_limited={row['removal_gas_limited']} expected 1"
            )
        if not (0.0 < float(row["frac_gas"]) < 1.0):
            errors.append(f"frac_gas={row['frac_gas']} expected 0 < frac_gas < 1")
        if int(row["cap_active"]) != 0:
            errors.append(f"cap_active={row['cap_active']} expected 0")
        assert_close(errors, "frac_cap", float(row["frac_cap"]), 1.0, 1e-12, 1e-12)

    metrics = {
        "seed_created_total": seed_total,
        "row": compact_row(row, D0C_ACC_REQUIRED_FIELDS),
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        f"frac_gas={float(row['frac_gas']):.8e} "
        f"dm_removed/dm_requested={float(row['dm_removed']) / float(row['dm_requested']):.8e}",
        metrics,
        str(case.log),
    )


def validate_d0c_f6(case):
    name = "D0c-acc F6 mismatch metric"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    rows = [
        row for row in bhaccr_rows(case.log)
        if is_finite_number(row.get("dm_requested")) and
        is_finite_number(row.get("dm_removed")) and
        float(row["dm_requested"]) > 0.0
    ]
    if not rows:
        return fail_result(
            name,
            "no A5 [BHACCR] rows with dm_requested > 0",
            {"bhaccr_rows": len(bhaccr_rows(case.log))},
            str(case.log),
        )
    denominator = sum(float(row["dm_requested"]) for row in rows)
    numerator = sum(float(row["dm_requested"]) - float(row["dm_removed"]) for row in rows)
    metric = numerator / denominator if denominator > 0.0 else math.nan
    row_failures = [
        compact_row(row, ("step", "bh_id", "dm_requested", "dm_removed", "frac_gas"))
        for row in rows
        if not (float(row["dm_removed"]) < float(row["dm_requested"]))
    ]
    errors = []
    if denominator <= 0.0:
        errors.append(f"denominator={denominator} expected > 0")
    if not (0.0 < metric < 1.0):
        errors.append(f"mismatch metric M={metric} expected 0 < M < 1")
    if metric <= D0C_F6_MIN_MISMATCH:
        errors.append(f"mismatch metric M={metric} expected > {D0C_F6_MIN_MISMATCH}")
    if row_failures:
        errors.append(f"rows without dm_removed < dm_requested: {row_failures[:3]}")
    metrics = {
        "formula": "sum(dm_requested - dm_removed) / sum(dm_requested)",
        "numerator": numerator,
        "denominator": denominator,
        "mismatch_metric": metric,
        "rows_used": len(rows),
        "representative_row": compact_row(
            rows[0],
            ("step", "bh_id", "dm_requested", "dm_removed", "frac_gas",
             "removal_gas_limited", "cap_active"),
        ),
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        f"M={metric:.10f} numerator={numerator:.8e} denominator={denominator:.8e}",
        metrics,
        str(case.log),
    )


def validate_d2_1_static_source():
    name = "D2-1 append-only realized mdot source audit"
    enzo_dir = ROOT / "src/enzo"
    typedefs = (enzo_dir / "typedefs.h").read_text()
    errors = []

    layout_match = re.search(
        r"#ifdef WINDS\s*"
        r"(const int PARTICLE_ATTRIBUTE_BHSEED_CHANNEL.*?"
        r"PARTICLE_ATTRIBUTE_BHFDBK_REQUIRED\s*=\s*.*?;)"
        r"\s*#else\s*"
        r"(const int PARTICLE_ATTRIBUTE_BHSEED_CHANNEL.*?"
        r"PARTICLE_ATTRIBUTE_BHFDBK_REQUIRED\s*=\s*.*?;)"
        r"\s*#endif",
        typedefs,
        re.S,
    )
    if layout_match is None:
        return fail_result(name, "could not split WINDS/non-WINDS typedefs blocks")
    winds_block = layout_match.group(1)
    nonwinds_block = layout_match.group(2)

    required_snippets = (
        (
            "WINDS realized index",
            "PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED = 23",
            winds_block,
        ),
        (
            "non-WINDS realized index",
            "PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED = 22",
            nonwinds_block,
        ),
        (
            "WINDS LAST_MDOT_ACTUAL preserved",
            "PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_ACTUAL = 20",
            winds_block,
        ),
        (
            "non-WINDS LAST_MDOT_ACTUAL preserved",
            "PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_ACTUAL = 17",
            nonwinds_block,
        ),
        (
            "non-WINDS feedback reservoir preserved",
            "PARTICLE_ATTRIBUTE_BHFDBK_ENERGY_RESERVOIR = 18",
            nonwinds_block,
        ),
        (
            "non-WINDS cumulative reservoir out preserved",
            "PARTICLE_ATTRIBUTE_BH_CUMUL_RESERVOIR_OUT  = 21",
            nonwinds_block,
        ),
        (
            "required count uses realized attribute",
            "PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED + 1",
            typedefs,
        ),
    )
    for label, snippet, text in required_snippets:
        if snippet not in text:
            errors.append(f"{label} missing snippet {snippet!r}")

    for fname in D2_1_IO_LABEL_FILES:
        text = (enzo_dir / fname).read_text()
        if "bhaccr_last_mdot_realized" not in text:
            errors.append(f"{fname} missing bhaccr_last_mdot_realized label")
        if fname != "Grid_WriteGridX.C" and "particle_attribute_22" in text:
            errors.append(f"{fname} still contains non-WINDS append placeholder particle_attribute_22")

    acc_text = (enzo_dir / "Grid_BHAccretionHandler.C").read_text()
    for snippet in (
            "PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED",
            "Mdot_realized=",
            "mdot_realized_code=",
            "realized_over_requested=",
            "realized_attr_written="):
        if snippet not in acc_text:
            errors.append(f"Grid_BHAccretionHandler.C missing {snippet}")

    feedback_text = (enzo_dir / "Grid_BHFeedbackHandler.C").read_text()
    for snippet in (
            "PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED",
            "feedback_mdot_basis=%s",
            "requested_actual",
            "mdot_realized_code=",
            "E_requested_realized_basis="):
        if snippet not in feedback_text:
            errors.append(f"Grid_BHFeedbackHandler.C missing {snippet}")

    if errors:
        return fail_result(name, errors[0], {"errors": errors[:12]})
    return pass_result(
        name,
        "non-WINDS index 22, WINDS index 23, labels and diagnostics present",
    )


def validate_d2_1_realized_rate(case):
    name = "D2-1 realized mdot accretion diagnostics"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    rows = bhaccr_rows(case.log)
    row, row_error = select_d0c_acc_row(
        rows,
        name,
        lambda item: item.get("removal_gas_limited") == 1 and
        is_finite_number(item.get("dm_requested")) and
        is_finite_number(item.get("dm_removed")) and
        float(item["dm_requested"]) > 0.0,
    )
    if row is None:
        return fail_result(name, row_error, {"bhaccr_rows": len(rows)}, str(case.log))

    errors = required_field_errors(row, D2_1_BHACCR_REQUIRED_FIELDS)
    ratio = math.nan
    if not errors:
        requested = float(row["Mdot_actual"])
        realized = float(row["Mdot_realized"])
        dm_requested = float(row["dm_requested"])
        dm_removed = float(row["dm_removed"])
        ratio = dm_removed / dm_requested
        if requested <= 0.0:
            errors.append(f"Mdot_actual={requested} expected > 0")
        if realized <= 0.0:
            errors.append(f"Mdot_realized={realized} expected > 0")
        if not (realized < requested):
            errors.append(
                f"Mdot_realized={realized:.8e} expected < Mdot_actual={requested:.8e}"
            )
        if int(row["realized_attr_written"]) != 1:
            errors.append(
                f"realized_attr_written={row['realized_attr_written']} expected 1"
            )
        assert_close(
            errors, "realized/requested from Mdot fields",
            realized / requested, ratio, 5e-6, 1e-12,
        )
        assert_close(
            errors, "realized_over_requested",
            float(row["realized_over_requested"]), ratio, 5e-6, 1e-12,
        )
        assert_close(
            errors, "mdot code ratio",
            float(row["mdot_realized_code"]) / float(row["mdot_actual_code"]),
            ratio, 5e-6, 1e-12,
        )

    metrics = {
        "row": compact_row(row, D2_1_BHACCR_REQUIRED_FIELDS),
        "dm_removed_over_requested": ratio,
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        f"Mdot_realized/Mdot_actual={float(row['Mdot_realized']) / float(row['Mdot_actual']):.8e}",
        metrics,
        str(case.log),
    )


def validate_d2_1_feedback_diagnostics(case):
    name = "D2-1 feedback diagnostic read keeps requested basis"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    acc_rows = {
        (row.get("step"), row.get("bh_id")): row
        for row in bhaccr_rows(case.log)
        if row.get("removal_gas_limited") == 1
    }
    rows = [
        row for row in thermal_bhfdbk_rows(case.log)
        if (row.get("step"), row.get("bh_id")) in acc_rows
    ]
    if not rows:
        return fail_result(
            name,
            "no THERMAL [BHFDBK] row matched gas-limited [BHACCR] row",
            {"bhaccr_gas_limited_rows": len(acc_rows)},
            str(case.log),
        )
    row = rows[0]
    acc_row = acc_rows[(row.get("step"), row.get("bh_id"))]

    missing = [field for field in D2_1_BHFDBK_REQUIRED_FIELDS if field not in row]
    non_finite = [
        field for field in D2_1_BHFDBK_REQUIRED_FIELDS
        if field != "feedback_mdot_basis" and
        field in row and not is_finite_number(row[field])
    ]
    errors = []
    if missing:
        errors.append(
            f"missing required [BHFDBK] fields {missing}: {row.get('__line', '')}"
        )
    if non_finite:
        errors.append(
            f"non-finite required [BHFDBK] fields {non_finite}: {row.get('__line', '')}"
        )
    if not errors:
        if row["feedback_mdot_basis"] != "requested_actual":
            errors.append(
                f"feedback_mdot_basis={row['feedback_mdot_basis']!r} expected requested_actual"
            )
        if int(row["realized_attr_present"]) != 1:
            errors.append(f"realized_attr_present={row['realized_attr_present']} expected 1")
        if int(row["realized_attr_anomalous"]) != 0:
            errors.append(f"realized_attr_anomalous={row['realized_attr_anomalous']} expected 0")
        assert_close(
            errors, "BHFDBK mdot_actual_code vs BHACCR",
            float(row["mdot_actual_code"]), float(acc_row["mdot_actual_code"]),
            1e-6, 1e-30,
        )
        assert_close(
            errors, "BHFDBK mdot_realized_code vs BHACCR",
            float(row["mdot_realized_code"]), float(acc_row["mdot_realized_code"]),
            1e-6, 1e-30,
        )
        if not (float(row["mdot_realized_code"]) < float(row["mdot_actual_code"])):
            errors.append(
                "expected gas-limited realized mdot below requested/capped mdot: "
                f"{row.get('__line', '')}"
            )
        if not (float(row["E_requested_realized_basis"]) < float(row["E_requested"])):
            errors.append(
                "feedback energy appears switched or comparison missing: "
                f"E_requested_realized_basis={row['E_requested_realized_basis']} "
                f"E_requested={row['E_requested']}"
            )
        if not (float(row["L_feedback_realized_basis"]) < float(row["L_feedback"])):
            errors.append(
                "L_feedback comparison does not preserve requested/capped basis"
            )

    metrics = {
        "bhaccr_row": compact_row(acc_row, D2_1_BHACCR_REQUIRED_FIELDS),
        "bhfdbk_row": compact_row(row, D2_1_BHFDBK_REQUIRED_FIELDS),
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        "feedback_mdot_basis=requested_actual with lower realized comparison",
        metrics,
        str(case.log),
    )


def dump_sort_key(path):
    match = re.search(r"(\d+)$", path.name)
    return int(match.group(1)) if match else -1


def parse_output_metadata(param_path):
    meta = {"path": str(param_path), "time": None, "cycle": None}
    text = param_path.read_text(errors="ignore")
    patterns = {
        "time": r"^InitialTime\s*=\s*([^\s#]+)",
        "cycle": r"^InitialCycleNumber\s*=\s*([^\s#]+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text, re.MULTILINE)
        if match:
            value = match.group(1)
            meta[key] = float(value) if key == "time" else int(float(value))
    return meta


def output_parameter_files(cdir, prefix="DD"):
    candidates = []
    for child in (sorted(cdir.iterdir()) if cdir.exists() else []):
        if not child.is_dir() or not child.name.startswith(prefix):
            continue
        for param in (child / child.name, *sorted(child.glob("data[0-9]*"))):
            if not param.is_file():
                continue
            if param.name == child.name or re.fullmatch(r"data\d+", param.name):
                candidates.append(param)
    return candidates


def final_output_metadata(cdir, prefix="DD"):
    files = output_parameter_files(cdir, prefix=prefix)
    if not files:
        return None
    metas = [parse_output_metadata(path) for path in files]
    return max(
        metas,
        key=lambda item: (
            item["time"] if item["time"] is not None else -math.inf,
            item["cycle"] if item["cycle"] is not None else -1,
            dump_sort_key(Path(item["path"])),
        ),
    )


def hdf5_files_for_output(param_path):
    param_path = Path(param_path)
    parent = param_path.parent
    base = param_path.name
    files = sorted(parent.glob(f"{base}.cpu*"))
    if files:
        return files
    return find_hdf5_candidates(parent)


def dataset_values(dataset):
    data = dataset[()]
    try:
        iterator = data.reshape(-1)
    except AttributeError:
        iterator = [data]
    values = []
    for value in iterator:
        if hasattr(value, "item"):
            value = value.item()
        values.append(value)
    return values


def maybe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def collect_group_particle_records(group, file_path, group_name):
    required_any = ("particle_index", "particle_mass", "particle_type")
    if not any(name in group for name in required_any):
        return [], []
    arrays = {}
    labels = {
        "particle_index": "ParticleID",
        "particle_type": "ParticleType",
        "particle_mass": "ParticleMass",
        "particle_position_x": "ParticlePosition_x",
        "particle_position_y": "ParticlePosition_y",
        "particle_position_z": "ParticlePosition_z",
        "particle_velocity_x": "ParticleVelocity_x",
        "particle_velocity_y": "ParticleVelocity_y",
        "particle_velocity_z": "ParticleVelocity_z",
    }
    for field in T7_RESTART_FIELDS:
        if field.startswith("Particle"):
            continue
        labels[field] = field
    errors = []
    for dataset_name, label in labels.items():
        if dataset_name not in group:
            continue
        try:
            arrays[label] = dataset_values(group[dataset_name])
        except Exception as err:
            errors.append(f"{file_path}:{group_name}/{dataset_name}: {err}")
    if "ParticleMass" not in arrays:
        return [], errors
    n_particles = len(arrays["ParticleMass"])
    records = []
    for index in range(n_particles):
        rec = {
            "__file": str(file_path),
            "__group": group_name or "/",
        }
        for label, values in arrays.items():
            if index < len(values):
                rec[label] = maybe_float(values[index])
        if all(key in rec for key in (
                "ParticlePosition_x", "ParticlePosition_y", "ParticlePosition_z")):
            rec["ParticlePosition"] = (
                float(rec.pop("ParticlePosition_x")),
                float(rec.pop("ParticlePosition_y")),
                float(rec.pop("ParticlePosition_z")),
            )
        if all(key in rec for key in (
                "ParticleVelocity_x", "ParticleVelocity_y", "ParticleVelocity_z")):
            rec["ParticleVelocity"] = (
                float(rec.pop("ParticleVelocity_x")),
                float(rec.pop("ParticleVelocity_y")),
                float(rec.pop("ParticleVelocity_z")),
            )
        ptype = rec.get("ParticleType")
        if ptype is not None and int(abs(float(ptype))) not in BH_PARTICLE_TYPES:
            continue
        if "ParticleID" in rec:
            rec["ParticleID"] = int(float(rec["ParticleID"]))
        records.append(rec)
    return records, errors


def collect_bh_particle_records_with_helper(files):
    helper_pythons = [
        DEFAULT_CONDA_BIN / "python",
        Path("/home/bkoh/miniconda3/envs/yt-env/bin/python"),
    ]
    code = r"""
import json
import sys
import h5py

fields = sys.argv[1].split(",")
files = sys.argv[2:]
bh_types = {6, 8}
labels = {
    "particle_index": "ParticleID",
    "particle_type": "ParticleType",
    "particle_mass": "ParticleMass",
    "particle_position_x": "ParticlePosition_x",
    "particle_position_y": "ParticlePosition_y",
    "particle_position_z": "ParticlePosition_z",
    "particle_velocity_x": "ParticleVelocity_x",
    "particle_velocity_y": "ParticleVelocity_y",
    "particle_velocity_z": "ParticleVelocity_z",
}
for field in fields:
    if not field.startswith("Particle"):
        labels[field] = field

def dataset_values(dataset):
    data = dataset[()]
    try:
        iterator = data.reshape(-1)
    except AttributeError:
        iterator = [data]
    values = []
    for value in iterator:
        if hasattr(value, "item"):
            value = value.item()
        try:
            values.append(float(value))
        except Exception:
            values.append(value)
    return values

def collect_group(group, file_path, group_name):
    if not any(name in group for name in ("particle_index", "particle_mass", "particle_type")):
        return [], []
    arrays = {}
    errors = []
    for dataset_name, label in labels.items():
        if dataset_name not in group:
            continue
        try:
            arrays[label] = dataset_values(group[dataset_name])
        except Exception as err:
            errors.append(f"{file_path}:{group_name}/{dataset_name}: {err}")
    if "ParticleMass" not in arrays:
        return [], errors
    n_particles = len(arrays["ParticleMass"])
    records = []
    for index in range(n_particles):
        rec = {"__file": file_path, "__group": group_name or "/"}
        for label, values in arrays.items():
            if index < len(values):
                rec[label] = values[index]
        if all(key in rec for key in ("ParticlePosition_x", "ParticlePosition_y", "ParticlePosition_z")):
            rec["ParticlePosition"] = [
                float(rec.pop("ParticlePosition_x")),
                float(rec.pop("ParticlePosition_y")),
                float(rec.pop("ParticlePosition_z")),
            ]
        if all(key in rec for key in ("ParticleVelocity_x", "ParticleVelocity_y", "ParticleVelocity_z")):
            rec["ParticleVelocity"] = [
                float(rec.pop("ParticleVelocity_x")),
                float(rec.pop("ParticleVelocity_y")),
                float(rec.pop("ParticleVelocity_z")),
            ]
        ptype = rec.get("ParticleType")
        if ptype is not None and int(abs(float(ptype))) not in bh_types:
            continue
        if "ParticleID" in rec:
            rec["ParticleID"] = int(float(rec["ParticleID"]))
        records.append(rec)
    return records, errors

records = []
errors = []
for path in files:
    try:
        with h5py.File(path, "r") as h5f:
            groups = [("", h5f)]
            def collect_groups(name, obj):
                if hasattr(obj, "keys"):
                    groups.append((name, obj))
            h5f.visititems(collect_groups)
            for group_name, group in groups:
                group_records, group_errors = collect_group(group, path, group_name)
                records.extend(group_records)
                errors.extend(group_errors)
    except Exception as err:
        errors.append(f"{path}: {err}")
print(json.dumps({"records": records, "errors": errors}))
"""
    for helper in helper_pythons:
        if not helper.exists():
            continue
        proc = run_cmd(
            [helper, "-c", code, ",".join(T7_RESTART_FIELDS), *files],
            cwd=ROOT,
        )
        if proc.returncode == 0:
            try:
                info = json.loads(proc.stdout)
            except json.JSONDecodeError as err:
                return [], [f"{helper}: failed decoding T7 HDF5 helper output: {err}"]
            return info.get("records", []), info.get("errors", [])
    return [], ["h5py unavailable for T7 particle comparison in python3 and helper envs"]


def collect_bh_particle_records(param_path):
    files = hdf5_files_for_output(param_path)
    if not files:
        return [], [f"no HDF5 files found for output parameter file {param_path}"]
    try:
        import h5py
    except ImportError as err:
        helper_records, helper_errors = collect_bh_particle_records_with_helper(files)
        if helper_errors == ["h5py unavailable for T7 particle comparison in python3 and helper envs"]:
            return [], [f"h5py unavailable for T7 particle comparison: {err}"]
        return helper_records, helper_errors
    records = []
    errors = []
    for path in files:
        try:
            with h5py.File(path, "r") as h5f:
                groups = [("", h5f)]

                def collect_groups(name, obj):
                    if hasattr(obj, "keys"):
                        groups.append((name, obj))

                h5f.visititems(collect_groups)
                for group_name, group in groups:
                    group_records, group_errors = collect_group_particle_records(
                        group, path, group_name)
                    records.extend(group_records)
                    errors.extend(group_errors)
        except Exception as err:
            errors.append(f"{path}: {err}")
    return records, errors


def d0c_seed_failure(name, expected, actual, row=None, line=None):
    raw = line if line is not None else (row.get("__line", "") if row else "")
    parts = [f"{name}: expected {expected}; actual {actual}"]
    if raw:
        parts.append(f"raw={raw}")
    if row is not None:
        parts.append(f"parsed={row}")
    return "; ".join(parts)


def d0c_seed_expect_equal(errors, name, row, field, expected):
    actual = row.get(field, "<missing>")
    if actual != expected:
        errors.append(d0c_seed_failure(
            name, f"{field} == {expected!r}", actual, row=row))


def d0c_seed_expect_at_least(errors, name, row, field, expected):
    actual = row.get(field, "<missing>")
    if not is_finite_number(actual) or float(actual) < float(expected):
        errors.append(d0c_seed_failure(
            name, f"{field} >= {expected!r}", actual, row=row))


def d0c_seed_expect_close(
        errors, name, row, field, expected, rel_tol=0.0, abs_tol=0.0):
    actual = row.get(field, "<missing>")
    if not is_finite_number(actual) or not math.isclose(
            float(actual), float(expected), rel_tol=rel_tol, abs_tol=abs_tol):
        errors.append(d0c_seed_failure(
            name,
            f"{field} ~= {expected!r} rel_tol={rel_tol} abs_tol={abs_tol}",
            actual,
            row=row,
        ))


def d0c_seed_particle_summary(record):
    return {
        key: record.get(key)
        for key in (
            "ParticleID", "ParticleType", "ParticleMass", "ParticlePosition",
            "__file", "__group",
        )
        if key in record
    }


def d0c_seed_is_bh_particle(record):
    try:
        return int(abs(float(record.get("ParticleType")))) in BH_PARTICLE_TYPES
    except (TypeError, ValueError):
        return False


def d0c_seed_final_bh_particles(case):
    meta = final_output_metadata(case.cdir)
    if meta is None:
        return [], [f"missing final DD output metadata under {case.cdir}"], None
    records, errors = collect_bh_particle_records(meta["path"])
    bh_records = [record for record in records if d0c_seed_is_bh_particle(record)]
    return bh_records, errors, meta


def d0c_seed_single_particle(name, case):
    records, h5_errors, meta = d0c_seed_final_bh_particles(case)
    errors = []
    if h5_errors:
        errors.append(d0c_seed_failure(
            name, "readable HDF5 BH particle records", h5_errors))
    if len(records) != 1:
        errors.append(d0c_seed_failure(
            name,
            "exactly one HDF5 BH particle with ParticleType 6 or 8",
            {
                "count": len(records),
                "records": [d0c_seed_particle_summary(record) for record in records[:4]],
            },
        ))
    record = records[0] if len(records) == 1 else None
    return record, records, meta, errors


def d0c_seed_particle_position(record):
    position = record.get("ParticlePosition")
    if not isinstance(position, (tuple, list)) or len(position) != 3:
        return None
    return tuple(float(value) for value in position)


def d0c_seed_validate_particle_type_mass(
        errors, name, record, expected_mass=None, mass_abs_tol=None):
    ptype = record.get("ParticleType", "<missing>")
    try:
        actual_type = int(abs(float(ptype)))
    except (TypeError, ValueError):
        actual_type = ptype
    if actual_type != D0C_SEED_EXPECTED_PARTICLE_TYPE:
        errors.append(d0c_seed_failure(
            name,
            f"ParticleType == {D0C_SEED_EXPECTED_PARTICLE_TYPE}",
            ptype,
        ))
    mass = record.get("ParticleMass")
    if expected_mass is None:
        expected_mass = D0C_SEED_EXPECTED_MASS
        close = (
            is_finite_number(mass) and
            math.isclose(float(mass), expected_mass, rel_tol=D0C_SEED_MASS_RTOL)
        )
        expected = (
            f"ParticleMass ~= {expected_mass:.8e} "
            f"rel_tol={D0C_SEED_MASS_RTOL}"
        )
    else:
        close = (
            is_finite_number(mass) and
            math.isclose(
                float(mass), float(expected_mass),
                rel_tol=0.0, abs_tol=mass_abs_tol or 0.0,
            )
        )
        expected = (
            f"ParticleMass matches S1 {float(expected_mass):.17e} "
            f"abs_tol={mass_abs_tol}"
        )
    if not close:
        errors.append(d0c_seed_failure(name, expected, mass))


def validate_d0c_seed_s0(case):
    name = "D0c-seed S0 default-off no-seed"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    errors = []
    seed_lines = lines(case.log, "[BHSEED]")
    seed_seed_lines = lines(case.log, "[BHSEED_SEED]")
    bhaccr_lines = diagnostic_lines(case.log, "[BHACCR]", warn_tag="[BHACCR_WARN]")
    bhfdbk_lines = diagnostic_lines(case.log, "[BHFDBK]", warn_tag="[BHFDBK_WARN]")
    if seed_lines:
        errors.append(d0c_seed_failure(
            name, "zero [BHSEED] lines", len(seed_lines), line=seed_lines[0]))
    if seed_seed_lines:
        errors.append(d0c_seed_failure(
            name, "zero [BHSEED_SEED] lines", len(seed_seed_lines),
            line=seed_seed_lines[0],
        ))
    if bhaccr_lines:
        errors.append(d0c_seed_failure(
            name, "zero non-WARN [BHACCR] lines", len(bhaccr_lines),
            line=bhaccr_lines[0],
        ))
    if bhfdbk_lines:
        errors.append(d0c_seed_failure(
            name, "zero non-WARN [BHFDBK] lines", len(bhfdbk_lines),
            line=bhfdbk_lines[0],
        ))
    bh_records, h5_errors, meta = d0c_seed_final_bh_particles(case)
    if h5_errors:
        errors.append(d0c_seed_failure(
            name, "readable final HDF5 particle inventory", h5_errors))
    if bh_records:
        errors.append(d0c_seed_failure(
            name,
            "zero HDF5 BH particles with ParticleType 6 or 8",
            {
                "count": len(bh_records),
                "records": [d0c_seed_particle_summary(record)
                            for record in bh_records[:4]],
            },
        ))
    metrics = {
        "run_dir": str(case.cdir),
        "seed_lines": len(seed_lines),
        "seed_seed_lines": len(seed_seed_lines),
        "bhaccr_lines": len(bhaccr_lines),
        "bhfdbk_lines": len(bhfdbk_lines),
        "hdf5_bh_particles": len(bh_records),
        "final_output": meta,
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        "default-off produced no seed/accretion/feedback diagnostics and no HDF5 BHs",
        metrics,
        str(case.log),
    )


def validate_d0c_seed_s1(case):
    name = "D0c-seed S1 explicit opt-in seed creation"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    errors = []
    seed_total = seed_created_total(case.log, min_step=0)
    seed_rows = bhseed_seed_rows(case.log, min_step=0)
    if seed_total != 1:
        errors.append(d0c_seed_failure(
            name, "seed_created_total(log, min_step=0) == 1", seed_total))
    if len(seed_rows) != 1:
        raw = lines(case.log, "[BHSEED_SEED]")
        errors.append(d0c_seed_failure(
            name, "exactly one [BHSEED_SEED] line", len(seed_rows),
            line=raw[0] if raw else "",
        ))
    seed_row = seed_rows[0] if len(seed_rows) == 1 else {}
    if seed_row:
        for field, expected in zip(("x", "y", "z"), D0C_SEED_EXPECTED_POSITION):
            d0c_seed_expect_close(
                errors, name, seed_row, field, expected,
                abs_tol=D0C_SEED_POSITION_ATOL,
            )
        d0c_seed_expect_equal(errors, name, seed_row, "channel", 0)
        d0c_seed_expect_equal(errors, name, seed_row, "accept_rank", 1)

    summary_rows = bhseed_rows(case.log, min_step=1)
    if not summary_rows:
        errors.append(d0c_seed_failure(
            name, "at least one [BHSEED] row with min_step=1",
            len(summary_rows),
        ))
        first_summary = {}
    else:
        first_summary = summary_rows[0]
        for field, expected in (
                ("ncand_global", 2),
                ("created", 1),
                ("nseeds_created", 1),
                ("total_mbh", 1),
                ("ncandidates_gathered", 2),
                ("ncandidates_dedup_rejected", 1),
        ):
            d0c_seed_expect_equal(errors, name, first_summary, field, expected)

    bhaccr_lines = diagnostic_lines(case.log, "[BHACCR]", warn_tag="[BHACCR_WARN]")
    bhfdbk_lines = diagnostic_lines(case.log, "[BHFDBK]", warn_tag="[BHFDBK_WARN]")
    if bhaccr_lines:
        errors.append(d0c_seed_failure(
            name, "zero non-WARN [BHACCR] rows", len(bhaccr_lines),
            line=bhaccr_lines[0],
        ))
    if bhfdbk_lines:
        errors.append(d0c_seed_failure(
            name, "zero non-WARN [BHFDBK] rows", len(bhfdbk_lines),
            line=bhfdbk_lines[0],
        ))

    particle, bh_records, meta, particle_errors = d0c_seed_single_particle(name, case)
    errors.extend(particle_errors)
    if particle is not None:
        d0c_seed_validate_particle_type_mass(errors, name, particle)

    metrics = {
        "run_dir": str(case.cdir),
        "seed_created_total": seed_total,
        "seed_row": seed_row,
        "first_bhseed_row": first_summary,
        "bhaccr_lines": len(bhaccr_lines),
        "bhfdbk_lines": len(bhfdbk_lines),
        "hdf5_bh_particles": len(bh_records),
        "bh_particles": [d0c_seed_particle_summary(record) for record in bh_records],
        "final_output": meta,
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        f"seed_created_total=1; seed=({seed_row.get('x')},"
        f"{seed_row.get('y')},{seed_row.get('z')}); "
        f"ParticleMass={particle.get('ParticleMass'):.8e}",
        metrics,
        str(case.log),
    )


def validate_d0c_seed_s2(case):
    name = "D0c-seed S2 temperature gate no-seed"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    errors = []
    # This is the one controlled negative gate used in D0c-seed1.
    summary_rows = bhseed_rows(case.log, min_step=0)
    if not summary_rows:
        errors.append(d0c_seed_failure(
            name,
            "at least one [BHSEED] row with min_step=0 confirming handler ran",
            len(summary_rows),
        ))
    for row in summary_rows:
        d0c_seed_expect_at_least(errors, name, row, "ngates_temp", 2)
        for field in ("created", "nseeds_created", "total_mbh"):
            d0c_seed_expect_equal(errors, name, row, field, 0)

    seed_rows = bhseed_seed_rows(case.log, min_step=0)
    if seed_rows:
        errors.append(d0c_seed_failure(
            name, "zero [BHSEED_SEED] lines", len(seed_rows), row=seed_rows[0]))
    bh_records, h5_errors, meta = d0c_seed_final_bh_particles(case)
    if h5_errors:
        errors.append(d0c_seed_failure(
            name, "readable final HDF5 particle inventory", h5_errors))
    if bh_records:
        errors.append(d0c_seed_failure(
            name,
            "zero HDF5 BH particles with ParticleType 6 or 8",
            {
                "count": len(bh_records),
                "records": [d0c_seed_particle_summary(record)
                            for record in bh_records[:4]],
            },
        ))
    representative = summary_rows[0] if summary_rows else {}
    metrics = {
        "run_dir": str(case.cdir),
        "bhseed_rows": len(summary_rows),
        "representative_bhseed_row": representative,
        "seed_seed_lines": len(seed_rows),
        "hdf5_bh_particles": len(bh_records),
        "final_output": meta,
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        f"temperature gate rejected candidates; ngates_temp={representative.get('ngates_temp')}",
        metrics,
        str(case.log),
    )


def validate_d0c_seed_s3(case_a, case_b):
    name = "D0c-seed S3 deterministic np1 repeat"
    errors = []
    for label, case in (("repeat_a", case_a), ("repeat_b", case_b)):
        run_error = validate_case_completed(name, case)
        if run_error:
            errors.append(f"{name} {label}: {run_error}")
    total_a = seed_created_total(case_a.log, min_step=0) if case_a.rc == 0 else None
    total_b = seed_created_total(case_b.log, min_step=0) if case_b.rc == 0 else None
    if total_a != 1:
        errors.append(d0c_seed_failure(
            name, "repeat_a seed_created_total(log, min_step=0) == 1", total_a))
    if total_b != 1:
        errors.append(d0c_seed_failure(
            name, "repeat_b seed_created_total(log, min_step=0) == 1", total_b))

    seed_a = bhseed_seed_rows(case_a.log, min_step=0)
    seed_b = bhseed_seed_rows(case_b.log, min_step=0)
    if len(seed_a) != 1:
        errors.append(d0c_seed_failure(
            name, "repeat_a exactly one [BHSEED_SEED] line", len(seed_a),
            line=lines(case_a.log, "[BHSEED_SEED]")[0] if lines(case_a.log, "[BHSEED_SEED]") else "",
        ))
    if len(seed_b) != 1:
        errors.append(d0c_seed_failure(
            name, "repeat_b exactly one [BHSEED_SEED] line", len(seed_b),
            line=lines(case_b.log, "[BHSEED_SEED]")[0] if lines(case_b.log, "[BHSEED_SEED]") else "",
        ))
    compared_fields = ("x", "y", "z", "channel", "patch_density_peak", "accept_rank")
    if len(seed_a) == 1 and len(seed_b) == 1:
        for field in compared_fields:
            if seed_a[0].get(field) != seed_b[0].get(field):
                errors.append(d0c_seed_failure(
                    name,
                    f"repeat [BHSEED_SEED] field {field} exactly matches",
                    {"repeat_a": seed_a[0].get(field), "repeat_b": seed_b[0].get(field)},
                    row=seed_a[0],
                ))

    particle_a, records_a, meta_a, particle_errors_a = d0c_seed_single_particle(name, case_a)
    particle_b, records_b, meta_b, particle_errors_b = d0c_seed_single_particle(name, case_b)
    errors.extend(f"repeat_a {err}" for err in particle_errors_a)
    errors.extend(f"repeat_b {err}" for err in particle_errors_b)
    hdf5_comparison = {}
    if particle_a is not None and particle_b is not None:
        mass_a = particle_a.get("ParticleMass")
        mass_b = particle_b.get("ParticleMass")
        if not (
                is_finite_number(mass_a) and is_finite_number(mass_b) and
                math.isclose(
                    float(mass_a), float(mass_b), rel_tol=0.0,
                    abs_tol=D0C_SEED_DETERMINISTIC_MASS_ATOL,
                )):
            errors.append(d0c_seed_failure(
                name,
                "repeat HDF5 ParticleMass matches within deterministic tolerance",
                {"repeat_a": mass_a, "repeat_b": mass_b},
            ))
        pos_a = d0c_seed_particle_position(particle_a)
        pos_b = d0c_seed_particle_position(particle_b)
        if pos_a is None or pos_b is None:
            errors.append(d0c_seed_failure(
                name,
                "repeat HDF5 ParticlePosition readable on both particles",
                {"repeat_a": pos_a, "repeat_b": pos_b},
            ))
        else:
            for dim, (left, right) in enumerate(zip(pos_a, pos_b)):
                if not math.isclose(
                        left, right, rel_tol=0.0,
                        abs_tol=D0C_SEED_DETERMINISTIC_POSITION_ATOL):
                    errors.append(d0c_seed_failure(
                        name,
                        f"repeat HDF5 ParticlePosition[{dim}] matches within "
                        f"{D0C_SEED_DETERMINISTIC_POSITION_ATOL}",
                        {"repeat_a": left, "repeat_b": right},
                    ))
        hdf5_comparison = {
            "mass_a": mass_a,
            "mass_b": mass_b,
            "position_a": pos_a,
            "position_b": pos_b,
            "mass_abs_tol": D0C_SEED_DETERMINISTIC_MASS_ATOL,
            "position_abs_tol": D0C_SEED_DETERMINISTIC_POSITION_ATOL,
        }

    metrics = {
        "run_dirs": [str(case_a.cdir), str(case_b.cdir)],
        "seed_created_total_a": total_a,
        "seed_created_total_b": total_b,
        "seed_rows": {
            "repeat_a": seed_a[0] if seed_a else {},
            "repeat_b": seed_b[0] if seed_b else {},
        },
        "compared_fields": compared_fields,
        "hdf5_bh_particles": [len(records_a), len(records_b)],
        "hdf5_comparison": hdf5_comparison,
        "final_outputs": [meta_a, meta_b],
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case_b.log))
    return pass_result(
        name,
        "np=1 repeat matched seed diagnostics and HDF5 mass/position",
        metrics,
        str(case_b.log),
    )


def validate_d0c_seed_s4(case, reference_particle, mpirun_path):
    name = "D0c-seed S4 optional MPI seed consistency"
    if case.timed_out or case.rc != 0:
        reason = (
            f"np=4 mpirun launch did not complete cleanly rc={case.rc} "
            f"timed_out={case.timed_out}"
        )
        if case.timed_out or mpi_environment_failure(case.log):
            return skip_result(
                name,
                reason,
                {
                    "run_dir": str(case.cdir),
                    "mpirun": mpirun_path,
                    "log": str(case.log),
                },
                str(case.log),
            )
        return fail_result(name, reason, {"mpirun": mpirun_path}, str(case.log))

    errors = []
    seed_total = seed_created_total(case.log, min_step=0)
    if seed_total != 1:
        errors.append(d0c_seed_failure(
            name, "seed_created_total(log, min_step=0) == 1", seed_total))
    summary_rows = bhseed_rows(case.log, min_step=1)
    if not summary_rows:
        errors.append(d0c_seed_failure(
            name, "at least one [BHSEED] row with min_step=1", len(summary_rows)))
        first_summary = {}
    else:
        first_summary = summary_rows[0]
        d0c_seed_expect_equal(errors, name, first_summary, "ncand_global", 2)

    particle, records, meta, particle_errors = d0c_seed_single_particle(name, case)
    errors.extend(particle_errors)
    ref_mass = reference_particle.get("ParticleMass")
    ref_pos = d0c_seed_particle_position(reference_particle)
    position_comparison = {}
    if particle is not None:
        d0c_seed_validate_particle_type_mass(
            errors,
            name,
            particle,
            expected_mass=ref_mass,
            mass_abs_tol=D0C_SEED_DETERMINISTIC_MASS_ATOL,
        )
        pos = d0c_seed_particle_position(particle)
        if pos is None or ref_pos is None:
            errors.append(d0c_seed_failure(
                name,
                "np=4 and S1 HDF5 ParticlePosition are readable",
                {"np4": pos, "s1": ref_pos},
            ))
        else:
            for dim, (actual, expected) in enumerate(zip(pos, ref_pos)):
                if not math.isclose(
                        actual, expected, rel_tol=0.0,
                        abs_tol=D0C_SEED_MPI_POSITION_ATOL):
                    errors.append(d0c_seed_failure(
                        name,
                        f"np=4 HDF5 ParticlePosition[{dim}] matches S1 "
                        f"within {D0C_SEED_MPI_POSITION_ATOL}",
                        {"np4": actual, "s1": expected},
                    ))
            position_comparison = {"np4": pos, "s1": ref_pos}

    metrics = {
        "run_dir": str(case.cdir),
        "mpirun": mpirun_path,
        "seed_created_total": seed_total,
        "first_bhseed_row": first_summary,
        "hdf5_bh_particles": len(records),
        "bh_particles": [d0c_seed_particle_summary(record) for record in records],
        "reference_particle": d0c_seed_particle_summary(reference_particle),
        "position_comparison": position_comparison,
        "final_output": meta,
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        "np=4 seed count/type/mass/position matched S1 within D0c-seed tolerances",
        metrics,
        str(case.log),
    )


def run_d0c_seed_s4(s1_result):
    name = "D0c-seed S4 optional MPI seed consistency"
    mpirun_path = str(DEFAULT_MPIRUN)
    if not executable_available(mpirun_path):
        return skip_result(name, f"mpirun unavailable: {mpirun_path}", {"mpirun": mpirun_path})
    reference_particles = s1_result.metrics.get("bh_particles", [])
    if s1_result.status != "PASS" or len(reference_particles) != 1:
        return skip_result(
            name,
            "S1 reference particle unavailable; optional S4 not attempted",
            {"s1_status": s1_result.status, "s1_particle_count": len(reference_particles)},
        )
    s4_dir = prepare_case(
        "D0c_seed_S4_optional_mpi_np4_seed_consistency",
        d0c_seed_opt_in_overrides(),
    )
    s4 = run_enzo(s4_dir, np_ranks=4, mpirun=mpirun_path)
    return validate_d0c_seed_s4(s4, reference_particles[0], mpirun_path)


def d0c_repos_failure(name, expected, actual, row=None, line=None):
    raw = line if line is not None else (row.get("__line", "") if row else "")
    parts = [f"{name}: expected {expected}; actual {actual}"]
    if raw:
        parts.append(f"raw={raw}")
    if row is not None:
        parts.append(f"parsed={row}")
    return "; ".join(parts)


def d0c_repos_expect_equal(errors, name, row, field, expected):
    actual = row.get(field, "<missing>")
    if actual != expected:
        errors.append(d0c_repos_failure(
            name, f"{field} == {expected!r}", actual, row=row))


def d0c_repos_expect_int_in(errors, name, row, field, allowed):
    actual = row.get(field, "<missing>")
    if not isinstance(actual, int) or actual not in allowed:
        errors.append(d0c_repos_failure(
            name, f"{field} in {sorted(allowed)}", actual, row=row))


def d0c_repos_expect_close(
        errors, name, row, field, expected, rel_tol=0.0,
        abs_tol=D0C_REPOS_ZERO_ATOL):
    actual = row.get(field, "<missing>")
    if not is_finite_number(actual) or not math.isclose(
            float(actual), float(expected), rel_tol=rel_tol, abs_tol=abs_tol):
        errors.append(d0c_repos_failure(
            name,
            f"{field} ~= {expected!r} rel_tol={rel_tol} abs_tol={abs_tol}",
            actual,
            row=row,
        ))


def d0c_repos_expect_finite(errors, name, row, field, minimum=None):
    actual = row.get(field, "<missing>")
    if not is_finite_number(actual):
        errors.append(d0c_repos_failure(name, f"{field} finite", actual, row=row))
        return
    if minimum is not None and float(actual) < float(minimum):
        errors.append(d0c_repos_failure(
            name, f"{field} >= {minimum}", actual, row=row))


def d0c_repos_validate_boolean_counter_safety(errors, name, row):
    for field in (
            "diag_peak_in_ghost", "active_target_exists",
            "reposition_occurred", "reposition_clamped",
            "newly_seeded_skip", "diagnostic_only",
            "search_kernel_truncated", "active_reposition_rejected"):
        if field in row:
            d0c_repos_expect_int_in(errors, name, row, field, {0, 1})
    if "method" in row:
        d0c_repos_expect_int_in(errors, name, row, "method", {0, 1, 2})
    for field in ("search_cells", "search_active_cells"):
        if field not in row:
            continue
        actual = row.get(field)
        if (not isinstance(actual, int) or actual < 0 or
                actual > D0C_REPOS_COUNTER_MAX_SAFE):
            errors.append(d0c_repos_failure(
                name,
                f"{field} non-negative initialized counter <= "
                f"{D0C_REPOS_COUNTER_MAX_SAFE}",
                actual,
                row=row,
            ))


def d0c_repos_validate_seed_formed(errors, name, case):
    seed_total = seed_created_total(case.log, min_step=0)
    seed_rows = bhseed_seed_rows(case.log, min_step=0)
    if seed_total != 1:
        errors.append(d0c_repos_failure(
            name, "seed_created_total(log, min_step=0) == 1", seed_total))
    if len(seed_rows) != 1:
        raw = lines(case.log, "[BHSEED_SEED]")
        errors.append(d0c_repos_failure(
            name, "exactly one [BHSEED_SEED] line", len(seed_rows),
            line=raw[0] if raw else "",
        ))
        return {}
    seed_row = seed_rows[0]
    for field, expected in zip(("x", "y", "z"), D0C_SEED_EXPECTED_POSITION):
        d0c_seed_expect_close(
            errors, name, seed_row, field, expected,
            abs_tol=D0C_SEED_POSITION_ATOL,
        )
    return seed_row


def d0c_repos_regular_rows(log):
    return [row for row in bhrepos_rows(log, min_step=0)
            if row.get("newly_seeded_skip") == 0]


def d0c_repos_skip_rows(log):
    return [row for row in bhrepos_rows(log, min_step=0)
            if row.get("newly_seeded_skip") == 1]


def d0c_repos_position_from_seed(seed_row):
    if not seed_row:
        return None
    if not all(field in seed_row and is_finite_number(seed_row[field])
               for field in ("x", "y", "z")):
        return None
    return (float(seed_row["x"]), float(seed_row["y"]), float(seed_row["z"]))


def d0c_repos_position_is_finite_in_unit_grid(position):
    return (
        isinstance(position, (tuple, list)) and len(position) == 3 and
        all(math.isfinite(float(value)) and 0.0 <= float(value) <= 1.0
            for value in position)
    )


def d0c_repos_expect_position_close(
        errors, name, actual, expected, label, abs_tol=D0C_REPOS_POSITION_ATOL):
    if actual is None or expected is None:
        errors.append(d0c_repos_failure(
            name, f"{label} position readable", {"actual": actual, "expected": expected}))
        return
    for dim, (left, right) in enumerate(zip(actual, expected)):
        if not math.isclose(
                float(left), float(right), rel_tol=0.0, abs_tol=abs_tol):
            errors.append(d0c_repos_failure(
                name,
                f"{label}[{dim}] matches within {abs_tol}",
                {"actual": left, "expected": right},
            ))


def d0c_repos_final_particle(name, case, errors):
    particle, records, meta, particle_errors = d0c_seed_single_particle(name, case)
    errors.extend(particle_errors)
    position = d0c_seed_particle_position(particle) if particle is not None else None
    if position is not None and not d0c_repos_position_is_finite_in_unit_grid(position):
        errors.append(d0c_repos_failure(
            name, "final HDF5 BH position finite and inside unit grid", position))
    return particle, records, meta, position


def d0c_repos_no_accr_feedback_errors(name, case):
    errors = []
    bhaccr_lines = diagnostic_lines(case.log, "[BHACCR]", warn_tag="[BHACCR_WARN]")
    bhfdbk_lines = diagnostic_lines(case.log, "[BHFDBK]", warn_tag="[BHFDBK_WARN]")
    if bhaccr_lines:
        errors.append(d0c_repos_failure(
            name, "zero non-WARN [BHACCR] rows", len(bhaccr_lines),
            line=bhaccr_lines[0],
        ))
    if bhfdbk_lines:
        errors.append(d0c_repos_failure(
            name, "zero non-WARN [BHFDBK] rows", len(bhfdbk_lines),
            line=bhfdbk_lines[0],
        ))
    return errors, bhaccr_lines, bhfdbk_lines


def d0c_repos_r4_cell_index(position, dim=50):
    return int(round(float(position) * dim - 0.5))


def d0c_repos_r4_hdf5_edit_code():
    return r"""
import json
import sys
import h5py

h5path = sys.argv[1]
initial_pos = [float(value) for value in sys.argv[2].split(",")]
target_pos = [float(value) for value in sys.argv[3].split(",")]
background_density = float(sys.argv[4])
target_density = float(sys.argv[5])

def cell_index(position, dim=50):
    return int(round(float(position) * dim - 0.5))

with h5py.File(h5path, "r+") as h5f:
    if "Grid00000001" not in h5f:
        raise RuntimeError("Grid00000001 missing from copied R4 restart")
    group = h5f["Grid00000001"]
    required = [
        "particle_type", "particle_position_x", "particle_position_y",
        "particle_position_z", "particle_velocity_x", "particle_velocity_y",
        "particle_velocity_z", "particle_index", "Density", "GasEnergy",
        "TotalEnergy", "x-velocity", "y-velocity", "z-velocity",
    ]
    missing = [name for name in required if name not in group]
    if missing:
        raise RuntimeError("missing R4 restart datasets: " + ",".join(missing))
    particle_types = group["particle_type"][()]
    bh_index = None
    for index, ptype in enumerate(particle_types):
        if int(abs(float(ptype))) in (6, 8):
            bh_index = index
            break
    if bh_index is None:
        raise RuntimeError("no type 6/8 BH particle in copied R4 restart")

    group["particle_position_x"][bh_index] = initial_pos[0]
    group["particle_position_y"][bh_index] = initial_pos[1]
    group["particle_position_z"][bh_index] = initial_pos[2]
    group["particle_velocity_x"][bh_index] = 0.0
    group["particle_velocity_y"][bh_index] = 0.0
    group["particle_velocity_z"][bh_index] = 0.0
    group["particle_index"][bh_index] = 1

    group["Density"][...] = background_density
    group["GasEnergy"][...] = 1.0
    group["TotalEnergy"][...] = 1.0
    group["x-velocity"][...] = 0.0
    group["y-velocity"][...] = 0.0
    group["z-velocity"][...] = 0.0

    i = cell_index(target_pos[0])
    j = cell_index(target_pos[1])
    k = cell_index(target_pos[2])
    group["Density"][k, j, i] = target_density

    info = {
        "h5path": h5path,
        "group": "Grid00000001",
        "bh_index": int(bh_index),
        "datasets_edited": required,
        "initial_position": initial_pos,
        "target_position": target_pos,
        "target_cell_ijk": [int(i), int(j), int(k)],
        "background_density": background_density,
        "target_density": target_density,
    }
print(json.dumps(info))
"""


def d0c_repos_r4_edit_hdf5(h5path):
    args = [
        str(h5path),
        ",".join(str(value) for value in D0C_REPOS_R4_INITIAL_POSITION),
        ",".join(str(value) for value in D0C_REPOS_R4_TARGET_POSITION),
        str(D0C_REPOS_R4_BACKGROUND_DENSITY),
        str(D0C_REPOS_R4_TARGET_DENSITY),
    ]
    code = d0c_repos_r4_hdf5_edit_code()
    helper_pythons = [
        Path(sys.executable),
        Path("/home/bkoh/miniconda3/envs/yt-env/bin/python"),
        DEFAULT_CONDA_BIN / "python",
    ]
    failures = []
    seen = set()
    for helper in helper_pythons:
        if helper in seen:
            continue
        seen.add(helper)
        if not helper.exists():
            continue
        proc = run_cmd([helper, "-c", code, *args], cwd=ROOT)
        if proc.returncode == 0:
            return json.loads(proc.stdout)
        failures.append(f"{helper}: rc={proc.returncode} {proc.stdout[-400:]}")
    raise RuntimeError(
        "h5py unavailable for R4 copied-restart edit in python3 and helper envs; "
        + " | ".join(failures)
    )


def d0c_repos_r4_prepare_case(name, base_case, method):
    meta = final_output_metadata(base_case.cdir)
    if meta is None:
        raise RuntimeError(f"R4 base case has no final output: {base_case.cdir}")
    base_param = Path(meta["path"])
    base_output_dir = base_param.parent
    cdir = OUT_ROOT / name
    if cdir.exists():
        shutil.rmtree(cdir)
    cdir.mkdir(parents=True, exist_ok=True)
    copied_output_dir = cdir / base_output_dir.name
    shutil.copytree(base_output_dir, copied_output_dir)
    param_path = copied_output_dir / base_param.name
    text = param_path.read_text()
    for key, value in {
            "StopTime": "0.03",
            "StopCycle": "3",
            "dtDataDump": "0.01",
            "CycleSkipDataDump": "0",
            "BHSeedingMethod": "0",
            "BHSeedVerbose": "0",
            "BHAccretionMethod": "0",
            "BHAccretionRunEveryTimestep": "1",
            "BHFeedbackMethod": "0",
            "BHRepositionMethod": str(method),
            "BHRepositionSearchRadius": f"{D0C_REPOS_R4_SEARCH_RADIUS_KPC:.6g}",
            "BHRepositionMaxDisplacement": "0.5",
            "BHRepositionDiagnosePotential": "0",
            "BHRepositionVerbose": "1",
            "SelfGravity": "0",
    }.items():
        text = set_param(text, key, value)
    param_path.write_text(text)

    h5path = copied_output_dir / f"{base_param.name}.cpu0000"
    edit_info = d0c_repos_r4_edit_hdf5(h5path)
    case = run_enzo(
        cdir,
        np_ranks=1,
        log_name="run.log",
        param_file=str(param_path.relative_to(cdir)),
        restart=True,
    )
    case.metrics = {
        "base_output": str(base_output_dir),
        "copied_output": str(copied_output_dir),
        "restart_param": str(param_path),
        "hdf5_edit": edit_info,
        "method": method,
    }
    return case


def validate_d0c_repos_r0(case):
    name = "D0c-repos R0 newly seeded skip diagnostic"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    errors = []
    seed_row = d0c_repos_validate_seed_formed(errors, name, case)
    skip_rows = d0c_repos_skip_rows(case.log)
    if not skip_rows:
        errors.append(d0c_repos_failure(
            name, "at least one [BHREPOS] row with newly_seeded_skip=1",
            len(skip_rows),
        ))
        representative = {}
    else:
        representative = skip_rows[0]
        d0c_repos_expect_equal(errors, name, representative, "reposition_occurred", 0)
        d0c_repos_expect_close(errors, name, representative, "displacement_kpc", 0.0)
        d0c_repos_expect_close(errors, name, representative, "displacement_cells", 0.0)
        d0c_repos_expect_int_in(
            errors, name, representative, "search_kernel_truncated", {0, 1})
        d0c_repos_expect_int_in(
            errors, name, representative, "active_reposition_rejected", {0, 1})
        d0c_repos_expect_equal(errors, name, representative, "search_cells", 0)
        d0c_repos_expect_equal(errors, name, representative, "search_active_cells", 0)
        d0c_repos_validate_boolean_counter_safety(errors, name, representative)

    metrics = {
        "run_dir": str(case.cdir),
        "seed_row": seed_row,
        "bhrepos_skip_rows": len(skip_rows),
        "representative_bhrepos_skip_row": representative,
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        "newly seeded skip row has initialized boolean/counter diagnostics",
        metrics,
        str(case.log),
    )


def validate_d0c_repos_r1(case):
    name = "D0c-repos R1 method0 diagnostic no movement"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    errors = []
    seed_row = d0c_repos_validate_seed_formed(errors, name, case)
    seed_position = d0c_repos_position_from_seed(seed_row)
    all_rows = bhrepos_rows(case.log, min_step=0)
    regular_rows = d0c_repos_regular_rows(case.log)
    warn_lines = [
        line for line in lines(case.log, "[BHREPOS_WARN]")
        if "BHRepositionMethod=0" in line and "diagnostics only" in line
    ]
    if not all_rows:
        errors.append(d0c_repos_failure(
            name, "[BHREPOS] rows are present", len(all_rows)))
    if not warn_lines:
        errors.append(d0c_repos_failure(
            name,
            "method-0 diagnostic-only [BHREPOS_WARN] line present",
            len(warn_lines),
        ))
    if not regular_rows:
        errors.append(d0c_repos_failure(
            name, "at least one non-newly-seeded [BHREPOS] row", len(regular_rows)))
    for row in regular_rows:
        d0c_repos_expect_equal(errors, name, row, "reposition_occurred", 0)
        d0c_repos_expect_close(errors, name, row, "displacement_kpc", 0.0)
        d0c_repos_expect_close(errors, name, row, "displacement_cells", 0.0)
        d0c_repos_validate_boolean_counter_safety(errors, name, row)

    particle, records, meta, final_position = d0c_repos_final_particle(name, case, errors)
    d0c_repos_expect_position_close(
        errors, name, final_position, seed_position, "final seeded BH position")
    no_source_errors, bhaccr_lines, bhfdbk_lines = d0c_repos_no_accr_feedback_errors(
        name, case)
    errors.extend(no_source_errors)

    representative = regular_rows[0] if regular_rows else (all_rows[0] if all_rows else {})
    metrics = {
        "run_dir": str(case.cdir),
        "seed_row": seed_row,
        "seed_position": seed_position,
        "final_position": final_position,
        "bhrepos_rows": len(all_rows),
        "bhrepos_regular_rows": len(regular_rows),
        "representative_bhrepos_row": representative,
        "method0_warning": warn_lines[0] if warn_lines else "",
        "bhaccr_lines": len(bhaccr_lines),
        "bhfdbk_lines": len(bhfdbk_lines),
        "bh_particles": [d0c_seed_particle_summary(record) for record in records],
        "final_output": meta,
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        f"method=0 diagnostics emitted without movement; final_position={final_position}",
        metrics,
        str(case.log),
    )


def validate_d0c_repos_r2(case):
    name = "D0c-repos R2 method2 safety rejection"
    run_error = validate_case_completed(name, case)
    if run_error:
        return fail_result(name, run_error, log=str(case.log))

    errors = []
    seed_row = d0c_repos_validate_seed_formed(errors, name, case)
    seed_position = d0c_repos_position_from_seed(seed_row)
    all_rows = bhrepos_rows(case.log, min_step=0)
    skip_rows = d0c_repos_skip_rows(case.log)
    regular_rows = d0c_repos_regular_rows(case.log)
    warn_rows = bhrepos_warn_rows(case.log, min_step=0)
    rejection_warns = [
        row for row in warn_rows
        if row.get("active_reposition_rejected") == 1 or
        "active_reposition_rejected=1" in row.get("__line", "")
    ]
    if not all_rows:
        errors.append(d0c_repos_failure(
            name, "[BHREPOS] rows are present", len(all_rows)))
    if skip_rows:
        d0c_repos_expect_equal(errors, name, skip_rows[0], "method", 2)
    else:
        errors.append(d0c_repos_failure(
            name, "newly_seeded_skip [BHREPOS] row reflects method 2", len(skip_rows)))
    for row in rejection_warns:
        if "bh_id" not in row:
            errors.append(d0c_repos_failure(
                name, "active_reposition_rejected warning associated with bh_id",
                "<missing>",
                row=row,
            ))
        raw = row.get("__line", "")
        if "search_truncated_by_grid" not in row and "search_truncated_by_grid" not in raw:
            errors.append(d0c_repos_failure(
                name,
                "active_reposition_rejected warning reports search_truncated_by_grid",
                row,
                row=row,
            ))

    moved_rows = [row for row in regular_rows if row.get("reposition_occurred") == 1]
    for row in regular_rows:
        d0c_repos_validate_boolean_counter_safety(errors, name, row)
    particle, records, meta, final_position = d0c_repos_final_particle(name, case, errors)

    if rejection_warns or not moved_rows:
        for row in regular_rows:
            d0c_repos_expect_equal(errors, name, row, "reposition_occurred", 0)
            d0c_repos_expect_close(errors, name, row, "displacement_kpc", 0.0)
            d0c_repos_expect_close(errors, name, row, "displacement_cells", 0.0)
        d0c_repos_expect_position_close(
            errors, name, final_position, seed_position,
            "final BH position after rejected/no method2 movement",
        )
        movement_status = "rejected" if rejection_warns else "no_target_or_no_movement"
    else:
        for row in moved_rows:
            d0c_repos_expect_finite(errors, name, row, "displacement_kpc", minimum=0.0)
            d0c_repos_expect_finite(errors, name, row, "displacement_cells", minimum=0.0)
            d0c_repos_expect_equal(errors, name, row, "active_target_exists", 1)
            if (is_finite_number(row.get("displacement_cells")) and
                    float(row["displacement_cells"]) > D0C_REPOS_COUNTER_MAX_SAFE):
                errors.append(d0c_repos_failure(
                    name,
                    "method2 movement displacement_cells bounded",
                    row.get("displacement_cells"),
                    row=row,
                ))
        movement_status = "finite_movement"

    metrics = {
        "run_dir": str(case.cdir),
        "seed_row": seed_row,
        "seed_position": seed_position,
        "final_position": final_position,
        "movement_status": movement_status,
        "bhrepos_rows": len(all_rows),
        "bhrepos_regular_rows": len(regular_rows),
        "bhrepos_rejection_warnings": len(rejection_warns),
        "representative_bhrepos_row": (
            moved_rows[0] if moved_rows else
            (regular_rows[0] if regular_rows else (all_rows[0] if all_rows else {}))
        ),
        "representative_bhrepos_warn": rejection_warns[0] if rejection_warns else {},
        "bh_particles": [d0c_seed_particle_summary(record) for record in records],
        "final_output": meta,
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(case.log))
    return pass_result(
        name,
        f"method=2 active with {movement_status}; final_position={final_position}",
        metrics,
        str(case.log),
    )


def validate_d0c_repos_r3(np1_case, np4_case, mpirun_path):
    name = "D0c-repos R3 optional MPI BHREPOS consistency"
    if np4_case.timed_out or np4_case.rc != 0:
        reason = (
            f"np=4 mpirun launch did not complete cleanly rc={np4_case.rc} "
            f"timed_out={np4_case.timed_out}"
        )
        if np4_case.timed_out or mpi_environment_failure(np4_case.log):
            return skip_result(
                name,
                reason,
                {"run_dir": str(np4_case.cdir), "mpirun": mpirun_path},
                str(np4_case.log),
            )
        return fail_result(name, reason, {"mpirun": mpirun_path}, str(np4_case.log))
    if np1_case.timed_out or np1_case.rc != 0:
        return skip_result(
            name,
            "np=1 reference run unavailable; optional R3 not attempted",
            {
                "np1_rc": np1_case.rc,
                "np1_timed_out": np1_case.timed_out,
                "mpirun": mpirun_path,
            },
            str(np1_case.log),
        )

    errors = []
    seed_total_1 = seed_created_total(np1_case.log, min_step=0)
    seed_total_4 = seed_created_total(np4_case.log, min_step=0)
    if seed_total_1 != seed_total_4 or seed_total_1 != 1:
        errors.append(d0c_repos_failure(
            name,
            "np=1 and np=4 seed counts both equal 1",
            {"np1": seed_total_1, "np4": seed_total_4},
        ))

    skip_1 = d0c_repos_skip_rows(np1_case.log)
    skip_4 = d0c_repos_skip_rows(np4_case.log)
    if not skip_1:
        errors.append(d0c_repos_failure(
            name, "np=1 newly_seeded_skip BHREPOS row parseable", len(skip_1)))
        skip_row_1 = {}
    else:
        skip_row_1 = skip_1[0]
        d0c_repos_validate_boolean_counter_safety(errors, name, skip_row_1)
    if not skip_4:
        errors.append(d0c_repos_failure(
            name, "np=4 newly_seeded_skip BHREPOS row parseable", len(skip_4)))
        skip_row_4 = {}
    else:
        skip_row_4 = skip_4[0]
        d0c_repos_validate_boolean_counter_safety(errors, name, skip_row_4)
    for field in ("newly_seeded_skip", "search_kernel_truncated",
                  "active_reposition_rejected"):
        if skip_row_1 and skip_row_4 and skip_row_1.get(field) != skip_row_4.get(field):
            errors.append(d0c_repos_failure(
                name,
                f"safe skip field {field} matches np=1 and np=4",
                {"np1": skip_row_1.get(field), "np4": skip_row_4.get(field)},
                row=skip_row_4,
            ))

    particle_1, records_1, meta_1, particle_errors_1 = d0c_seed_single_particle(
        name, np1_case)
    particle_4, records_4, meta_4, particle_errors_4 = d0c_seed_single_particle(
        name, np4_case)
    errors.extend(f"np=1 {err}" for err in particle_errors_1)
    errors.extend(f"np=4 {err}" for err in particle_errors_4)
    pos_1 = d0c_seed_particle_position(particle_1) if particle_1 is not None else None
    pos_4 = d0c_seed_particle_position(particle_4) if particle_4 is not None else None
    if particle_1 is not None and particle_4 is not None:
        if int(abs(float(particle_1.get("ParticleType")))) != int(abs(float(particle_4.get("ParticleType")))):
            errors.append(d0c_repos_failure(
                name,
                "np=1 and np=4 HDF5 BH ParticleType match",
                {"np1": particle_1.get("ParticleType"), "np4": particle_4.get("ParticleType")},
            ))
        for label, particle, position in (
                ("np=1", particle_1, pos_1), ("np=4", particle_4, pos_4)):
            if not is_finite_number(particle.get("ParticleMass")):
                errors.append(d0c_repos_failure(
                    name, f"{label} ParticleMass finite", particle.get("ParticleMass")))
            if not d0c_repos_position_is_finite_in_unit_grid(position):
                errors.append(d0c_repos_failure(
                    name, f"{label} ParticlePosition finite and in unit grid", position))
        if (is_finite_number(particle_1.get("ParticleMass")) and
                is_finite_number(particle_4.get("ParticleMass")) and
                not math.isclose(
                    float(particle_1["ParticleMass"]), float(particle_4["ParticleMass"]),
                    rel_tol=0.0, abs_tol=D0C_SEED_DETERMINISTIC_MASS_ATOL)):
            errors.append(d0c_repos_failure(
                name,
                f"np=4 ParticleMass matches np=1 within "
                f"{D0C_SEED_DETERMINISTIC_MASS_ATOL}",
                {"np1": particle_1.get("ParticleMass"),
                 "np4": particle_4.get("ParticleMass")},
            ))
        d0c_repos_expect_position_close(
            errors, name, pos_4, pos_1, "np=4 vs np=1 final BH position",
            abs_tol=D0C_SEED_MPI_POSITION_ATOL,
        )

    metrics = {
        "run_dirs": {"np1": str(np1_case.cdir), "np4": str(np4_case.cdir)},
        "mpirun": mpirun_path,
        "seed_created_total": {"np1": seed_total_1, "np4": seed_total_4},
        "skip_rows": {"np1": skip_row_1, "np4": skip_row_4},
        "safe_fields_compared": (
            "newly_seeded_skip",
            "search_kernel_truncated",
            "active_reposition_rejected",
        ),
        "hdf5_bh_particles": {"np1": len(records_1), "np4": len(records_4)},
        "bh_particles": {
            "np1": [d0c_seed_particle_summary(record) for record in records_1],
            "np4": [d0c_seed_particle_summary(record) for record in records_4],
        },
        "final_outputs": {"np1": meta_1, "np4": meta_4},
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(np4_case.log))
    return pass_result(
        name,
        "np=1 and np=4 matched seed count, finite BH HDF5 fields, and safe skip diagnostics",
        metrics,
        str(np4_case.log),
    )


def run_d0c_repos_r3(np1_case):
    name = "D0c-repos R3 optional MPI BHREPOS consistency"
    mpirun_path = str(DEFAULT_MPIRUN)
    if not executable_available(mpirun_path):
        return skip_result(name, f"mpirun unavailable: {mpirun_path}", {"mpirun": mpirun_path})
    np4_dir = prepare_case(
        "D0c_repos_R3_optional_mpi_np4_bhrepos_consistency",
        d0c_repos_base_overrides(BHRepositionMethod=0, BHRepositionVerbose=1),
    )
    np4_case = run_enzo(np4_dir, np_ranks=4, mpirun=mpirun_path)
    return validate_d0c_repos_r3(np1_case, np4_case, mpirun_path)


def d0c_repos_r4_representative_row(case):
    rows = [
        row for row in d0c_repos_regular_rows(case.log)
        if row.get("bh_id") == 1
    ]
    if rows:
        return rows[0]
    rows = d0c_repos_regular_rows(case.log)
    return rows[0] if rows else {}


def d0c_repos_r4_expect_target(errors, name, row):
    d0c_repos_expect_equal(errors, name, row, "active_target_exists", 1)
    for axis, expected in zip(("x", "y", "z"), D0C_REPOS_R4_TARGET_POSITION):
        d0c_repos_expect_close(
            errors,
            name,
            row,
            f"active_peak_pos_{axis}",
            expected,
            abs_tol=D0C_REPOS_R4_POSITION_ATOL,
        )


def validate_d0c_repos_r4(method0_case, method2_case):
    name = "D0c-repos R4 actual movement ON/OFF"
    errors = []
    metrics = {
        "run_dirs": {
            "method0": str(method0_case.cdir),
            "method2": str(method2_case.cdir),
        },
        "case_setup": {
            "initial_position": D0C_REPOS_R4_INITIAL_POSITION,
            "target_position": D0C_REPOS_R4_TARGET_POSITION,
            "search_radius_kpc": D0C_REPOS_R4_SEARCH_RADIUS_KPC,
            "background_density": D0C_REPOS_R4_BACKGROUND_DENSITY,
            "target_density": D0C_REPOS_R4_TARGET_DENSITY,
            "method0": getattr(method0_case, "metrics", {}),
            "method2": getattr(method2_case, "metrics", {}),
        },
    }
    for label, case in (("method0", method0_case), ("method2", method2_case)):
        run_error = validate_case_completed(f"{name} {label}", case)
        if run_error:
            errors.append(run_error)
        no_source_errors, bhaccr_lines, bhfdbk_lines = d0c_repos_no_accr_feedback_errors(
            f"{name} {label}", case)
        errors.extend(no_source_errors)
        metrics.setdefault("diagnostic_counts", {})[label] = {
            "bhaccr_lines": len(bhaccr_lines),
            "bhfdbk_lines": len(bhfdbk_lines),
        }

    row0 = d0c_repos_r4_representative_row(method0_case)
    row2 = d0c_repos_r4_representative_row(method2_case)
    if not row0:
        errors.append(d0c_repos_failure(
            name, "method0 has at least one regular [BHREPOS] row", 0))
    if not row2:
        errors.append(d0c_repos_failure(
            name, "method2 has at least one regular [BHREPOS] row", 0))

    if row0:
        d0c_repos_r4_expect_target(errors, f"{name} method0", row0)
        d0c_repos_expect_equal(errors, f"{name} method0", row0, "reposition_occurred", 0)
        d0c_repos_expect_close(errors, f"{name} method0", row0, "displacement_cells", 0.0)
        d0c_repos_expect_close(errors, f"{name} method0", row0, "displacement_kpc", 0.0)
        d0c_repos_validate_boolean_counter_safety(errors, f"{name} method0", row0)
    if row2:
        d0c_repos_r4_expect_target(errors, f"{name} method2", row2)
        for field in ("search_kernel_truncated", "active_reposition_rejected"):
            if field in row2:
                d0c_repos_expect_equal(errors, f"{name} method2", row2, field, 0)
        d0c_repos_expect_equal(errors, f"{name} method2", row2, "reposition_occurred", 1)
        d0c_repos_expect_equal(errors, f"{name} method2", row2, "reposition_clamped", 0)
        d0c_repos_expect_close(
            errors, f"{name} method2", row2, "displacement_cells", 2.0,
            abs_tol=D0C_REPOS_R4_DISPLACEMENT_ATOL,
        )
        d0c_repos_expect_close(
            errors, f"{name} method2", row2, "displacement_kpc", 4.0,
            abs_tol=D0C_REPOS_R4_DISPLACEMENT_ATOL,
        )
        d0c_repos_validate_boolean_counter_safety(errors, f"{name} method2", row2)

    particle0, records0, meta0, final0 = d0c_repos_final_particle(
        f"{name} method0", method0_case, errors)
    particle2, records2, meta2, final2 = d0c_repos_final_particle(
        f"{name} method2", method2_case, errors)
    d0c_repos_expect_position_close(
        errors,
        f"{name} method0",
        final0,
        D0C_REPOS_R4_INITIAL_POSITION,
        "method0 final BH position",
        abs_tol=D0C_REPOS_R4_POSITION_ATOL,
    )
    d0c_repos_expect_position_close(
        errors,
        f"{name} method2",
        final2,
        D0C_REPOS_R4_TARGET_POSITION,
        "method2 final BH position",
        abs_tol=D0C_REPOS_R4_POSITION_ATOL,
    )
    if row2:
        logged_target = tuple(float(row2.get(f"active_peak_pos_{axis}"))
                              for axis in ("x", "y", "z"))
        d0c_repos_expect_position_close(
            errors,
            f"{name} method2",
            final2,
            logged_target,
            "method2 final BH position vs logged active peak",
            abs_tol=D0C_REPOS_R4_POSITION_ATOL,
        )

    metrics.update({
        "representative_bhrepos": {"method0": row0, "method2": row2},
        "initial_position": D0C_REPOS_R4_INITIAL_POSITION,
        "final_positions": {"method0": final0, "method2": final2},
        "bh_particles": {
            "method0": [d0c_seed_particle_summary(record) for record in records0],
            "method2": [d0c_seed_particle_summary(record) for record in records2],
        },
        "final_outputs": {"method0": meta0, "method2": meta2},
    })
    if errors:
        detail = errors[0]
        if row0 or row2:
            detail += (
                f"; method0_row={row0}; method2_row={row2}; "
                f"initial={D0C_REPOS_R4_INITIAL_POSITION}; "
                f"final0={final0}; final2={final2}"
            )
        return fail_result(name, detail, {"errors": errors[:10], **metrics},
                           str(method2_case.log))
    return pass_result(
        name,
        "method0 saw the target without moving; method2 moved to the logged active target",
        metrics,
        str(method2_case.log),
    )


def run_d0c_repos_r4(base_case):
    method0_case = d0c_repos_r4_prepare_case(
        "D0c_repos_R4A_method0_no_movement",
        base_case,
        method=0,
    )
    method2_case = d0c_repos_r4_prepare_case(
        "D0c_repos_R4B_method2_actual_movement",
        base_case,
        method=2,
    )
    return validate_d0c_repos_r4(method0_case, method2_case)


def unique_id_map(records):
    out = {}
    for record in records:
        pid = record.get("ParticleID")
        if pid is None or pid in out:
            return None
        out[pid] = record
    return out


def nearest_position_matches(left_records, right_records):
    remaining = list(right_records)
    matches = []
    for left in left_records:
        if "ParticlePosition" not in left:
            return None
        best = None
        best_dist = math.inf
        for right in remaining:
            if "ParticlePosition" not in right:
                return None
            dist = math.sqrt(sum(
                (float(left["ParticlePosition"][dim]) -
                 float(right["ParticlePosition"][dim])) ** 2
                for dim in range(3)
            ))
            if dist < best_dist:
                best_dist = dist
                best = right
        if best is None:
            return None
        remaining.remove(best)
        matches.append((f"nearest_pos_{len(matches)}", left, best))
    return matches


def match_bh_records(continuous_records, restart_records):
    left_ids = unique_id_map(continuous_records)
    right_ids = unique_id_map(restart_records)
    if left_ids is not None and right_ids is not None:
        if set(left_ids) != set(right_ids):
            missing_restart = sorted(set(left_ids) - set(right_ids))
            missing_cont = sorted(set(right_ids) - set(left_ids))
            return None, (
                f"BH particle ID mismatch: missing_restart={missing_restart} "
                f"missing_continuous={missing_cont}"
            )
        return [(pid, left_ids[pid], right_ids[pid]) for pid in sorted(left_ids)], ""
    if len(continuous_records) != len(restart_records):
        return None, (
            f"BH particle count mismatch continuous={len(continuous_records)} "
            f"restart={len(restart_records)}"
        )
    matches = nearest_position_matches(continuous_records, restart_records)
    if matches is None:
        return None, "could not match BHs by particle ID or nearest position"
    return matches, ""


def value_components(value):
    if isinstance(value, (tuple, list)):
        return [float(v) for v in value]
    return [float(value)]


def relative_diff(a, b):
    denom = max(abs(a), abs(b), TINY_ABS_TOL)
    return abs(a - b) / denom


def restart_tolerances_for_field(field, c_value, r_value):
    # Continuous-vs-restart evolved BH floating fields show about 1e-5
    # relative roundoff in TS3_wrap. Identity, position, velocity, seed, and
    # expected-zero reservoir-out fields retain the stricter tolerance.
    rtol = EVOLVED_BH_FLOAT_RTOL if field in T7_EVOLVED_BH_FLOAT_FIELDS else RESTART_FLOAT_RTOL
    if field == "bh_cumul_reservoir_out" and c_value == 0.0 and r_value == 0.0:
        return RESTART_FLOAT_RTOL, TINY_ABS_TOL
    return rtol, RESTART_FLOAT_ATOL


def compare_restart_field(field, key, continuous, restarted):
    if field not in continuous or field not in restarted:
        return {
            "field": field,
            "key": key,
            "error": f"missing field continuous={field in continuous} restart={field in restarted}",
            "continuous_file": continuous.get("__file"),
            "restart_file": restarted.get("__file"),
        }
    if field in T7_CATEGORICAL_FIELDS:
        if int(round(float(continuous[field]))) != int(round(float(restarted[field]))):
            return {
                "field": field,
                "key": key,
                "continuous": continuous[field],
                "restart": restarted[field],
                "error": "categorical mismatch",
                "continuous_file": continuous.get("__file"),
                "restart_file": restarted.get("__file"),
            }
        return None
    c_values = value_components(continuous[field])
    r_values = value_components(restarted[field])
    if len(c_values) != len(r_values):
        return {
            "field": field,
            "key": key,
            "continuous": continuous[field],
            "restart": restarted[field],
            "error": "component count mismatch",
            "continuous_file": continuous.get("__file"),
            "restart_file": restarted.get("__file"),
        }
    for c_value, r_value in zip(c_values, r_values):
        abs_delta = abs(c_value - r_value)
        rel_delta = relative_diff(c_value, r_value)
        rtol, atol = restart_tolerances_for_field(field, c_value, r_value)
        if abs_delta > atol and rel_delta > rtol:
            return {
                "field": field,
                "key": key,
                "continuous": continuous[field],
                "restart": restarted[field],
                "abs_diff": abs_delta,
                "rel_diff": rel_delta,
                "rtol": rtol,
                "atol": atol,
                "continuous_file": continuous.get("__file"),
                "restart_file": restarted.get("__file"),
            }
    return None


def update_field_stats(stats, field, continuous, restarted):
    if field not in continuous or field not in restarted:
        return
    if field in T7_CATEGORICAL_FIELDS:
        stat = stats.setdefault(
            field,
            {"max_abs_diff": 0.0, "max_rel_diff": 0.0, "categorical_mismatches": 0},
        )
        mismatch = int(round(float(continuous[field]))) != int(round(float(restarted[field])))
        stat["categorical_mismatches"] += 1 if mismatch else 0
        return
    c_values = value_components(continuous[field])
    r_values = value_components(restarted[field])
    stat = stats.setdefault(field, {"max_abs_diff": 0.0, "max_rel_diff": 0.0})
    for c_value, r_value in zip(c_values, r_values):
        stat["max_abs_diff"] = max(stat["max_abs_diff"], abs(c_value - r_value))
        stat["max_rel_diff"] = max(stat["max_rel_diff"], relative_diff(c_value, r_value))


def run_t7_restart_cases():
    continuous_dir = prepare_case(
        "D0b_T7_continuous_no_burst",
        overrides=d0b_no_burst_overrides(
            StopTime="0.02",
            dtDataDump="0.01",
            CycleSkipDataDump=1,
            CycleSkipRestartDump=0,
            DataDumpName="DD",
            DataDumpDir="DD",
        ),
    )
    continuous = run_enzo(continuous_dir, log_name="continuous.log")

    restart_dir = prepare_case(
        "D0b_T7_restart_no_burst",
        overrides=d0b_no_burst_overrides(
            StopTime="0.02",
            dtDataDump="0.01",
            CycleSkipDataDump=1,
            CycleSkipRestartDump=1,
            DataDumpName="DD",
            DataDumpDir="DD",
            RestartDumpName="RS",
            RestartDumpDir="RS",
        ),
    )
    first = run_enzo(restart_dir, log_name="run_to_checkpoint.log")
    restart_param_file = restart_dir / "RestartParamFile"
    restart_param = None
    if restart_param_file.exists():
        restart_param = restart_param_file.read_text(errors="ignore").strip()
    second = None
    if restart_param:
        second = run_enzo(
            restart_dir,
            log_name="restart_continue.log",
            param_file=restart_param,
            restart=True,
        )
    return continuous, first, second, restart_param


def validate_d0b_t7(continuous, first, second, restart_param):
    name = "D0b-lite T7 restart preservation"
    run_errors = []
    if continuous.timed_out or continuous.rc != 0:
        run_errors.append(f"continuous run rc={continuous.rc} timed_out={continuous.timed_out}")
    if first.timed_out or first.rc != 0:
        run_errors.append(f"checkpoint run rc={first.rc} timed_out={first.timed_out}")
    if not restart_param:
        run_errors.append(f"RestartParamFile missing after checkpoint run in {first.cdir}")
    if second is None:
        run_errors.append("restart continuation was not launched")
    elif second.timed_out or second.rc != 0:
        run_errors.append(f"restart continuation rc={second.rc} timed_out={second.timed_out}")
    if run_errors:
        return fail_result(
            name,
            run_errors[0],
            {"errors": run_errors},
            str(continuous.log),
        )

    cont_meta = final_output_metadata(continuous.cdir)
    rest_meta = final_output_metadata(second.cdir)
    if cont_meta is None or rest_meta is None:
        return fail_result(
            name,
            f"missing final DD output metadata continuous={cont_meta} restart={rest_meta}",
            {"continuous_log": str(continuous.log), "restart_log": str(second.log)},
            str(second.log),
        )
    meta_errors = []
    if cont_meta["cycle"] != rest_meta["cycle"]:
        meta_errors.append(
            f"final cycle mismatch continuous={cont_meta['cycle']} restart={rest_meta['cycle']}"
        )
    if cont_meta["time"] is None or rest_meta["time"] is None:
        meta_errors.append(
            f"final time missing continuous={cont_meta['time']} restart={rest_meta['time']}"
        )
    elif not math.isclose(
            float(cont_meta["time"]), float(rest_meta["time"]),
            rel_tol=RESTART_FLOAT_RTOL, abs_tol=RESTART_FLOAT_ATOL):
        meta_errors.append(
            f"final time mismatch continuous={cont_meta['time']} restart={rest_meta['time']}"
        )
    cont_records, cont_errors = collect_bh_particle_records(cont_meta["path"])
    rest_records, rest_errors = collect_bh_particle_records(rest_meta["path"])
    errors = meta_errors + cont_errors + rest_errors
    if not cont_records or not rest_records:
        errors.append(
            f"BH records missing continuous={len(cont_records)} restart={len(rest_records)}"
        )
    matches, match_error = match_bh_records(cont_records, rest_records)
    if match_error:
        errors.append(match_error)
    field_stats = {}
    compared_fields = []
    if matches:
        for key, continuous_record, restart_record in matches:
            for field in T7_RESTART_FIELDS:
                if field == "ParticleVelocity" and (
                        field not in continuous_record and field not in restart_record):
                    continue
                failure = compare_restart_field(field, key, continuous_record, restart_record)
                if failure:
                    errors.append(
                        f"field={field} key={key} continuous={failure.get('continuous')} "
                        f"restart={failure.get('restart')} abs_diff={failure.get('abs_diff')} "
                        f"rel_diff={failure.get('rel_diff')} files="
                        f"{failure.get('continuous_file')} {failure.get('restart_file')}"
                    )
                    if len(errors) >= 8:
                        break
                else:
                    compared_fields.append(field)
                    update_field_stats(field_stats, field, continuous_record, restart_record)
            if len(errors) >= 8:
                break
    metrics = {
        "continuous_final": cont_meta,
        "restart_final": rest_meta,
        "restart_param": restart_param,
        "bh_count_continuous": len(cont_records),
        "bh_count_restart": len(rest_records),
        "fields_compared": sorted(set(compared_fields)),
        "field_differences": field_stats,
        "logs": {
            "continuous": str(continuous.log),
            "checkpoint": str(first.log),
            "restart": str(second.log) if second else "",
        },
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(second.log))
    return pass_result(
        name,
        f"compared {len(matches)} BHs at cycle={cont_meta['cycle']} time={cont_meta['time']}",
        metrics,
        str(second.log),
    )


def mpi_environment_failure(log_path):
    text = log_path.read_text(errors="ignore") if log_path and log_path.exists() else ""
    markers = (
        "Permission denied",
        "socket",
        "ORTE_ERROR_LOG",
        "PMIX ERROR",
        "MPI_Init",
        "not enough slots",
        "There are not enough slots",
        "orted",
        "mpirun was unable",
        "orterun was unable",
    )
    return any(marker in text for marker in markers)


def diagnostic_key(row, tag):
    if tag == "BHFDBK":
        return (tag, row.get("step"), row.get("bh_id"), row.get("feedback_mode"))
    return (tag, row.get("step"), row.get("bh_id"))


def diagnostic_map(rows, tag):
    out = {}
    duplicates = []
    for row in rows:
        key = diagnostic_key(row, tag)
        if key in out:
            duplicates.append(key)
        out[key] = row
    return out, duplicates


def compare_numeric_diagnostic_field(tag, key, field, aliases, left, right):
    left_value = row_value(left, aliases)
    right_value = row_value(right, aliases)
    if left_value is None and right_value is None:
        return None, None
    if left_value is None or right_value is None:
        return {
            "tag": tag,
            "key": key,
            "field": field,
            "np1": left_value,
            "npN": right_value,
            "error": "field present on only one side",
        }, None
    if not is_finite_number(left_value) or not is_finite_number(right_value):
        return {
            "tag": tag,
            "key": key,
            "field": field,
            "np1": left_value,
            "npN": right_value,
            "error": "non-finite value",
        }, None
    left_float = float(left_value)
    right_float = float(right_value)
    abs_delta = abs(left_float - right_float)
    rel_delta = relative_diff(left_float, right_float)
    if abs_delta > T10_ABS_TOL and rel_delta > T10_SHORT_RTOL:
        return {
            "tag": tag,
            "key": key,
            "field": field,
            "np1": left_float,
            "npN": right_float,
            "abs_diff": abs_delta,
            "rel_diff": rel_delta,
            "rtol": T10_SHORT_RTOL,
            "atol": T10_ABS_TOL,
        }, None
    return None, {
        "field": field,
        "abs_diff": abs_delta,
        "rel_diff": rel_delta,
    }


def compare_categorical_diagnostic_field(tag, key, field, left, right):
    if field not in left and field not in right:
        return None
    if field not in left or field not in right or left[field] != right[field]:
        return {
            "tag": tag,
            "key": key,
            "field": field,
            "np1": left.get(field),
            "npN": right.get(field),
            "error": "categorical mismatch",
        }
    return None


def update_diag_stats(stats, tag, key, compared):
    if compared is None:
        return
    stat_key = f"{tag}.{compared['field']}"
    stat = stats.setdefault(stat_key, {"max_abs_diff": 0.0, "max_rel_diff": 0.0})
    stat["max_abs_diff"] = max(stat["max_abs_diff"], compared["abs_diff"])
    stat["max_rel_diff"] = max(stat["max_rel_diff"], compared["rel_diff"])


def gas_present_burst_keys(rows):
    keys = set()
    summaries = []
    for row in rows:
        if (
            row.get("feedback_mode") == "THERMAL" and
            row.get("kernel_gas_zero") == 0 and
            row.get("burst_diag") == 1 and
            row.get("burst_occurred") == 1 and
            is_finite_number(row.get("E_deposited")) and
            float(row["E_deposited"]) > 0.0
        ):
            key = diagnostic_key(row, "BHFDBK")
            keys.add(key)
            summaries.append(compact_row(row, (
                "step", "level", "bh_id", "feedback_mode", "E_requested",
                "E_deposited", "reservoir_before", "reservoir_after_accum",
                "reservoir_final", "burst_diag", "burst_occurred",
                "kernel_gas_zero", "cumul_reservoir_in_cgs",
                "cumul_reservoir_out_cgs", "conservation_residual_cgs",
            )))
    return keys, summaries


def compare_diagnostic_sets(tag, left_rows, right_rows, numeric_fields, categorical_fields):
    left_map, left_duplicates = diagnostic_map(left_rows, tag)
    right_map, right_duplicates = diagnostic_map(right_rows, tag)
    errors = []
    stats = {}
    compared_fields = set()
    if left_duplicates:
        errors.append(f"{tag} duplicate np=1 keys: {left_duplicates[:5]}")
    if right_duplicates:
        errors.append(f"{tag} duplicate np={T10_NP} keys: {right_duplicates[:5]}")
    if not left_map or not right_map:
        errors.append(f"{tag} diagnostics missing np1={len(left_map)} np{T10_NP}={len(right_map)}")
    if set(left_map) != set(right_map):
        errors.append(
            f"{tag} diagnostic key mismatch "
            f"only_np1={sorted(set(left_map) - set(right_map), key=str)[:5]} "
            f"only_np{T10_NP}={sorted(set(right_map) - set(left_map), key=str)[:5]}"
        )
    for key in sorted(set(left_map) & set(right_map), key=str):
        left = left_map[key]
        right = right_map[key]
        for field, aliases in numeric_fields:
            failure, compared = compare_numeric_diagnostic_field(
                tag, key, field, aliases, left, right)
            if failure:
                errors.append(
                    f"tag={tag} key={key} field={field} np1={failure.get('np1')} "
                    f"np{T10_NP}={failure.get('npN')} abs_diff={failure.get('abs_diff')} "
                    f"rel_diff={failure.get('rel_diff')} error={failure.get('error')}"
                )
                if len(errors) >= 8:
                    break
            elif compared:
                compared_fields.add(field)
                update_diag_stats(stats, tag, key, compared)
        if len(errors) >= 8:
            break
        for field in categorical_fields:
            failure = compare_categorical_diagnostic_field(tag, key, field, left, right)
            if failure:
                errors.append(
                    f"tag={tag} key={key} field={field} np1={failure.get('np1')} "
                    f"np{T10_NP}={failure.get('npN')} error={failure.get('error')}"
                )
                if len(errors) >= 8:
                    break
            elif field in left or field in right:
                compared_fields.add(field)
        if len(errors) >= 8:
            break
    return errors, stats, compared_fields, len(left_map), len(right_map)


def validate_d0b_t10(np1_case, npn_case, config_name):
    name = "D0b-lite T10 MPI diagnostic comparison"
    run_errors = []
    if np1_case.timed_out or np1_case.rc != 0:
        run_errors.append(f"np=1 {config_name} run rc={np1_case.rc} timed_out={np1_case.timed_out}")
    if npn_case.timed_out or npn_case.rc != 0:
        run_errors.append(f"np={T10_NP} {config_name} run rc={npn_case.rc} timed_out={npn_case.timed_out}")
    if run_errors:
        return fail_result(
            name,
            run_errors[0],
            {
                "errors": run_errors,
                "config": config_name,
                "np1_log": str(np1_case.log),
                "npN_log": str(npn_case.log),
            },
            str(npn_case.log),
        )
    np1_bhaccr = bhaccr_rows(np1_case.log)
    npn_bhaccr = bhaccr_rows(npn_case.log)
    np1_bhfdbk = bhfdbk_rows(np1_case.log)
    npn_bhfdbk = bhfdbk_rows(npn_case.log)
    errors = []
    stats = {}
    compared_fields = set()
    np1_burst_keys, np1_burst_summaries = gas_present_burst_keys(np1_bhfdbk)
    npn_burst_keys, npn_burst_summaries = gas_present_burst_keys(npn_bhfdbk)
    matched_burst_keys = sorted(np1_burst_keys & npn_burst_keys, key=str)
    if not matched_burst_keys:
        errors.append(
            "T10 runs did not include matching gas-present burst rows "
            f"np1={len(np1_burst_keys)} np{T10_NP}={len(npn_burst_keys)}"
        )
    for tag, left_rows, right_rows, numeric, categorical in (
        ("BHACCR", np1_bhaccr, npn_bhaccr,
         T10_BHACCR_NUMERIC_FIELDS, T10_BHACCR_CATEGORICAL_FIELDS),
        ("BHFDBK", np1_bhfdbk, npn_bhfdbk,
         T10_BHFDBK_NUMERIC_FIELDS, T10_BHFDBK_CATEGORICAL_FIELDS),
    ):
        tag_errors, tag_stats, tag_fields, left_count, right_count = compare_diagnostic_sets(
            tag, left_rows, right_rows, numeric, categorical)
        errors.extend(tag_errors)
        stats.update(tag_stats)
        compared_fields.update(f"{tag}.{field}" for field in tag_fields)
        if left_count == 0 or right_count == 0:
            errors.append(f"{tag} matched diagnostics absent for config={config_name}")
    metrics = {
        "config": config_name,
        "np1_log": str(np1_case.log),
        "npN_log": str(npn_case.log),
        "np1_bhaccr_rows": len(np1_bhaccr),
        "npN_bhaccr_rows": len(npn_bhaccr),
        "np1_bhfdbk_rows": len(np1_bhfdbk),
        "npN_bhfdbk_rows": len(npn_bhfdbk),
        "matched_gas_present_burst_keys": [str(key) for key in matched_burst_keys],
        "np1_gas_present_burst_rows": np1_burst_summaries,
        "npN_gas_present_burst_rows": npn_burst_summaries,
        "fields_compared": sorted(compared_fields),
        "field_differences": stats,
        "rtol": T10_SHORT_RTOL,
        "atol": T10_ABS_TOL,
    }
    if errors:
        return fail_result(name, errors[0], {"errors": errors[:8], **metrics}, str(npn_case.log))
    max_abs = max((item["max_abs_diff"] for item in stats.values()), default=0.0)
    max_rel = max((item["max_rel_diff"] for item in stats.values()), default=0.0)
    return pass_result(
        name,
        f"{config_name} np=1 vs np={T10_NP}; max_abs_diff={max_abs:.8e} "
        f"max_rel_diff={max_rel:.8e}",
        metrics,
        str(npn_case.log),
    )


def run_t10_cases():
    np1 = run_named("D0b_T10_np1_burst", d0b_burst_overrides())
    npn_dir = prepare_case("D0b_T10_npN_burst", overrides=d0b_burst_overrides())
    npn = run_enzo(npn_dir, np_ranks=T10_NP)
    if np1.rc == 0 and npn.rc == 0:
        return np1, npn, "burst"
    if mpi_environment_failure(npn.log):
        fallback_np1 = run_named("D0b_T10_np1_no_burst", d0b_no_burst_overrides())
        fallback_npn_dir = prepare_case(
            "D0b_T10_npN_no_burst",
            overrides=d0b_no_burst_overrides(),
        )
        fallback_npn = run_enzo(fallback_npn_dir, np_ranks=T10_NP)
        return fallback_np1, fallback_npn, "no-burst fallback after MPI/environment burst failure"
    return np1, npn, "burst"


def validate_expected_failure(name, case, expected_substring):
    text = case.log.read_text(errors="ignore") if case.log.exists() else ""
    if case.timed_out:
        return fail_result(name, f"run timed out after {CASE_TIMEOUT}s", log=str(case.log))
    if case.rc == 0:
        return fail_result(
            name,
            "run exited rc=0 but expected nonzero validation failure",
            log=str(case.log),
        )
    if expected_substring not in text:
        return fail_result(
            name,
            f"run exited rc={case.rc} but missing expected message: "
            f"{expected_substring!r}",
            log=str(case.log),
        )
    return pass_result(
        name,
        f"run exited rc={case.rc} and emitted expected message",
        {"expected_message": expected_substring},
        str(case.log),
    )


def take_hygiene_snapshot():
    branch = git_output("branch", "--show-current").strip()
    head = git_output("rev-parse", "HEAD").strip()
    origin_main = git_output("rev-parse", "origin/main").strip()
    status_all = git_output("status", "--short")
    stash = git_output("stash", "list")
    diff_src = git_output("diff", "--", "src/enzo/")
    diff_uuid = git_output("diff", "--", "src/enzo/uuid/gen_uuid.c")
    protected_status = git_output("status", "--short", "--", *PROTECTED_PATHS)
    return {
        "branch": branch,
        "head": head,
        "origin_main": origin_main,
        "d0b_baseline_in_head": git_is_ancestor(D0B_BASELINE, "HEAD"),
        "d0a_commit_in_head": git_is_ancestor(D0A_COMMIT, "HEAD"),
        "status_all": status_all,
        "stash": stash,
        "diff_src": diff_src,
        "diff_uuid": diff_uuid,
        "protected_status": protected_status,
    }


def validate_h1(start_snapshot):
    name = "H1 repository hygiene"
    end_snapshot = take_hygiene_snapshot()
    errors = []
    if end_snapshot["diff_src"] != start_snapshot["diff_src"]:
        errors.append("git diff -- src/enzo/ changed during harness run")
    if end_snapshot["diff_uuid"].strip():
        errors.append("git diff -- src/enzo/uuid/gen_uuid.c is non-empty")
    if PRESERVED_STASH_SUBSTRING not in end_snapshot["stash"]:
        errors.append(f"stash missing preserved entry containing {PRESERVED_STASH_SUBSTRING!r}")
    if end_snapshot["stash"] != start_snapshot["stash"]:
        errors.append("git stash list changed during harness run")
    if end_snapshot["protected_status"] != start_snapshot["protected_status"]:
        errors.append("protected-path git status changed during harness run")
    metrics = {"start": start_snapshot, "end": end_snapshot}
    if errors:
        return fail_result(name, errors[0], {"errors": errors, **metrics})
    return pass_result(
        name,
        "src/enzo diff unchanged; uuid diff empty; preserved stash present; "
        "protected-path statuses unchanged",
        metrics,
    )


def preflight_results():
    missing = []
    if not TS3_PARAM.exists():
        missing.append(f"TS3 parameter file missing: {TS3_PARAM}")
    if not ENZO_EXE.exists():
        missing.append(f"ENZO_EXE missing: {ENZO_EXE}")
    elif not os.access(ENZO_EXE, os.X_OK):
        missing.append(f"ENZO_EXE is not executable: {ENZO_EXE}")
    if not executable_available(MPIRUN):
        missing.append(f"MPIRUN unavailable: {MPIRUN}")
    return missing


def run_optional_phaseb_lite(include_phaseb):
    name = "T12-lite existing Phase B matrix compatibility"
    if not include_phaseb:
        return skip_result(
            name,
            "set --include-phaseb-lite or BH_FEEDBACK_PHASEC_RUN_PHASEB=1 to run",
        )
    script = ROOT / "run/BHSeed/feedback_phaseB_tests/run_feedback_phaseB_matrix.py"
    if not script.exists():
        return skip_result(name, f"Phase B matrix missing: {script}")
    missing = preflight_results()
    if missing:
        return skip_result(name, "; ".join(missing))

    phaseb_out = OUT_ROOT / "T12_lite_phaseB_matrix"
    env = os.environ.copy()
    env["BH_FEEDBACK_PHASEB_OUT_ROOT"] = str(phaseb_out)
    if DEFAULT_CONDA_BIN.exists():
        env["PATH"] = f"{DEFAULT_CONDA_BIN}:{env.get('PATH', '')}"
    log_path = OUT_ROOT / "T12_lite_phaseB_matrix.log"
    with open(log_path, "w") as log:
        try:
            proc = subprocess.run(
                [sys.executable, str(script)],
                cwd=ROOT,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                timeout=PHASEB_TIMEOUT,
            )
        except subprocess.TimeoutExpired:
            log.write(f"\nTIMEOUT after {PHASEB_TIMEOUT} seconds\n")
            return skip_result(
                name,
                f"Phase B matrix timed out after {PHASEB_TIMEOUT}s",
                log=str(log_path),
            )
    text = log_path.read_text(errors="ignore")
    environment_failures = (
        "ImportError",
        "No such file or directory",
        "No module named",
        "orterun was unable",
        "mpirun was unable",
    )
    if proc.returncode == 0 and re.search(r"SUMMARY .*failed=0", text):
        return pass_result(name, "Phase B matrix reported failed=0", log=str(log_path))
    if any(token in text for token in environment_failures):
        return skip_result(
            name,
            "environment prevented Phase B matrix run",
            {"rc": proc.returncode},
            str(log_path),
        )
    return optional_fail_result(
        name,
        f"Phase B matrix ran but returned rc={proc.returncode} or reported failures",
        {"rc": proc.returncode},
        str(log_path),
    )


def print_repo_snapshot(snapshot):
    print("Repository snapshot:")
    print(f"  root={ROOT}")
    print(f"  branch={snapshot['branch']}")
    print(f"  HEAD={snapshot['head']}")
    print(f"  origin/main={snapshot['origin_main']}")
    print(f"  d0b_baseline={D0B_BASELINE} in_HEAD={snapshot['d0b_baseline_in_head']}")
    print(f"  d0a_commit={D0A_COMMIT} in_HEAD={snapshot['d0a_commit_in_head']}")
    print("  protected_status:")
    status = snapshot["protected_status"].strip()
    print(f"    {status if status else '<clean>'}")
    print("  stash:")
    stash = snapshot["stash"].strip()
    print(f"    {stash if stash else '<empty>'}")


def print_environment():
    print("Environment:")
    print(f"  OUT_ROOT={OUT_ROOT}")
    print(f"  TS3_PARAM={TS3_PARAM}")
    print(f"  ENZO_EXE={ENZO_EXE}")
    print(f"  MPIRUN={MPIRUN}")
    print(f"  T10_NP={T10_NP}")
    print(f"  H5LS={H5LS}")
    print(f"  H5DUMP={H5DUMP}")
    print(f"  CASE_TIMEOUT={CASE_TIMEOUT}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--include-phaseb-lite",
        action="store_true",
        help="run optional T12-lite Phase B matrix compatibility check",
    )
    args = parser.parse_args(argv)
    include_phaseb = args.include_phaseb_lite or (
        os.environ.get("BH_FEEDBACK_PHASEC_RUN_PHASEB") == "1"
    )

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    start_snapshot = take_hygiene_snapshot()
    print_repo_snapshot(start_snapshot)
    print_environment()

    results = []
    if not start_snapshot["d0b_baseline_in_head"]:
        results.append(fail_result(
            "Repository D0b baseline ancestry guard",
            f"HEAD={start_snapshot['head']} does not contain D0b baseline {D0B_BASELINE}",
        ))
    else:
        results.append(pass_result(
            "Repository D0b baseline ancestry guard",
            f"HEAD contains D0b baseline {D0B_BASELINE}",
        ))
    if not start_snapshot["d0a_commit_in_head"]:
        results.append(fail_result(
            "Repository D0a ancestry guard",
            f"HEAD={start_snapshot['head']} does not contain D0a commit {D0A_COMMIT}",
        ))
    else:
        results.append(pass_result(
            "Repository D0a ancestry guard",
            f"HEAD contains D0a commit {D0A_COMMIT}",
        ))
    results.append(pass_result(
        "Repository origin/main observation",
        f"origin/main currently resolves to {start_snapshot['origin_main']}",
    ))
    if PRESERVED_STASH_SUBSTRING not in start_snapshot["stash"]:
        results.append(fail_result(
            "Repository preserved stash guard",
            f"stash missing entry containing {PRESERVED_STASH_SUBSTRING!r}",
        ))
    else:
        results.append(pass_result(
            "Repository preserved stash guard",
            "preserved gen_uuid.c stash present at script start",
        ))

    results.append(validate_d2_1_static_source())

    missing = preflight_results()
    if missing:
        reason = "; ".join(missing)
        for name in (
            "T1 default Phase C smoke",
            "T2 BHFeedbackEddingtonFactor = 0.5",
            "T3 BHFeedbackEddingtonFactor = 2.0",
            "T4 cumulative reservoir counters",
            "T5 conservation_residual diagnostic",
            "T6 HDF5 label regression guard",
            "D0b-lite F1 gas-present burst",
            "D0b-lite F3 conservation audit through burst",
            "D0b-lite T7 restart preservation",
            "D0b-lite T10 MPI diagnostic comparison",
            "D0c-acc A1 Bondi formula/unit",
            "D0c-acc A2 Eddington cap activation",
            "D0c-acc A5 gas-limited mismatch",
            "D0c-acc F6 mismatch metric",
            "D2-1 realized mdot accretion diagnostics",
            "D2-1 feedback diagnostic read keeps requested basis",
            "D0c-seed S0 default-off no-seed",
            "D0c-seed S1 explicit opt-in seed creation",
            "D0c-seed S2 temperature gate no-seed",
            "D0c-seed S3 deterministic np1 repeat",
            "D0c-repos R0 newly seeded skip diagnostic",
            "D0c-repos R1 method0 diagnostic no movement",
            "D0c-repos R2 method2 safety rejection",
            "D0c-repos R4 actual movement ON/OFF",
            "T8 negative BHFeedbackEddingtonFactor validation",
            "T9 efficiency product > 1.0 validation",
            "D2-0 config guard MBHAccretion/BHAccretionMethod mutual exclusion",
        ):
            results.append(fail_result(name, f"preflight failed: {reason}"))
        results.append(skip_result(
            "D0c-seed S4 optional MPI seed consistency",
            f"preflight unavailable: {reason}",
        ))
        results.append(skip_result(
            "D0c-repos R3 optional MPI BHREPOS consistency",
            f"preflight unavailable: {reason}",
        ))
    else:
        t1 = run_named(
            "T1_default_phaseC_smoke",
            baseline_overrides(),
            remove_keys=("BHFeedbackEddingtonFactor",),
        )
        results.append(validate_t1(t1))
        t1_bhaccr = bhaccr_rows(t1.log) if t1.rc == 0 else []

        t2 = run_named(
            "T2_edd_factor_0p5",
            baseline_overrides(BHFeedbackEddingtonFactor="0.5"),
        )
        results.append(validate_edd_factor_case(
            "T2 BHFeedbackEddingtonFactor = 0.5",
            t2,
            0.5,
            t1_bhaccr,
        ))

        t3 = run_named(
            "T3_edd_factor_2p0",
            baseline_overrides(BHFeedbackEddingtonFactor="2.0"),
        )
        results.append(validate_edd_factor_case(
            "T3 BHFeedbackEddingtonFactor = 2.0",
            t3,
            2.0,
            t1_bhaccr,
        ))

        results.append(validate_t4(t1))
        results.append(validate_t5(t1))
        results.append(validate_t6(t1))

        f1 = run_named(
            "D0b_F1_F3_gas_present_burst",
            d0b_burst_overrides(),
        )
        f1_result = validate_d0b_f1(f1)
        results.append(f1_result)
        results.append(validate_d0b_f3(f1, f1_result))

        t8 = run_named(
            "T8_negative_edd_factor",
            baseline_overrides(BHFeedbackEddingtonFactor="-1.0"),
        )
        results.append(validate_expected_failure(
            "T8 negative BHFeedbackEddingtonFactor validation",
            t8,
            EDD_FACTOR_FATAL,
        ))

        t9 = run_named(
            "T9_efficiency_product_gt_one",
            baseline_overrides(
                BHAccretionMethod=0,
                BHAccretionRadiativeEfficiency="2.0",
                BHFeedbackThermalEfficiency="0.6",
            ),
            remove_keys=("BHFeedbackEddingtonFactor",),
        )
        results.append(validate_expected_failure(
            "T9 efficiency product > 1.0 validation",
            t9,
            EFFICIENCY_PRODUCT_FATAL,
        ))

        d2_guard = run_named(
            "D2_0_config_guard_dual_mbh_bhaccretion",
            baseline_overrides(
                MBHAccretion=1,
                BHAccretionMethod=1,
            ),
            remove_keys=("BHFeedbackEddingtonFactor",),
        )
        results.append(validate_expected_failure(
            "D2-0 config guard MBHAccretion/BHAccretionMethod mutual exclusion",
            d2_guard,
            DUAL_MBH_BHACCRETION_FATAL,
        ))

        t7_continuous, t7_first, t7_second, t7_restart_param = run_t7_restart_cases()
        results.append(validate_d0b_t7(
            t7_continuous,
            t7_first,
            t7_second,
            t7_restart_param,
        ))

        t10_np1, t10_npn, t10_config = run_t10_cases()
        results.append(validate_d0b_t10(t10_np1, t10_npn, t10_config))

        a1 = run_named("D0c_A1_bondi_unit", d0c_acc_base_overrides())
        results.append(validate_d0c_a1(a1))

        a2 = run_named(
            "D0c_A2_eddington_cap",
            d0c_acc_base_overrides(
                BHSeedMass="1e8",
                BHSeedMinEnclosedMass="1e8",
                BHFeedbackEddingtonFactor="1.0",
            ),
        )
        results.append(validate_d0c_a2(a2))

        a5 = run_named(
            "D0c_A5_gas_limited_mismatch",
            d0c_acc_base_overrides(
                BHSeedMass="1e10",
                BHSeedMinEnclosedMass="1e10",
                BHFeedbackEddingtonFactor="50000",
                BHAccretionRemovalMode=1,
                TestStarParticleDensity=51,
            ),
        )
        results.append(validate_d0c_a5(a5))
        results.append(validate_d0c_f6(a5))
        results.append(validate_d2_1_realized_rate(a5))
        results.append(validate_d2_1_feedback_diagnostics(a5))

        # The TS3_wrap initializer does not plant overdense candidate cells when
        # BHSeedingMethod is absent/default-off, so this validates the
        # master-switch early-exit path.
        s0 = run_named(
            "D0c_seed_S0_default_off_no_seed",
            {
                "BHAccretionMethod": 0,
                "BHFeedbackMethod": 0,
            },
            remove_keys=("BHSeedingMethod",),
        )
        results.append(validate_d0c_seed_s0(s0))

        s1 = run_named(
            "D0c_seed_S1_explicit_opt_in_seed_creation",
            d0c_seed_opt_in_overrides(),
        )
        s1_result = validate_d0c_seed_s1(s1)
        results.append(s1_result)

        s2 = run_named(
            "D0c_seed_S2_temperature_gate_no_seed",
            d0c_seed_opt_in_overrides(BHSeedTemperatureThreshold="100.0"),
        )
        results.append(validate_d0c_seed_s2(s2))

        # This is a fixture-level determinism check, not a production-level
        # stochastic/AMR/MPI determinism proof.
        s3_a = run_named(
            "D0c_seed_S3_deterministic_np1_repeat_a",
            d0c_seed_opt_in_overrides(),
        )
        s3_b = run_named(
            "D0c_seed_S3_deterministic_np1_repeat_b",
            d0c_seed_opt_in_overrides(),
        )
        results.append(validate_d0c_seed_s3(s3_a, s3_b))
        results.append(run_d0c_seed_s4(s1_result))

        repos_r0_r1 = run_named(
            "D0c_repos_R0_R1_method0_diagnostic",
            d0c_repos_base_overrides(BHRepositionMethod=0, BHRepositionVerbose=1),
        )
        results.append(validate_d0c_repos_r0(repos_r0_r1))
        results.append(validate_d0c_repos_r1(repos_r0_r1))

        repos_r2 = run_named(
            "D0c_repos_R2_method2_safety_rejection",
            d0c_repos_base_overrides(BHRepositionMethod=2, BHRepositionVerbose=1),
        )
        results.append(validate_d0c_repos_r2(repos_r2))
        results.append(run_d0c_repos_r3(repos_r0_r1))
        # R4 uses a scratch-only copied mini-restart derived from the R0/R1
        # TS3 output. The HDF5 edits place a non-new BH in the interior and
        # inject one offset dense cell so method 0 can see a target without
        # moving while method 2 moves to the logged active peak.
        results.append(run_d0c_repos_r4(repos_r0_r1))

    results.append(run_optional_phaseb_lite(include_phaseb))
    results.append(validate_h1(start_snapshot))

    counts = {}
    for result in results:
        counts[result.status] = counts.get(result.status, 0) + 1
        log_note = f" log={result.log}" if result.log else ""
        print(f"{result.status} {result.name}: {result.details}{log_note}")

    summary = {
        "out_root": str(OUT_ROOT),
        "counts": counts,
        "results": [result.__dict__ for result in results],
    }
    (OUT_ROOT / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    print(
        "SUMMARY "
        f"PASS={counts.get('PASS', 0)} "
        f"FAIL={counts.get('FAIL', 0)} "
        f"SKIP={counts.get('SKIP', 0)} "
        f"OPTIONAL_FAIL={counts.get('OPTIONAL_FAIL', 0)} "
        f"out={OUT_ROOT}"
    )
    return 1 if counts.get("FAIL", 0) else 0


if __name__ == "__main__":
    raise SystemExit(main())
