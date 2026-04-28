#!/usr/bin/env python3
"""Conservative Phase C regression harness for Enzo BH feedback D0a.

This script is validation-only. It copies TS3_wrap parameter files into an
output directory, applies per-case overrides there, and never edits source
inputs in place.

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

EXPECTED_HEAD = "05483a26e0a5628be92b68ce64b17290a7914691"
PHASEC_TIP = "41d3eda7e1acb8b44ea3d89b5818317467ced520"
PRESERVED_STASH_SUBSTRING = "Preserve gen_uuid.c outside Phase B merge"

EDD_FACTOR_FATAL = "BHFeedbackEddingtonFactor must be >= 0"
EFFICIENCY_PRODUCT_FATAL = (
    "Product BHAccretionRadiativeEfficiency * "
    "BHFeedbackThermalEfficiency exceeds 1.0"
)

HDF5_EXPECTED_LABELS = (
    "bh_cumul_reservoir_in",
    "bh_cumul_reservoir_out",
)
HDF5_STALE_LABELS = (
    "particle_attribute_20",
    "particle_attribute_21",
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


def prepare_case(name, overrides=None, remove_keys=()):
    overrides = overrides or {}
    cdir = OUT_ROOT / name
    if cdir.exists():
        shutil.rmtree(cdir)
    cdir.mkdir(parents=True, exist_ok=True)
    text = make_text(TS3_PARAM, overrides, remove_keys=remove_keys)
    (cdir / "case.enzo").write_text(text)
    return cdir


def run_enzo(cdir, np_ranks=1, log_name="run.log"):
    env = os.environ.copy()
    env["HDF5_DISABLE_VERSION_CHECK"] = "2"
    if DEFAULT_CONDA_BIN.exists():
        env["PATH"] = f"{DEFAULT_CONDA_BIN}:{env.get('PATH', '')}"
    cmd = [MPIRUN, "-n", str(np_ranks), str(ENZO_EXE), "-d", "case.enzo"]
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
    merge_base_phasec = git_output("merge-base", "HEAD", PHASEC_TIP).strip()
    status_all = git_output("status", "--short")
    stash = git_output("stash", "list")
    diff_src = git_output("diff", "--", "src/enzo/")
    diff_uuid = git_output("diff", "--", "src/enzo/uuid/gen_uuid.c")
    protected_status = git_output("status", "--short", "--", *PROTECTED_PATHS)
    return {
        "branch": branch,
        "head": head,
        "origin_main": origin_main,
        "merge_base_phasec": merge_base_phasec,
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
    if end_snapshot["diff_src"].strip():
        errors.append("git diff -- src/enzo/ is non-empty")
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
        "src/enzo diffs empty; uuid diff empty; preserved stash present; "
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
    print(f"  phaseC_tip_merge_base={snapshot['merge_base_phasec']}")
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
    if start_snapshot["head"] != EXPECTED_HEAD:
        results.append(fail_result(
            "Repository HEAD guard",
            f"HEAD={start_snapshot['head']} expected {EXPECTED_HEAD}",
        ))
    else:
        results.append(pass_result(
            "Repository HEAD guard",
            f"HEAD matches expected Phase C merge {EXPECTED_HEAD}",
        ))
    if start_snapshot["origin_main"] != EXPECTED_HEAD:
        results.append(fail_result(
            "Repository origin/main guard",
            f"origin/main={start_snapshot['origin_main']} expected {EXPECTED_HEAD}",
        ))
    else:
        results.append(pass_result(
            "Repository origin/main guard",
            f"origin/main matches expected Phase C merge {EXPECTED_HEAD}",
        ))
    if start_snapshot["merge_base_phasec"] != PHASEC_TIP:
        results.append(fail_result(
            "Repository Phase C ancestry guard",
            f"merge-base HEAD {PHASEC_TIP}={start_snapshot['merge_base_phasec']}",
        ))
    else:
        results.append(pass_result(
            "Repository Phase C ancestry guard",
            f"HEAD contains Phase C tip {PHASEC_TIP}",
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
            "T8 negative BHFeedbackEddingtonFactor validation",
            "T9 efficiency product > 1.0 validation",
        ):
            results.append(fail_result(name, f"preflight failed: {reason}"))
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
