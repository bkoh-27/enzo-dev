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

ROOT = Path(__file__).resolve().parents[3]
TS3_DIR = ROOT / "run/BHSeed/TS3_wrap"
OUT_ROOT = ROOT / f"run/BHSeed/feedback_phaseA_tests/matrix_{time.strftime('%Y%m%d_%H%M%S')}"
ENZO_EXE = ROOT / "src/enzo/enzo.exe"

MODULE_INIT = "source /etc/profile.d/modules.sh && module load openmpi/4.1.8"
LD_PATH = ":".join([
    "/mnt/sw/nix/store/8r0bpn8x0ljxndwzf7jjkp15ir96qpfg-hdf5-1.10.11/lib",
    "/mnt/sw/nix/store/cprcg9nf8xpr7gljz5d3vgg97blyzwfx-hdf5-1.12.1/lib",
    "/mnt/home/boh10/bh_proj/local/lib",
    "/mnt/home/boh10/ltu_proj/local/lib",
])

# TS3_wrap units (non-comoving)
CLIGHT = 2.99792458e10
MSUN = 1.989e33
YR_S = 3.1557e7


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


def parse_kv_line(line: str):
    data = {}
    for tok in line.strip().split():
        if "=" not in tok:
            continue
        key, val = tok.split("=", 1)
        if key.startswith("["):
            continue
        try:
            if val.lower().startswith("nan"):
                data[key] = float("nan")
            elif val.lower().startswith("inf"):
                data[key] = float("inf")
            elif re.match(r"^-?\\d+$", val):
                data[key] = int(val)
            else:
                data[key] = float(val)
        except Exception:
            data[key] = val
    return data


def read_lines(log_path: Path, tag: str):
    lines = []
    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if tag in line:
                lines.append(line.strip())
    return lines


def run_enzo(case_dir: Path, np_ranks: int, log_name: str, param_name: str):
    cmd = (
        f"cd {case_dir} && "
        f"{MODULE_INIT} && "
        f"export LD_LIBRARY_PATH={LD_PATH}:$LD_LIBRARY_PATH && "
        f"mpirun -np {np_ranks} {ENZO_EXE} {param_name} > {log_name} 2>&1"
    )
    proc = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
    return proc.returncode, case_dir / log_name, proc.stdout + proc.stderr


def prepare_case(case_name: str, param_text: str):
    case_dir = OUT_ROOT / case_name
    ensure_clean_dir(case_dir)

    # Copy TS3_wrap IC and helper script for reproducibility
    shutil.copytree(TS3_DIR / "DD0000", case_dir / "DD0000")
    shutil.copy2(TS3_DIR / "make_ts3wrap_ic.py", case_dir / "make_ts3wrap_ic.py")

    (case_dir / "case.enzo").write_text(param_text)
    return case_dir


def make_case_text(base_file: Path, overrides: dict):
    text = base_file.read_text()
    for key, value in overrides.items():
        text = set_param(text, key, value)
    return text


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

def test0_ts3_wrap_baseline():
    base_text = (TS3_DIR / "bhseed_ts3wrap.enzo").read_text()
    case_dir = prepare_case("case0_ts3_wrap_baseline", base_text)

    rc1, log1, out1 = run_enzo(case_dir, 1, "ts3wrap_np1.log", "case.enzo")
    rc4, log4, out4 = run_enzo(case_dir, 4, "ts3wrap_np4.log", "case.enzo")
    if rc1 != 0 or rc4 != 0:
        return TestResult("0_ts3_wrap_baseline", False, f"run failed rc1={rc1} rc4={rc4}", {})

    def check_log(log_path):
        lines = read_lines(log_path, "[BHSEED]")
        if len(lines) < 2:
            return False, f"missing BHSEED lines in {log_path.name}"
        step1 = lines[0]
        step2 = lines[1]
        ok = (
            "created=1" in step1 and "ncand_global=2" in step1 and
            "created=0" in step2 and "dist_blocked=2" in step2
        )
        return ok, f"step1={step1} step2={step2}"

    ok1, d1 = check_log(log1)
    ok4, d4 = check_log(log4)
    ok = ok1 and ok4
    details = f"np1: {d1}; np4: {d4}"
    return TestResult("0_ts3_wrap_baseline", ok, details, {})


def test_caseA():
    overrides = {
        "StopTime": "0.06",
        "BHAccretionMethod": "1",
        "BHAccretionRunEveryTimestep": "1",
        "BHAccretionVerbose": "1",
        "BHFeedbackMethod": "1",
        "BHFeedbackModeThreshold": "0.01",
        "BHFeedbackKernelRadius": "1.0",
        "BHFeedbackThermalEfficiency": "0.02",
        "BHFeedbackMinEnergyBurst": "1e50",
        "BHFeedbackKineticEfficiency": "0.1",
        "BHFeedbackWindVelocity": "1e4",
        "BHFeedbackKineticGeometry": "0",
        "BHFeedbackVerbose": "1",
    }
    text = make_case_text(TS3_DIR / "bhseed_ts3wrap.enzo", overrides)
    case_dir = prepare_case("caseA_feedback_baseline", text)
    rc, log_path, out = run_enzo(case_dir, 1, "np1.log", "case.enzo")
    if rc != 0:
        return None, TestResult("caseA", False, f"run failed rc={rc}", {})
    return log_path, None


def test_caseB():
    overrides = {
        "StopTime": "0.06",
        "BHSeedMass": "1e6",
        "BHAccretionMethod": "1",
        "BHAccretionRunEveryTimestep": "1",
        "BHAccretionVerbose": "1",
        "BHFeedbackMethod": "1",
        "BHFeedbackModeThreshold": "0.01",
        "BHFeedbackKernelRadius": "1.0",
        "BHFeedbackThermalEfficiency": "0.02",
        "BHFeedbackMinEnergyBurst": "1e50",
        "BHFeedbackKineticEfficiency": "0.1",
        "BHFeedbackWindVelocity": "1e4",
        "BHFeedbackKineticGeometry": "0",
        "BHFeedbackVerbose": "1",
    }
    text = make_case_text(TS3_DIR / "bhseed_ts3wrap.enzo", overrides)
    case_dir = prepare_case("caseB_feedback_thermal", text)
    rc, log_path, out = run_enzo(case_dir, 1, "np1.log", "case.enzo")
    if rc != 0:
        return None, TestResult("caseB", False, f"run failed rc={rc}", {})
    return log_path, None


def test_caseC():
    overrides = {
        "StopTime": "0.06",
        "BHAccretionMethod": "1",
        "BHAccretionRunEveryTimestep": "1",
        "BHAccretionVerbose": "1",
        "BHFeedbackMethod": "1",
        "BHFeedbackModeThreshold": "0.01",
        "BHFeedbackKernelRadius": "1.0",
        "BHFeedbackThermalEfficiency": "0.02",
        "BHFeedbackMinEnergyBurst": "1e41",
        "BHFeedbackKineticEfficiency": "0.1",
        "BHFeedbackWindVelocity": "1e4",
        "BHFeedbackKineticGeometry": "0",
        "BHFeedbackVerbose": "1",
    }
    text = make_case_text(TS3_DIR / "bhseed_ts3wrap.enzo", overrides)
    case_dir = prepare_case("caseC_feedback_burst", text)
    rc, log_path, out = run_enzo(case_dir, 1, "np1.log", "case.enzo")
    if rc != 0:
        return None, TestResult("caseC", False, f"run failed rc={rc}", {})
    return log_path, None


def test_caseD_multi():
    text = (TS3_DIR / "bhseed_ts3wrap_feedback_multi.enzo").read_text()
    case_dir = prepare_case("caseD_feedback_multi_bh", text)
    rc, log_path, out = run_enzo(case_dir, 1, "np1.log", "case.enzo")
    if rc != 0:
        return None, TestResult("caseD", False, f"run failed rc={rc}", {})
    return log_path, None


def test_caseA_np4():
    overrides = {
        "StopTime": "0.06",
        "BHAccretionMethod": "1",
        "BHAccretionRunEveryTimestep": "1",
        "BHAccretionVerbose": "1",
        "BHFeedbackMethod": "1",
        "BHFeedbackModeThreshold": "0.01",
        "BHFeedbackKernelRadius": "1.0",
        "BHFeedbackThermalEfficiency": "0.02",
        "BHFeedbackMinEnergyBurst": "1e50",
        "BHFeedbackKineticEfficiency": "0.1",
        "BHFeedbackWindVelocity": "1e4",
        "BHFeedbackKineticGeometry": "0",
        "BHFeedbackVerbose": "1",
    }
    text = make_case_text(TS3_DIR / "bhseed_ts3wrap.enzo", overrides)
    case_dir = prepare_case("caseA_feedback_baseline_np4", text)
    rc, log_path, out = run_enzo(case_dir, 4, "np4.log", "case.enzo")
    if rc != 0:
        return None, TestResult("caseA_np4", False, f"run failed rc={rc}", {})
    return log_path, None


def compute_l_feedback_from_mdot(mdot_msun_per_yr, epsilon_r):
    mdot_cgs = mdot_msun_per_yr * (MSUN / YR_S)
    return epsilon_r * mdot_cgs * (CLIGHT ** 2)


def get_last_bhfdbk(log_path: Path, bh_id: int = 1):
    lines = [l for l in read_lines(log_path, "[BHFDBK]")]
    if not lines:
        raise RuntimeError("no BHFDBK lines")
    parsed = [parse_kv_line(l) for l in lines]
    chosen = None
    for d, l in zip(parsed, lines):
        if int(d.get("bh_id", -99999)) == bh_id:
            chosen = (d, l)
    if chosen is None:
        chosen = (parsed[-1], lines[-1])
    return chosen[0], lines


def get_last_bhaccr(log_path: Path, bh_id: int = 1):
    lines = read_lines(log_path, "[BHACCR]")
    if not lines:
        raise RuntimeError("no BHACCR lines")
    parsed = [parse_kv_line(l) for l in lines]
    chosen = None
    for d, l in zip(parsed, lines):
        if int(d.get("bh_id", -99999)) == bh_id:
            chosen = (d, l)
    if chosen is None:
        chosen = (parsed[-1], lines[-1])
    return chosen[0], lines


def main():
    results = []

    # Test 0: TS3_wrap baseline
    results.append(test0_ts3_wrap_baseline())

    # Case A baseline
    logA, errA = test_caseA()
    if errA:
        results.append(errA)
        logA = None

    # Case B thermal
    logB, errB = test_caseB()
    if errB:
        results.append(errB)
        logB = None

    # Case C burst
    logC, errC = test_caseC()
    if errC:
        results.append(errC)
        logC = None

    # Case D multi-BH
    logD, errD = test_caseD_multi()
    if errD:
        results.append(errD)
        logD = None

    # Case A np4
    logA4, errA4 = test_caseA_np4()
    if errA4:
        results.append(errA4)
        logA4 = None

    # Test 1: feedback luminosity analytic
    if logA:
        bh, _ = get_last_bhfdbk(logA, bh_id=1)
        bhaccr, _ = get_last_bhaccr(logA, bh_id=1)
        epsilon_r = 0.1
        mdot = bhaccr.get("Mdot_actual", 0.0)
        l_pred = compute_l_feedback_from_mdot(mdot, epsilon_r)
        l_fb = bh.get("L_feedback", 0.0)
        rel_err = abs(l_fb - l_pred) / max(abs(l_pred), 1.0)
        ok = rel_err < 5e-3
        results.append(TestResult("1", ok, f"rel_err={rel_err:.3e}", {"rel_err": rel_err}))

    # Test 2: mode classification thermal
    if logB:
        bh, _ = get_last_bhfdbk(logB, bh_id=1)
        mode = bh.get("feedback_mode", "")
        ok = (mode == "THERMAL")
        results.append(TestResult("2", ok, f"mode={mode} f_Edd={bh.get('f_Edd', -1):.8e}", {}))

    # Test 3: mode classification kinetic
    if logA:
        bh, _ = get_last_bhfdbk(logA, bh_id=1)
        mode = bh.get("feedback_mode", "")
        ok = (mode == "KINETIC")
        results.append(TestResult("3", ok, f"mode={mode} f_Edd={bh.get('f_Edd', -1):.8e}", {}))

    # Test 4: reservoir accumulation diagnostic
    if logA:
        lines = read_lines(logA, "[BHFDBK]")
        parsed = [parse_kv_line(l) for l in lines]
        seq = [d for d in parsed if int(d.get("bh_id", -1)) == 1]
        ra = [d.get("reservoir_after_accum", d.get("reservoir_after", 0.0)) for d in seq[0:3]]
        rb = [d.get("reservoir_before", 0.0) for d in seq[0:3]]
        ok = all(rb_i == 0.0 for rb_i in rb) and (ra[2] > ra[1] > ra[0] > 0.0)
        results.append(TestResult("4", ok, f"ra={ra} rb={rb}", {}))

    # Test 5: burst threshold diagnostic
    if logC:
        bh, _ = get_last_bhfdbk(logC, bh_id=1)
        burst = int(bh.get("burst_diag", 0))
        e_dep = bh.get("E_deposited", 0.0)
        ok = (burst == 1 and abs(e_dep) < 1e-20)
        results.append(TestResult("5", ok, f"burst={burst} Edep={e_dep:.8e}", {}))

    # Test 6: feedback kernel geometry
    if logA:
        bh, _ = get_last_bhfdbk(logA, bh_id=1)
        cells = int(bh.get("feedback_kernel_cells", -1))
        temp = bh.get("T_before_mean", 0.0)
        ok = (cells == 1 and abs(temp - 2.0e4) / 2.0e4 < 1e-3)
        results.append(TestResult("6", ok, f"cells={cells} T={temp:.6e}", {}))

    # Test 7: newly-seeded BH skip
    if logA:
        lines = read_lines(logA, "[BHFDBK]")
        first = parse_kv_line(lines[0]) if lines else {}
        bh_id = int(first.get("bh_id", 0))
        ok = (bh_id == -99999 and int(first.get("newly_seeded_skip", 0)) == 1)
        results.append(TestResult("7", ok, f"bh_id={bh_id}", {}))

    # Test 8: MPI determinism
    if logA and logA4:
        a_line = [l for l in read_lines(logA, "[BHFDBK]") if "bh_id=1" in l][-1]
        b_line = [l for l in read_lines(logA4, "[BHFDBK]") if "bh_id=1" in l][-1]
        def strip_wall(line):
            parts = [p for p in line.split() if not p.startswith("feedback_wall_ms=")]
            return " ".join(parts)
        ok = strip_wall(a_line) == strip_wall(b_line)
        results.append(TestResult("8", ok, "diff=none" if ok else "diff=found", {}))

    # Test 9: full lifecycle regression tags present
    if logA:
        text = Path(logA).read_text()
        tags = ["[BHSEED]", "[BHREPOS]", "[BHACCR]", "[BHFDBK]"]
        ok = all(t in text for t in tags)
        results.append(TestResult("9", ok, "tags present" if ok else "missing tags", {}))

    # Test 10: multi-BH diagnostics
    if logD:
        lines = read_lines(logD, "[BHFDBK]")
        parsed = [parse_kv_line(l) for l in lines]
        by_id = {}
        for d in parsed:
            bh_id = int(d.get("bh_id", -99999))
            if bh_id <= 0:
                continue
            by_id.setdefault(bh_id, []).append(d)
        ids = sorted(by_id.keys())
        ok = len(ids) >= 2
        details = f"ids={ids}"
        if ok:
            id1, id2 = ids[0], ids[1]
            seq1 = by_id[id1]
            seq2 = by_id[id2]
            ra1 = [d.get("reservoir_after_accum", d.get("reservoir_after", 0.0)) for d in seq1[:3]]
            ra2 = [d.get("reservoir_after_accum", d.get("reservoir_after", 0.0)) for d in seq2[:3]]
            mode1 = seq1[-1].get("feedback_mode", "")
            mode2 = seq2[-1].get("feedback_mode", "")
            ok = (ra1[2] > ra1[1] > ra1[0] > 0.0 and ra2[2] > ra2[1] > ra2[0] > 0.0 and mode1 != "" and mode2 != "")
            if ok:
                ok = (mode1 != mode2) or (abs(seq1[-1].get("f_Edd", 0.0) - seq2[-1].get("f_Edd", 0.0)) > 1e-8)
            details = f"ids={ids} mode1={mode1} mode2={mode2}"
        results.append(TestResult("10", ok, details, {}))

    # Test 11: performance
    if logA:
        lines = read_lines(logA, "[BHFDBK]")
        vals = []
        for line in lines:
            d = parse_kv_line(line)
            if int(d.get("bh_id", -1)) != 1:
                continue
            v = d.get("feedback_wall_ms", 0.0)
            if v > 0:
                vals.append(v)
        avg = sum(vals) / max(len(vals), 1)
        ok = len(vals) >= 3
        results.append(TestResult("11", ok, f"avg_ms={avg:.4f} n={len(vals)}", {}))

    # Emit summary and JSON
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summary_lines = []
    payload = {"results": []}
    for r in results:
        summary_lines.append(f"{r.name}: {'PASS' if r.passed else 'FAIL'} :: {r.details}")
        payload["results"].append({
            "name": r.name,
            "passed": r.passed,
            "details": r.details,
            "metrics": r.metrics,
        })

    (OUT_ROOT / "feedback_phaseA_matrix_summary.txt").write_text("\n".join(summary_lines) + "\n")
    (OUT_ROOT / "feedback_phaseA_matrix_results.json").write_text(json.dumps(payload, indent=2))

    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
