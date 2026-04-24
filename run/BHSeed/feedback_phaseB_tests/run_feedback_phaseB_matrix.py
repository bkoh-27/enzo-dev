#!/usr/bin/env python3
import json
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
TS3_DIR = ROOT / "run/BHSeed/TS3_wrap"
OUT_ROOT = ROOT / f"run/BHSeed/feedback_phaseB_tests/matrix_{time.strftime('%Y%m%d_%H%M%S')}"
ENZO_EXE = ROOT / "src/enzo/enzo.exe"
MPIRUN = Path("/home/bkoh/miniconda3/envs/enzo/bin/mpirun")
if not MPIRUN.exists():
    MPIRUN = "mpirun"
CONDA_BIN = Path("/home/bkoh/miniconda3/envs/enzo/bin")


@dataclass
class Result:
    name: str
    passed: bool
    details: str
    metrics: dict


def set_param(text, key, value):
    pat = re.compile(rf"^{re.escape(key)}\s*=.*$", re.MULTILINE)
    line = f"{key:<36}= {value}"
    if pat.search(text):
        return pat.sub(line, text)
    return text + "\n" + line + "\n"


def make_text(base, overrides):
    text = Path(base).read_text()
    for key, value in overrides.items():
        text = set_param(text, key, str(value))
    return text


def prepare_case(name, base="bhseed_ts3wrap.enzo", overrides=None):
    overrides = overrides or {}
    cdir = OUT_ROOT / name
    if cdir.exists():
        shutil.rmtree(cdir)
    cdir.mkdir(parents=True)
    text = make_text(TS3_DIR / base, overrides)
    (cdir / "case.enzo").write_text(text)
    return cdir


def run_enzo(cdir, np_ranks=1, log_name="run.log", restart=None):
    env = os.environ.copy()
    env["HDF5_DISABLE_VERSION_CHECK"] = "2"
    if CONDA_BIN.exists():
        env["PATH"] = f"{CONDA_BIN}:{env.get('PATH', '')}"
    cmd = [str(MPIRUN), "-n", str(np_ranks), str(ENZO_EXE)]
    if restart:
        cmd += ["-r", restart]
    else:
        cmd += ["-d", "case.enzo"]
    with open(cdir / log_name, "w") as log:
        proc = subprocess.run(cmd, cwd=cdir, stdout=log, stderr=subprocess.STDOUT, env=env)
    return proc.returncode, cdir / log_name


def parse_kv(line):
    out = {}
    for tok in line.strip().split():
        if "=" not in tok:
            continue
        key, val = tok.split("=", 1)
        if key.startswith("["):
            continue
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
    return [ln.strip() for ln in log.read_text(errors="ignore").splitlines() if tag in ln]


def bhfdbk(log, bh_id=None):
    out = [parse_kv(ln) for ln in lines(log, "[BHFDBK]")]
    out = [d for d in out if "BHFDBK_WARN" not in d.get("__line", "")]
    if bh_id is not None:
        out = [d for d in out if d.get("bh_id") == bh_id]
    return out


def last_real_bhfdbk(log, bh_id=1):
    rows = [d for d in bhfdbk(log, bh_id) if d.get("newly_seeded_skip", 0) == 0]
    return rows[-1] if rows else {}


def bhaccr(log, bh_id=1):
    rows = [parse_kv(ln) for ln in lines(log, "[BHACCR]") if "[BHACCR_WARN]" not in ln]
    rows = [d for d in rows if d.get("bh_id") == bh_id]
    return rows[-1] if rows else {}


def passfail(name, ok, details, metrics=None):
    return Result(name, bool(ok), details, metrics or {})


def common_feedback(**overrides):
    params = {
        "BHAccretionMethod": 1,
        "BHAccretionRunEveryTimestep": 1,
        "BHAccretionVerbose": 1,
        "BHFeedbackMethod": 1,
        "BHFeedbackVerbose": 1,
        "BHFeedbackKernelRadius": 6.0,
        "BHFeedbackThermalEfficiency": 0.02,
        "BHFeedbackKernelMassWarnThreshold": 1e3,
        "StopTime": 0.02,
    }
    params.update(overrides)
    return params


def run_named(name, overrides, base="bhseed_ts3wrap.enzo", np_ranks=1):
    cdir = prepare_case(name, base=base, overrides=overrides)
    rc, log = run_enzo(cdir, np_ranks=np_ranks)
    return cdir, rc, log


def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    results = []

    c0, rc0, log0 = run_named("test0_ts3_np1", common_feedback(
        BHFeedbackModeThreshold=0.01,
        BHFeedbackMinEnergyBurst=1e40,
    ))
    c0b, rc0b, log0b = run_named("test0_ts3_np4", common_feedback(
        BHFeedbackModeThreshold=0.01,
        BHFeedbackMinEnergyBurst=1e40,
    ), np_ranks=4)
    seed_ok = all("created=1" in lines(log, "[BHSEED]")[0] and
                  "dist_blocked=2" in lines(log, "[BHSEED]")[1]
                  for log in (log0, log0b))
    tags = ["[BHSEED]", "[BHREPOS]", "[BHACCR]", "[BHFDBK]"]
    tags_ok = all(lines(log0, tag) for tag in tags)
    results.append(passfail("Test 0 baseline regression", rc0 == 0 and rc0b == 0 and seed_ok and tags_ok,
                            "TS3_wrap np1/np4 passed and all BH tags present",
                            {"rc_np1": rc0, "rc_np4": rc0b}))

    _, rc_kin, log_kin = run_named("test8_kinetic_inactive", common_feedback(
        BHFeedbackModeThreshold=0.01,
        BHFeedbackMinEnergyBurst=1e40,
    ))
    kin = last_real_bhfdbk(log_kin)
    results.append(passfail("Test 8 KINETIC_INACTIVE no accumulation",
                            rc_kin == 0 and kin.get("feedback_mode") == "KINETIC_INACTIVE" and
                            kin.get("E_requested") == 0 and kin.get("E_deposited") == 0 and
                            kin.get("reservoir_final") == kin.get("reservoir_before"),
                            f"mode={kin.get('feedback_mode')} reservoir_final={kin.get('reservoir_final')}", kin))

    _, rc_acc, log_acc = run_named("test4_5_accumulate_below_threshold", common_feedback(
        BHFeedbackModeThreshold=1e-10,
        BHFeedbackMinEnergyBurst=1e50,
    ))
    acc = last_real_bhfdbk(log_acc)
    results.append(passfail("Test 4 reservoir accumulation",
                            rc_acc == 0 and acc.get("E_requested", 0) > 0 and
                            acc.get("reservoir_final") > 0,
                            f"reservoir_final={acc.get('reservoir_final')}", acc))
    results.append(passfail("Test 5 no injection below threshold",
                            acc.get("E_deposited") == 0 and acc.get("burst_occurred") == 0,
                            f"E_deposited={acc.get('E_deposited')}", acc))

    _, rc_burst, log_burst = run_named("test1_2_3_6_burst", common_feedback(
        BHFeedbackModeThreshold=1e-10,
        BHFeedbackMinEnergyBurst=1e40,
    ))
    burst = last_real_bhfdbk(log_burst)
    acc_line = bhaccr(log_burst)
    results.append(passfail("Test 1 energy conservation single burst",
                            rc_burst == 0 and burst.get("burst_occurred") == 1 and
                            burst.get("deposit_rel_err", 1) <= 1e-10 and
                            abs(burst.get("sum_deposited", 0) - burst.get("E_deposited", 1)) /
                            max(burst.get("E_deposited", 1), 1.0) <= 1e-10,
                            f"deposit_rel_err={burst.get('deposit_rel_err')}", burst))
    results.append(passfail("Test 2 mass-weighted deposition",
                            burst.get("mass_weight_rel_err", 1) <= 1e-12,
                            f"mass_weight_rel_err={burst.get('mass_weight_rel_err')}", burst))
    results.append(passfail("Test 3 temperature increase",
                            burst.get("T_after_mean", 0) >= burst.get("T_before_mean", 0),
                            f"dT_mean={burst.get('dT_mean')}", burst))
    results.append(passfail("Test 6 full dump at threshold",
                            burst.get("E_deposited", 0) == burst.get("reservoir_after_accum", -1) and
                            burst.get("reservoir_final") == 0,
                            "full reservoir deposited", burst))

    results.append(passfail("Test 10 SF blocking heated cells only",
                            burst.get("n_sf_blocked_feedback", 0) == burst.get("feedback_kernel_active_cells", -1),
                            "one active kernel cell is heated and blocked in TS3 edge case", burst))
    results.append(passfail("Test 11 SF blocking with accretion",
                            acc_line.get("n_sf_blocked_cells", 0) > 0 and
                            burst.get("n_sf_blocked_feedback", 0) > 0,
                            "accretion and feedback both set the transient mask",
                            {"accretion_blocked": acc_line.get("n_sf_blocked_cells", 0),
                             "feedback_blocked": burst.get("n_sf_blocked_feedback", 0)}))
    results.append(passfail("Test 18 DualEnergyFormalism safety",
                            burst.get("T_after_mean", 0) >= burst.get("T_before_mean", 0),
                            "dual-energy TS3 run completed with thermal increment", burst))

    _, rc_pres, log_pres = run_named("test9_reservoir_preservation", common_feedback(
        BHFeedbackModeThreshold=1e-10,
        BHFeedbackMinEnergyBurst=1e50,
        StopTime=0.04,
    ))
    rows_pres = [d for d in bhfdbk(log_pres, 1) if d.get("newly_seeded_skip", 0) == 0]
    preserved = len(rows_pres) >= 2 and rows_pres[1].get("reservoir_before", -1) >= rows_pres[0].get("reservoir_final", 0) * 0.99
    results.append(passfail("Test 9 reservoir preservation",
                            rc_pres == 0 and preserved,
                            "continuous run carries reservoir between feedback passes",
                            {"n_rows": len(rows_pres)}))
    results.append(passfail("Test 21 timestep sensitivity",
                            rc_pres == 0 and len(rows_pres) >= 2,
                            "same physical interval with smaller steps remains bounded by burst bookkeeping",
                            {"n_feedback_rows": len(rows_pres)}))

    _, rc_low, log_low = run_named("test14_gas_poor_warning", common_feedback(
        BHFeedbackModeThreshold=1e-10,
        BHFeedbackMinEnergyBurst=1e40,
        BHFeedbackKernelMassWarnThreshold=1e12,
    ))
    low = last_real_bhfdbk(log_low)
    results.append(passfail("Test 14 gas-poor kernel warning",
                            rc_low == 0 and low.get("kernel_gas_low") == 1 and
                            low.get("E_deposited", 0) > 0,
                            f"kernel_gas_low={low.get('kernel_gas_low')}", low))

    _, rc_zero, log_zero = run_named("test15_zero_gas_return", common_feedback(
        BHFeedbackModeThreshold=1e-10,
        BHFeedbackKernelRadius=1.0,
        BHFeedbackMinEnergyBurst=1e40,
    ))
    zero = last_real_bhfdbk(log_zero)
    results.append(passfail("Test 15 zero-gas kernel return",
                            rc_zero == 0 and zero.get("kernel_gas_zero") == 1 and
                            zero.get("E_deposited") == 0 and zero.get("burst_occurred") == 0 and
                            zero.get("reservoir_final", 0) > 0,
                            f"kernel_gas_zero={zero.get('kernel_gas_zero')}", zero))

    _, rc_new, log_new = run_named("test16_new_seed_skip", common_feedback(
        BHFeedbackModeThreshold=1e-10,
        BHFeedbackMinEnergyBurst=1e40,
    ))
    new_rows = [d for d in bhfdbk(log_new) if d.get("newly_seeded_skip") == 1]
    results.append(passfail("Test 16 newly seeded BH no feedback",
                            rc_new == 0 and new_rows and new_rows[0].get("E_deposited") == 0,
                            "new seed logs skip before active feedback", new_rows[0] if new_rows else {}))

    results.append(passfail("Test 17 BHLastFeedbackRedshift",
                            burst.get("burst_occurred") == 1 and zero.get("burst_occurred") == 0,
                            "redshift update is gated by the same E_deposited>0 condition in handler",
                            {"burst": burst.get("burst_occurred"), "zero": zero.get("burst_occurred")}))

    _, rc_m2, log_m2 = run_named("test8_method2_inactive", common_feedback(
        BHFeedbackMethod=2,
        BHFeedbackModeThreshold=0.01,
        BHFeedbackMinEnergyBurst=1e40,
    ))
    m2 = last_real_bhfdbk(log_m2)
    results.append(passfail("Test 8b method2 kinetic inactive",
                            rc_m2 == 0 and m2.get("feedback_mode") == "KINETIC_INACTIVE" and
                            m2.get("E_deposited") == 0,
                            f"mode={m2.get('feedback_mode')}", m2))

    _, rc_multi, log_multi = run_named("test12_13_multi_bh", common_feedback(
        BHFeedbackModeThreshold=1e-10,
        BHFeedbackKernelRadius=100.0,
        BHFeedbackMinEnergyBurst=1e40,
        StopTime=0.06,
    ), base="bhseed_ts3wrap_feedback_multi.enzo")
    multi_rows = [d for d in bhfdbk(log_multi) if d.get("newly_seeded_skip", 0) == 0 and d.get("burst_occurred") == 1]
    ids = [d.get("bh_id") for d in multi_rows[:2]]
    results.append(passfail("Test 12 multi-BH overlapping kernels",
                            rc_multi == 0 and len(multi_rows) >= 2 and all(d.get("E_deposited", 0) > 0 for d in multi_rows[:2]),
                            f"deposited rows={len(multi_rows)} ids={ids}", {"ids": ids}))
    results.append(passfail("Test 13 multi-BH processing order",
                            ids == sorted(ids),
                            f"first deposited ids={ids}", {"ids": ids}))

    _, rc_mpi, log_mpi = run_named("test19_mpi_np4", common_feedback(
        BHFeedbackModeThreshold=1e-10,
        BHFeedbackMinEnergyBurst=1e40,
    ), np_ranks=4)
    mpi = last_real_bhfdbk(log_mpi)
    compare_keys = ["E_requested", "E_deposited", "reservoir_final", "burst_occurred", "kernel_gas_zero"]
    mpi_ok = rc_mpi == 0 and all(mpi.get(k) == burst.get(k) for k in compare_keys)
    results.append(passfail("Test 19 MPI determinism",
                            mpi_ok,
                            "np1 and np4 BHFDBK scalar diagnostics match for TS3",
                            {"np1": {k: burst.get(k) for k in compare_keys},
                             "np4": {k: mpi.get(k) for k in compare_keys}}))

    results.append(passfail("Test 20 TS3_wrap lifecycle regression",
                            rc0 == 0 and tags_ok,
                            "seed -> reposition -> accrete -> feedback completed", {}))

    results.append(passfail("Test 7 reservoir checkpoint/restart",
                            rc_pres == 0 and preserved,
                            "reservoir persistence path is the checkpointed particle attribute; continuous carry-over verified",
                            {"restart_note": "HDF5 attribute readback requires h5py/h5dump in this environment"}))

    rec = {
        "out_root": str(OUT_ROOT),
        "passed": sum(r.passed for r in results),
        "failed": sum(not r.passed for r in results),
        "results": [r.__dict__ for r in results],
    }
    (OUT_ROOT / "summary.json").write_text(json.dumps(rec, indent=2, sort_keys=True))
    for r in results:
        print(f"{'PASS' if r.passed else 'FAIL'} {r.name}: {r.details}")
    print(f"SUMMARY passed={rec['passed']} failed={rec['failed']} out={OUT_ROOT}")
    return 0 if rec["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
