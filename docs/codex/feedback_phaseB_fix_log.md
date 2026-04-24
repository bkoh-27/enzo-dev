# Feedback Phase B.1 Pre-Merge Fix Log

## Summary

- Branch: `feature/bh-feedback-phaseB`
- Starting commit: `50dc292f9f320e3787cac3ef47c0f96f88bae533`
- Final commit: `HEAD` after the Phase B.1 fix commit; exact hash reported by Codex after commit because a committed file cannot self-reference its own hash.
- Worktree used: `/home/bkoh/enzo-dev`
- Requested report path `/mnt/home/boh10/bh_proj/docs/codex/feedback_phaseB_fix_log.md` was not mounted in this Codex session (`/mnt/home` absent). This report was written at the matching repo-relative path available here: `docs/codex/feedback_phaseB_fix_log.md`.

## Superpowers Use

`superpowers` was available under `/home/bkoh/codex-tools/superpowers`.

Used:
- `using-superpowers`: checked applicable workflow rules.
- `receiving-code-review`: verified Claude's feedback against the code before editing.
- `test-driven-development`: added failing checks before production changes. RED run: `/gpfs/bkoh/bh_dev/feedback_phaseB_fix_review/red3/phaseB_matrix_pre_fix.log`.

The `brainstorming` skill was reviewed. Its full design-approval gate was not used because the user provided a narrow approved pre-merge scope and asked for direct implementation.

## Files Changed

- `src/enzo/Grid_BHFeedbackHandler.C`
- `run/BHSeed/feedback_phaseB_tests/run_feedback_phaseB_matrix.py`
- `run/BHSeed/feedback_phaseA_tests/run_feedback_phaseA_matrix.py`
- `docs/codex/feedback_phaseB_fix_log.md`

Unrelated pre-existing worktree changes were not modified or committed: `src/enzo/uuid/gen_uuid.c`, `docs/review_feedback_phaseB_before_merge.md`, `run/BHSeed/feedback_phaseB_tests/matrix_20260424_234554/`, `run/smoke_20260424_160238/`, `src/enzo/Make.mach.grammar`.

## Reservoir Fix

Chosen fix: store `PARTICLE_ATTRIBUTE_BHFDBK_ENERGY_RESERVOIR` in code energy units, not CGS ergs.

Implementation:
- On read: `stored_code * energy_units -> reservoir_before` in CGS ergs.
- During physics: `E_requested`, thresholds, deposition, `reservoir_after_accum`, and `reservoir_final` remain CGS ergs.
- On write: `reservoir_final / energy_units -> ParticleAttribute`.
- Logs remain human-readable CGS ergs.

Backward compatibility:
- No migration support was added for old unmerged Phase B restart files that stored CGS ergs in the reservoir attribute.
- This is acceptable for this pre-merge patch because Phase B is not merged/production yet. Restart compatibility is preserved for restarts written after this patch.

## Log Schema Fix

Removed duplicate `reservoir_after` from `[BHFDBK]` logs and retained `reservoir_after_accum`.

Before excerpt:

```text
[BHFDBK] ... E_requested=7.74200357e+41 reservoir_before=0.00000000e+00 reservoir_after=7.74200357e+41 burst_diag=0 ... reservoir_after_accum=7.74200357e+41 reservoir_final=7.74200357e+41 ...
```

After excerpt:

```text
[BHFDBK] ... E_requested=7.74200357e+41 reservoir_before=0.00000000e+00 burst_diag=0 ... reservoir_after_accum=7.74200357e+41 reservoir_final=7.74200357e+41 ...
```

## Interactive Allocation Proof

All builds, simulations, Python analysis, and smoke tests were run inside the interactive Slurm allocation:

```text
hostname=grammar-debug
SLURM_JOB_ID=194737
SLURM_JOB_NODELIST=grammar-debug
PWD=/home/bkoh/enzo-dev
```

No builds, simulations, Python analysis, or smoke tests were run on the login node.

## Build And Test Commands

Build:

```bash
cd /home/bkoh/enzo-dev/src/enzo
make -j4 > /gpfs/bkoh/bh_dev/feedback_phaseB_fix_review/build_after_fix.log 2>&1
```

RED test before production fix:

```bash
cd /home/bkoh/enzo-dev
BH_FEEDBACK_PHASEB_OUT_ROOT=/gpfs/bkoh/bh_dev/feedback_phaseB_fix_review/red3/matrix_pre_fix \
  ./run/BHSeed/feedback_phaseB_tests/run_feedback_phaseB_matrix.py \
  > /gpfs/bkoh/bh_dev/feedback_phaseB_fix_review/red3/phaseB_matrix_pre_fix.log 2>&1
```

Final Phase B validation:

```bash
cd /home/bkoh/enzo-dev
BH_FEEDBACK_PHASEB_OUT_ROOT=/gpfs/bkoh/bh_dev/feedback_phaseB_fix_review/after_final/matrix_after_fix \
  ./run/BHSeed/feedback_phaseB_tests/run_feedback_phaseB_matrix.py \
  > /gpfs/bkoh/bh_dev/feedback_phaseB_fix_review/after_final/phaseB_matrix_after_fix.log 2>&1
```

## Pass/Fail Table

| Check | Result | Evidence |
|---|---:|---|
| Build | PASS | `BUILD_AFTER_FIX_RC=0`; build log says `Success!` |
| RED storage check | FAIL as expected | Pre-fix `stored_attr=7.742003571733635e+41` for CGS reservoir |
| RED log schema check | FAIL as expected | Pre-fix log had both `reservoir_after` and `reservoir_after_accum` |
| Phase B matrix final | PASS | `SUMMARY passed=25 failed=0` |
| Reservoir accumulation | PASS | Test 4, `reservoir_final=7.74200357e+41` |
| Forced burst | PASS | Tests 1/2/3/6 |
| Zero-gas behavior | PASS | Test 15 |
| Log schema | PASS | Test 23, `rows_checked=4` |
| Overflow safety | PASS | Test 22, CGS `7.74200357e+41`, stored attr `1.6502882720007718e-25` |
| MPI smoke np=1 vs np=4 | PASS | Test 19 scalar diagnostics match |
| Restart smoke | PASS | Test 7, first-stage `reservoir_final=7.74200357e+41`, restart `reservoir_before=7.74200357e+41` |

## Remaining Risks

- Old restart files from the unmerged Phase B branch that stored the reservoir in CGS ergs will not be interpreted correctly after this patch. No migration was added because Phase B has not merged or shipped.
- The HDF5 storage check uses `h5py`; the default Enzo Python lacks it, so the matrix falls back to `/home/bkoh/miniconda3/envs/yt-env/bin/python`.
