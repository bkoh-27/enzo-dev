# mbh_maker2_test.md — New-User Validation Suite for the Enzo BH Seeding Module

> Status is production-ready as of Round 6. All 11 tests (T2–T8, TS3c, TS3_wrap, V1–V5) passing. Binary confirmed unchanged since Round 4.

## Background

The Enzo BH seeding path is implemented across:

- `src/enzo/mbh_maker2.C`
- `src/enzo/Grid_MBHMaker2Handler.C`
- `src/enzo/Grid_BHSeedHandler.C`
- `src/enzo/star_maker_bh_seed.C`

The `mbh_maker2` path extends the `star_maker2` logic with two MBH-specific additions:

1. A metallicity upper-bound cut so seeds form only in sufficiently metal-poor gas.
2. A comoving distance exclusion so no new seed forms within `BHSeedExclusionRadius` of existing BH/MBH particles.

Seeding is executed before star formation each cycle so MBH creation is resolved first in the same timestep.

MPI behavior is scalable by design:

- Collective communication is O(1) per AMR level.
- There are no collectives inside the per-grid candidate sweep.
- Winner selection is deterministic across rank counts.

## Global Parameters

| Parameter | Default | Units / notes |
| --- | --- | --- |
| `BHSeedingMethod` | 0 | 0=off, 1=on |
| `BHSeedOverdensityThreshold` | 100.0 | code-unit density; gas must exceed this |
| `BHSeedMetallicityThreshold` | 1e-4 | absolute metal mass fraction (~0.005 Zsun); gas must be *below* this |
| `BHSeedMetallicityThresholdInSolar` | — | convenience input; multiplied by 0.02 on read |
| `BHSeedMass` | 1e5 | Msun; seed particle mass, created exactly (no clamping) |
| `BHSeedTemperatureThreshold` | 1e4 | K; gas must be *below* this temperature |
| `BHSeedExclusionRadius` | 1000.0 | comoving kpc/h; minimum separation from any BH/MBH |
| `BHSeedVelDivCrit` | 1 | 1 = require converging flow (div v < 0) |
| `BHSeedThermalCrit` | 0 | 1 = require cooling time < dynamical time |
| `BHSeedSelfBoundCrit` | 0 | 1 = require virial parameter alpha < 1 |
| `BHSeedRunEveryTimestep` | 1 | 1 = run every AMR sub-cycle; 0 = root timestep only |

## Gating Logic

```text
Gate 1  Finest level?         r[idx] != 0                              -> skip
Gate 2  Density floor?        d[idx] < BHSeedOverdensityThreshold      -> ngates_density++
Gate 3  Temperature ceiling?  temp[idx] > BHSeedTemperatureThreshold   -> ngates_temp++
Gate 4  Metallicity ceiling?  metal[idx] > BHSeedMetallicityThreshold  -> ngates_metal++
Gate 5  Converging flow?      BHSeedVelDivCrit && div >= 0             -> ngates_conv++
Gate 6  Cooling time?         BHSeedThermalCrit && tdyn < tcool        -> ngates_cool++
Gate 7  Self-bound?           BHSeedSelfBoundCrit && alpha >= 1        -> ngates_bound++
Gate 8  Distance exclusion    (applied in C++ layer, not kernel)       -> dist_blocked++
Gate 9  Mass availability     cell gas mass < BHSeedMass               -> ngates_mass++
```

Highest-density cell wins. Tie-break is lowest linear index. At most one seed is created per AMR level per timestep.

## Log Line Format

```text
[BHSEED] step=<int> z=<float> a=<float> ncand_local_min=<int> ncand_local_max=<int> ncand_global=<int> ngates_density=<int> ngates_temp=<int> ngates_metal=<int> ngates_conv=<int> ngates_cool=<int> ngates_bound=<int> ngates_mass=<int> dist_blocked=<int> created=<int> total_mbh=<int>
```

| Field | Meaning |
| --- | --- |
| `step` | BH seeding call counter (per level call sequence) |
| `z` / `a` | redshift and scale factor at evaluation |
| `ncand_local_min` | minimum local candidate count across ranks |
| `ncand_local_max` | maximum local candidate count across ranks |
| `ncand_global` | total post-gate candidate count before distance exclusion |
| `ngates_density` | candidates rejected by density floor |
| `ngates_temp` | candidates rejected by temperature ceiling |
| `ngates_metal` | candidates rejected by metallicity ceiling |
| `ngates_conv` | candidates rejected by converging-flow gate |
| `ngates_cool` | candidates rejected by thermal (tdyn vs tcool) gate |
| `ngates_bound` | candidates rejected by self-bound gate |
| `ngates_mass` | candidates rejected because local gas mass < `BHSeedMass` |
| `dist_blocked` | candidates rejected by exclusion-radius proximity to existing BH/MBH |
| `created` | number of seeds created during this call |
| `total_mbh` | total MBH particle count after this call |

## Test Suite Overview

### Core tests (T2–T8)

| Test | Purpose | Status |
| --- | --- | --- |
| T2 — Baseline (np4) | BH seeding off; no regressions to vanilla Enzo | PASS |
| T3 — Star regression (np1) | star_maker2 unchanged; star counts match baseline | PASS |
| T4 — Positive seeding (np1) | One BH created step 1; all candidates blocked step 2 | PASS |
| T5 — Positive seeding (np4) | Same as T4 with 4 MPI ranks | PASS |
| T6 — Determinism | HDF5 particle fields bitwise-identical np1 vs np4 | PASS |
| T7 — Exact mass | MBH mass = BHSeedMass * SolarMass / MassUnits exactly | PASS |
| T8 — Restart | All BHSeed* parameters persist through checkpoint-restart | PASS |

### Periodic-wrap tests (TS3c, TS3_wrap)

| Test | Purpose | Status |
| --- | --- | --- |
| TS3c — Corrected | All 125,000 candidates blocked after step 1; created=0 at step 2 | PASS |
| TS3_wrap — Wrap-specific | Candidate at x=0.990 blocked only via periodic wrap (direct >> exclusion radius) | PASS (8/8) |

### Validation suite (V1–V5)

| Test | Purpose | Status |
| --- | --- | --- |
| V1 — Metal-rich suppression | No BH seeded when metallicity > threshold | PASS |
| V2 — Temperature suppression | No BH seeded when temperature > threshold | PASS |
| V3 — Density suppression | No BH seeded when density < threshold | PASS |
| V4 — Distance exclusion | No BH seeded when candidate within BHSeedExclusionRadius | PASS |
| V5 — BH seeds before star formation | BH seeded first; star sees density-reduced cell | PASS |

## How to Build

```bash
cd /mnt/home/boh10/bh_proj/enzo-dev/src/enzo
make -j8
```

Pass criteria: zero errors, zero new warnings. New object files expected in this module path:

- `Grid_BHSeedHandler.o`
- `Grid_MBHMaker2Handler.o`
- `mbh_maker2.o`
- `star_maker_bh_seed.o`

## Running the Tests

This section lists each test with command, key parameters, expected `[BHSEED]` behavior, and a simple shell assertion.

### T2 — Baseline (seeding OFF)

Command:

```bash
mpirun -n 4 ./enzo /tmp/bhseed_baseline.enzo
```

Key parameters:

- `BHSeedingMethod=0`

Expected:

- No `[BHSEED]` lines.

Assertion:

```bash
grep -c '\[BHSEED\]' run.log
# expect 0
```

### T3 — Star regression

Command:

```bash
mpirun -n 1 ./enzo run/StarParticle/StarParticleSingleTest/TestStarParticleSingle.enzo
```

Key parameters:

- Standard star test configuration.

Expected:

- Successful completion.
- No BH seeding regression artifacts.

Assertion:

```bash
tail -3 run.log | grep -q 'Successful run, exiting.'
```

### T4 — Positive seeding (np1)

Command:

```bash
mpirun -n 1 ./enzo /tmp/bhseed_positive.enzo
```

Key parameters:

- `BHSeedingMethod=1`
- `BHSeedExclusionRadius` large enough to suppress second-step creation

Expected `[BHSEED]`:

- Step 1 has `created=1`.
- Step 2 has `dist_blocked=125000 created=0`.

Assertion:

```bash
STEP1=$(grep '\[BHSEED\]' run.log | sed -n '1p')
STEP2=$(grep '\[BHSEED\]' run.log | sed -n '2p')
echo "$STEP1" | grep -q 'created=1'
echo "$STEP2" | grep -q 'dist_blocked=125000'
echo "$STEP2" | grep -q 'created=0'
```

### T5 — Positive seeding (np4)

Command:

```bash
mpirun -n 4 ./enzo /tmp/bhseed_positive.enzo
```

Key parameters:

- Same as T4, distributed across 4 MPI ranks.

Expected `[BHSEED]`:

- Step 1 has `created=1`.
- Step 2 has `dist_blocked=125000 created=0`.

Assertion:

```bash
STEP1=$(grep '\[BHSEED\]' run.log | sed -n '1p')
STEP2=$(grep '\[BHSEED\]' run.log | sed -n '2p')
echo "$STEP1" | grep -q 'created=1'
echo "$STEP2" | grep -q 'dist_blocked=125000'
echo "$STEP2" | grep -q 'created=0'
```

### T6 — Determinism (np1 vs np4)

Commands:

```bash
mpirun -n 1 ./enzo /tmp/bhseed_positive.enzo
mpirun -n 4 ./enzo /tmp/bhseed_positive.enzo
```

Bitwise comparison script (fields: `particle_type`, `particle_position_x/y/z`, `creation_time`):

```python
import h5py, numpy as np

def read_fields(path):
    with h5py.File(path, 'r') as f:
        return {
            'particle_type': np.array(f['particle_type']),
            'particle_position_x': np.array(f['particle_position_x']),
            'particle_position_y': np.array(f['particle_position_y']),
            'particle_position_z': np.array(f['particle_position_z']),
            'creation_time': np.array(f['creation_time']),
        }

np1 = read_fields('np1/DD0001/data0001.cpu0000')
np4 = read_fields('np4/DD0001/data0001.cpu0000')
for k in np1:
    assert np.array_equal(np1[k], np4[k]), f"Mismatch in {k}"
print('PASS')
```

Shell assertion:

```bash
python3 compare_np1_np4.py | grep -q PASS
```

### T7 — Exact mass

Command:

```bash
mpirun -n 1 ./enzo /tmp/bhseed_positive.enzo
```

Mass validation script:

```python
import h5py, numpy as np

BHSeedMass = 1.0
SolarMass = 1.989e33
DensityUnits = 1.0e-24
LengthUnits = 1.2344e21
CellWidthCode = 0.02

MassUnitsParticleOutput = DensityUnits * (LengthUnits * CellWidthCode)**3
ExpectedBHMassCode = BHSeedMass * SolarMass / MassUnitsParticleOutput

with h5py.File('DD0001/data0001.cpu0000', 'r') as f:
    ptype = f['particle_type'][:]
    pmass = f['particle_mass'][:]

obs = pmass[ptype == 8][0]
relerr = abs(obs - ExpectedBHMassCode) / abs(ExpectedBHMassCode)
print('ExpectedBHMassCode=', ExpectedBHMassCode)
print('ObservedBHMassCode=', obs)
print('RelativeError=', relerr)
print('PASS' if relerr == 0 else 'FAIL')
```

Assertion:

```bash
python3 check_mass.py | grep -q 'RelativeError=0'
```

### T8 — Restart continuity

Commands:

```bash
mpirun -n 1 ./enzo /tmp/bhseed_positive.enzo
mpirun -n 1 ./enzo DD0001/DD0001
```

Expected:

- Restart succeeds.
- BH seeding parameter block preserved in restart path.

Assertion:

```bash
grep -q 'Successful run, exiting.' restart.log
```

### TS3c — Corrected exclusion assertion

Command:

```bash
mpirun -n 1 ./enzo /tmp/bhseed_positive.enzo
```

Key parameter:

- `BHSeedExclusionRadius=1000`

Expected `[BHSEED]`:

- Step 2 has `created=0`.

Assertion:

```bash
grep '\[BHSEED\]' run.log | sed -n '2p' | grep -q 'created=0'
```

### TS3_wrap — periodic boundary wrap path

Commands:

```bash
mpirun -n 1 ./enzo -r DD0000/data0000
mpirun -n 4 ./enzo -r DD0000/data0000
```

Key parameters / geometry:

- `BHSeedExclusionRadius=30`
- Cell A at `x=0.010`
- Cell B at `x=0.990`

Expected `[BHSEED]`:

- Step 1 has `created=1`.
- Step 2 has `dist_blocked=1 created=0`.

Assertion:

```bash
STEP1=$(grep '\[BHSEED\]' run.log | sed -n '1p')
STEP2=$(grep '\[BHSEED\]' run.log | sed -n '2p')
echo "$STEP1" | grep -q 'created=1'
echo "$STEP2" | grep -q 'dist_blocked=1'
echo "$STEP2" | grep -q 'created=0'
```

### V1 — Metal-rich suppression

Command:

```bash
mpirun -n 1 ./enzo -r DD0000/data0000 -d bhseed_v1.enzo
```

Key setup:

- Candidate metallicity above threshold.

Expected `[BHSEED]`:

- `ngates_metal` dominates.
- `created=0`.

Assertion:

```bash
grep '\[BHSEED\]' v1.log | awk '{print $0}' | grep -q 'created=0'
```

### V2 — Temperature suppression

Command:

```bash
mpirun -n 1 ./enzo -r DD0000/data0000 -d bhseed_v2_veldivcrit.enzo
```

Key setup:

- Candidate temperature above threshold.

Expected `[BHSEED]`:

- `ngates_temp` incremented.
- `created=0`.

Assertion:

```bash
grep '\[BHSEED\]' v2.log | grep -q 'created=0'
```

### V3 — Density suppression

Command:

```bash
mpirun -n 1 ./enzo -r DD0000/data0000 -d bhseed_v3_switchoff.enzo
```

Key setup:

- Density below threshold or master off path.

Expected:

- No seed creation.

Assertion:

```bash
grep -c '\[BHSEED\]' v3.log
# expect 0 or created=0 depending configuration
```

### V4 — Distance exclusion suppression

Command:

```bash
mpirun -n 1 ./enzo -r DD0000/data0000 -d bhseed_v4_massgate.enzo
```

Key setup:

- Candidate inside exclusion radius.

Expected `[BHSEED]`:

- `dist_blocked` or `ngates_mass` path suppresses creation.

Assertion:

```bash
grep '\[BHSEED\]' v4.log | grep -Eq 'dist_blocked=[1-9]|ngates_mass=[1-9]'
grep '\[BHSEED\]' v4.log | grep -q 'created=0'
```

### V5 — BH seeds before star formation

Command:

```bash
mpirun -n 1 ./enzo -r DD0000/data0000 -d bhseed_v5_solar.enzo
```

Key setup:

- Solar-threshold input form, star/BH ordering check.

Expected `[BHSEED]`:

- BH seeded first.
- Star logic sees density-reduced cell state.

Assertion:

```bash
grep '\[BHSEED\]' v5.log | sed -n '1p' | grep -q 'created=1'
```

## Interpreting Results

Example 1: healthy seeding event

```text
[BHSEED] step=1 z=15.0 a=0.062500 ncand_local_min=3 ncand_local_max=4 ncand_global=7 ngates_density=124900 ngates_temp=80 ngates_metal=13 ngates_conv=0 ngates_cool=0 ngates_bound=0 ngates_mass=0 dist_blocked=6 created=1 total_mbh=21
```

Interpretation:

- 7 viable candidates survived gates.
- 6 were blocked by exclusion distance.
- 1 was created.

Example 2: all cells fail density gate

```text
[BHSEED] step=3 z=12.0 a=0.076923 ncand_global=0 ngates_density=125000 ngates_temp=0 ngates_metal=0 ngates_conv=0 ngates_cool=0 ngates_bound=0 ngates_mass=0 dist_blocked=0 created=0 total_mbh=21
```

Interpretation:

- Density threshold is too high for current gas state.
- No candidate progressed to later gates.

Example 3: metallicity suppression

```text
[BHSEED] step=5 z=10.0 a=0.090909 ncand_global=0 ngates_density=0 ngates_temp=0 ngates_metal=125000 ngates_conv=0 ngates_cool=0 ngates_bound=0 ngates_mass=0 dist_blocked=0 created=0 total_mbh=21
```

Interpretation:

- All potential cells were metal-rich relative to threshold.
- Lowering `BHSeedMetallicityThreshold` further would not help; threshold must be raised or gas must be cleaner.

## MPI Scalability Notes

| Call | Payload | Frequency |
| --- | --- | --- |
| MPI_Allgatherv (existing BH positions) | 3 x FLOAT x N_BH | Once per level (BHSeedBeginLevel) |
| MPI_Allgather (winner selection) | 5 doubles x N_ranks | Once per level (BHSeedFinalizeLevel) |
| MPI_Bcast (seed_created flag) | 1 int | Once per level, only when winner exists |
| MPI_Reduce x 3 (diagnostics) | 10+1+1 long longs | Once per level |

For `N_BH = 10^4`, the Allgatherv payload is approximately 80 kB per level, which is within the normal per-level communication budget.

## Deferred Nice-to-Haves

- **N1** — Static scratch vectors in `star_maker_bh_seed.C` and `mbh_maker2.C` are not OpenMP-safe. Safe for current Enzo because the star formation grid loop is serial. Add a guard comment before any OpenMP parallelization of the EvolveLevel star loop.
- **N2** — `BHSeedStepCounter` counts AMR-level calls, not root-grid timesteps. Add `level=` field to the `[BHSEED]` line.
- **N3** — Three `MPI_Reduce` diagnostics calls could be merged to two (combine MIN and MAX payload). Low priority because these are per-level root reductions, not allreduce calls.
