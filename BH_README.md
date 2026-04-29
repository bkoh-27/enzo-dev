# BH_README.md

## What This Document Is For
This guide is for users who are new to this BH seeding module and want to:
- understand what each BH seeding parameter does,
- run a known-good test from scratch,
- confirm the module is working with simple pass/fail checks.

You do not need any prior context (no internal test names required).

## What Was Added
This branch adds a standalone MBH seeding path to Enzo.

Main code files:
- `src/enzo/star_maker_bh_seed.C` (local cell gating)
- `src/enzo/mbh_maker2.C` (candidate kernel evaluation and local candidate list build)
- `src/enzo/Grid_MBHMaker2Handler.C` (grid wrapper)
- `src/enzo/Grid_BHSeedHandler.C` (global gather/sort/accept walk, multi-seed creation, logging)

The BH seeding pass runs before normal star formation when `BHSeedingMethod = 1`.

## How Seeding Works (Simple Version)
For each level pass:
1. Collect existing BH/MBH positions.
2. If `BHSeedRequireFinestLevel = 1`, reject cells covered by a finer AMR grid.
3. Apply local gates in order:
   - density threshold,
   - local 26-neighbor peak check (`BHSeedRequireLocalPeak`),
   - temperature,
   - metallicity,
   - optional extra checks (converging flow / thermal / self-bound).
4. Run the legacy per-cell mass check:
   - `BHSeedLegacyCellMassGate = 0` (default): shadow diagnostic only,
   - `BHSeedLegacyCellMassGate = 1`: hard reject.
5. Reject candidates too close to existing BHs using a linked-cell spatial hash.
6. Evaluate spherical kernel properties around each surviving candidate:
   - enclosed gas mass,
   - mass-weighted metallicity,
   - kernel density peak,
   - kernel completeness flag.
7. Reject candidates with enclosed mass below `BHSeedMinEnclosedMass`.
8. Gather all surviving candidates across MPI ranks.
9. Sort candidates with deterministic lexicographic ranking (`BHSeedRankingOrder` + deterministic tiebreak).
10. Walk sorted candidates and accept up to `BHSeedMaxPerPass` with:
   - exclusion against pre-existing BHs (`BHSeedExclusionMode`),
   - within-pass de-duplication (`BHSeedMinCandidateSeparation`).
11. For accepted candidates, create BHs in global acceptance order and perform kernel-distributed active-zone gas removal.
12. Mark created-seed cells so star formation skips those exact cells in the same pass.
13. Print `[BHSEED]` and verbose `[BHSEED_SEED]` diagnostics.

## BH Parameters (Beginner Table)
Defaults come from `src/enzo/SetDefaultGlobalValues.C`.

| Parameter | Default | What it means in plain language | Units |
| --- | --- | --- | --- |
| `BHSeedingMethod` | `0` | `int`, **Phase 1 active**. Main on/off switch (`0` = disabled, `1` = enabled). Seeding is opt-in; omitted parameters do not run seed creation. | none |
| `BHSeedOverdensityThreshold` | `1000.0` | `float`, **Phase 1 active**. Gas density must be at least this value to be considered for a seed. | code density units |
| `BHSeedTemperatureThreshold` | `1e4` | `float`, **Phase 1 active**. Gas temperature must be below this value. | K |
| `BHSeedMetallicityThreshold` | `1e-4` | `float`, **Phase 1 active**. Gas metallicity must be below this value. Lower value = more metal-poor requirement. | absolute metal mass fraction |
| `BHSeedMetallicityThresholdInSolar` | `FLOAT_UNDEFINED` | `float`, **Phase 1 active input conversion**. Convenience input in solar units (`Z/Zsun`); code converts using `Zsun = 0.02`. | solar metallicity units |
| `BHSeedVelDivCrit` | `1` | `int (bool-style)`, **Phase 1 active**. If `1`, require converging flow (`div(v) < 0`). | boolean-style int |
| `BHSeedThermalCrit` | `0` | `int (bool-style)`, **Phase 1 active**. If `1`, apply cooling-time criterion (`tcool` vs `tdyn`). | boolean-style int |
| `BHSeedSelfBoundCrit` | `0` | `int (bool-style)`, **Phase 1 active**. If `1`, require self-bound gas (`alpha < 1`). | boolean-style int |
| `BHSeedLegacyCellMassGate` | `0` | `int (bool-style)`, **Phase 2 active**. `0` = shadow-only legacy cell-mass diagnostic, `1` = hard reject when cell gas mass `< BHSeedMass`. | boolean-style int |
| `BHSeedRequireFinestLevel` | `1` | `int (bool-style)`, **Phase 1 active**. If `1`, reject cells covered by finer grids. | boolean-style int |
| `BHSeedRequireLocalPeak` | `1` | `int (bool-style)`, **Phase 1 active**. If `1`, require candidate density to be a 26-neighbor local peak. | boolean-style int |
| `BHSeedPatchRadius` | `3.0` | `float`, **Phase 2 active**. Radius of spherical enclosed-property kernel on the candidate grid. | code length |
| `BHSeedMinEnclosedMass` | `1e6` | `float`, **Phase 2 active**. Minimum enclosed gas mass required by the kernel gate. | Msun |
| `BHSeedMass` | `1e5` | `float`, **Phase 1 active**. Mass of each seeded BH particle. | Msun |
| `BHSeedChannel` | `0` | `int`, **Phase 1 active (metadata)**. Stored on each created BH in seed metadata (`bhseed_channel`). | channel id |
| `BHSeedExclusionMode` | `2` | `int`, **Phase 3 active**. Exclusion radius mode: `0` physical kpc, `1` comoving kpc/h, `2` resolution-scaled (`BHSeedExclusionCells * dx`). | mode id |
| `BHSeedExclusionRadius` | `100.0` | `float`, **Phase 1 active**. Minimum allowed distance from any existing BH/MBH. If candidate is closer, it is blocked. | **physical kpc** |
| `BHSeedExclusionCells` | `16` | `int`, **Phase 3 active**. Used when `BHSeedExclusionMode = 2`: candidate-centric exclusion radius `R = BHSeedExclusionCells * dx_candidate`. | cells |
| `BHSeedMinCandidateSeparation` | `3.0` | `float`, **Phase 3 active**. Minimum within-pass separation between accepted candidates. | physical kpc |
| `BHSeedMaxPerPass` | `10` | `int`, **Phase 3 active**. Maximum number of accepted candidates processed per level pass. | count |
| `BHSeedRunEveryTimestep` | `0` | `int (bool-style)`, **Phase 1 active**. If `1`, run every sub-cycle; if `0`, run with root-grid cadence logic. | boolean-style int |
| `BHSeedRankingOrder` | `0` | `int`, **Phase 3 active**. Ranking key order: `0` enclosed-mass-first, `1` density-peak-first. | enum id |
| `BHSeedVerbose` | `1` | `int`, **Phase 2 active**. Controls BH seeding log verbosity (`[BHSEED]` + `[BHSEED_SEED]` kernel metadata line). | verbosity level |
| `BHSeedDeterministicTiebreak` | `1` | `int (bool-style)`, **Phase 3 active**. Enables deterministic position-based final tie-break in global sort. | boolean-style int |

### Important Unit Note for Exclusion Radius Modes
- Mode 0 (`BHSeedExclusionMode = 0`): `BHSeedExclusionRadius` is interpreted as **physical kpc**.
- Mode 1 (`BHSeedExclusionMode = 1`): `BHSeedExclusionRadius` is interpreted as **comoving kpc/h**.
  - physical conversion used in logs: `R_phys_kpc = R_com_kpc/h * a_phys / h`.
- Mode 2 (`BHSeedExclusionMode = 2`): radius is candidate-centric and computed as
  - `R_excl = BHSeedExclusionCells * dx_candidate`.

The runtime log prints both:
- `excl_phys_kpc`
- `excl_com_kpch`

## Detailed Parameter Explanations (Added in v1.5)
This section expands the quick parameter table with practical behavior notes for each added BH seeding control.

### 1) Enable and cadence controls
- `BHSeedingMethod`
  - Master on/off switch for the MBH seeding module.
  - `0`: seeding code path is skipped.
  - `1`: seeding path runs before star formation on each eligible level pass.
- `BHSeedRunEveryTimestep`
  - Controls cadence under subcycling.
  - `0` (default): follows root-grid cadence logic.
  - `1`: force seeding evaluation every sub-step.

### 2) Core local gas gates
- `BHSeedOverdensityThreshold`
  - Minimum local gas density gate. Cells below this are rejected early (`ngates_density`).
  - Raising this value strongly reduces candidate count.
- `BHSeedTemperatureThreshold`
  - Maximum local gas temperature gate (`ngates_temp`).
  - Lower values prefer colder gas and reduce seeds.
- `BHSeedMetallicityThreshold`
  - Maximum local metallicity fraction gate (`ngates_metal`).
  - Lower values restrict seeding to more metal-poor gas.
- `BHSeedMetallicityThresholdInSolar`
  - Convenience input in solar units (`Z/Zsun`).
  - Converted internally to absolute fraction (using `Zsun = 0.02`) and mapped to `BHSeedMetallicityThreshold`.

### 3) Structural and optional physics gates
- `BHSeedRequireFinestLevel`
  - If `1`, reject cells covered by finer AMR subgrids (`ngates_finestlevel`).
  - Helps avoid parent-level duplicates.
- `BHSeedRequireLocalPeak`
  - If `1`, require candidate to be a 26-neighbor local density peak (`ngates_peak`).
  - If `0`, many more cells can pass downstream gates.
- `BHSeedVelDivCrit`
  - If `1`, require converging flow (`div(v) < 0`, tracked by `ngates_conv`).
- `BHSeedThermalCrit`
  - If `1`, apply cooling-time/dynamical-time filter (`ngates_cool`).
- `BHSeedSelfBoundCrit`
  - If `1`, apply self-bound criterion (`alpha < 1`, tracked by `ngates_bound`).

### 4) Kernel and mass constraints
- `BHSeedPatchRadius`
  - Radius of the spherical kernel used for enclosed-property evaluation.
  - Larger values sample broader environment but can increase truncation near boundaries.
- `BHSeedMinEnclosedMass`
  - Minimum enclosed kernel gas mass required to pass (`ngates_enclosedmass`).
  - Runtime hard check: must be `>= BHSeedMass`.
- `BHSeedMass`
  - Fixed mass assigned to each created BH seed.
  - Also the target mass removed from gas during kernel-distributed removal.
- `BHSeedLegacyCellMassGate`
  - Legacy single-cell mass check behavior:
  - `0` (default): shadow-only (`nlegacy_cellmass_would_fail` increments, no rejection).
  - `1`: hard reject at legacy gate (`ngates_mass`).

### 5) Exclusion geometry controls
- `BHSeedExclusionMode`
  - Chooses exclusion-radius interpretation for candidate-vs-existing-BH blocking.
  - `0`: fixed physical kpc (`BHSeedExclusionRadius` interpreted as physical).
  - `1`: fixed comoving kpc/h (`BHSeedExclusionRadius` interpreted as comoving).
  - `2`: resolution-scaled (`BHSeedExclusionCells * dx_candidate`).
- `BHSeedExclusionRadius`
  - Radius input used by mode 0 or mode 1 (units depend on mode).
- `BHSeedExclusionCells`
  - Integer cell multiplier used only in mode 2.

### 6) Ranking, de-duplication, and multiplicity controls
- `BHSeedRankingOrder`
  - Global lexicographic sort priority among kernel-qualified candidates.
  - `0`: enclosed-mass first, then peak density, then metallicity, then tiebreak.
  - `1`: peak-density first, then enclosed mass, then metallicity, then tiebreak.
- `BHSeedDeterministicTiebreak`
  - Enables deterministic position-based tie resolution in global sorting.
  - Keeps acceptance order reproducible across MPI layouts.
- `BHSeedMinCandidateSeparation`
  - Within-pass de-duplication radius in physical kpc.
  - If two accepted-walk candidates are too close, lower-ranked candidate is rejected (`ncandidates_dedup_rejected`).
- `BHSeedMaxPerPass`
  - Hard cap on accepted candidates processed per level pass.
  - Acceptance walk stops when this limit is reached (`walk_stopped_at_max=1`).

### 7) Metadata and diagnostics controls
- `BHSeedChannel`
  - User channel tag stored in BH metadata (`bhseed_channel`) for downstream analysis.
- `BHSeedVerbose`
  - Logging verbosity:
  - `0`: only core structured pass-level logs.
  - `1`: include per-seed `[BHSEED_SEED]` metadata lines.
  - `>=2`: include debug conservation diagnostics (`[BHSEED_DEBUG]`).

### 8) Practical tuning effects (quick guide)
- If `ncand_global` is very high:
  - raise `BHSeedOverdensityThreshold`, enable/keep `BHSeedRequireLocalPeak=1`, tighten `BHSeedTemperatureThreshold` and `BHSeedMetallicityThreshold`.
- If almost everything is blocked by exclusion (`dist_blocked` high):
  - reduce exclusion radius (or choose a different exclusion mode consistent with your science intent).
- If many candidates fail enclosed mass (`ngates_enclosedmass` high):
  - increase local gas availability or reduce `BHSeedMinEnclosedMass` while keeping `BHSeedMinEnclosedMass >= BHSeedMass`.
- If many accepted candidates are skipped (`nseeds_skipped_insufficient_gas` high):
  - reduce `BHSeedMass`, reduce `BHSeedPatchRadius`, or improve active-zone gas support near boundaries.

## Runtime Validation (Phase 3)
At startup, Enzo aborts if:
- `BHSeedMinEnclosedMass < BHSeedMass`

This protects fixed-mass kernel removal from impossible parameter combinations.

## Recommended Explicit Settings (Do Not Rely on Defaults)
For reproducibility, set BH flags explicitly in your parameter file:

```ini
BHSeedingMethod                   = 1
BHSeedOverdensityThreshold        = 1000
BHSeedTemperatureThreshold        = 1e4
BHSeedMetallicityThreshold        = 1e-4
BHSeedMetallicityThresholdInSolar = 1e-4
BHSeedPatchRadius                 = 3.0
BHSeedMinEnclosedMass             = 1e6
BHSeedMass                        = 1e5
BHSeedChannel                     = 0
BHSeedExclusionMode               = 2
BHSeedExclusionRadius             = 100.0
BHSeedExclusionCells              = 16
BHSeedMinCandidateSeparation      = 3.0
BHSeedMaxPerPass                  = 10
BHSeedRunEveryTimestep            = 0
BHSeedRankingOrder                = 0
BHSeedVelDivCrit                  = 1
BHSeedLegacyCellMassGate          = 0
BHSeedThermalCrit                 = 0
BHSeedSelfBoundCrit               = 0
BHSeedRequireFinestLevel          = 1
BHSeedRequireLocalPeak            = 1
BHSeedVerbose                     = 1
BHSeedDeterministicTiebreak       = 1
```

## Understanding `[BHSEED]` Log Output
You will see lines like:

```text
[BHSEED] step=... level=... z=... a_phys=... excl_phys_kpc=... excl_com_kpch=...
cell_code=... nbins=... ncand_local_min=... ncand_local_max=... ncand_global=...
ngates_density=... ngates_temp=... ngates_metal=... ngates_conv=... ngates_cool=...
ngates_bound=... ngates_mass=... dist_blocked=... created=... total_mbh=... pre_cache_bh=...
ngates_finestlevel=... ngates_peak=...
nkernel_evaluated=... nkernel_truncated=... nlegacy_cellmass_would_fail=...
ngates_enclosedmass=... seeding_wall_ms=...
nseeds_created=... ncandidates_gathered=... ncandidates_dedup_rejected=...
ncandidates_exclusion_rejected=... walk_stopped_at_max=...
nseeds_skipped_insufficient_gas=...
```

Key fields:
- `created`: number of BH seeds created in this level pass.
- `dist_blocked`: candidates rejected because they are too close to an existing BH.
- `ngates_mass`: candidates rejected because cell gas mass is less than `BHSeedMass`.
- `ncand_global`: total candidates after local gates.
- `total_mbh`: BH count in global cache after this pass.
- `pre_cache_bh`: BH count already in cache before this pass.
- `ngates_finestlevel`: candidates rejected because cell is covered by finer AMR grid.
- `ngates_peak`: candidates rejected by 26-neighbor local density peak check.
- `nkernel_evaluated`: candidates that reached kernel evaluation.
- `nkernel_truncated`: kernel evaluations marked incomplete (`BHSeedKernelComplete=0`).
- `nlegacy_cellmass_would_fail`: shadow count for legacy cell-mass check.
- `ngates_enclosedmass`: candidates rejected by `BHSeedMinEnclosedMass`.
- `seeding_wall_ms`: wall-clock time for the full level seeding pass (ms).
- `nseeds_created`: authoritative Phase 3 created-seed count.
- `ncandidates_gathered`: global candidate count before ranking/acceptance.
- `ncandidates_dedup_rejected`: candidates rejected by within-pass separation.
- `ncandidates_exclusion_rejected`: candidates rejected by exclusion during acceptance walk.
- `walk_stopped_at_max`: `1` if walk stopped at `BHSeedMaxPerPass`.
- `nseeds_skipped_insufficient_gas`: accepted candidates skipped at creation time because active-zone removable mass is insufficient.

Log compatibility note:
- New fields are appended at the end of the existing `[BHSEED]` line after all legacy and Phase 1 fields.

When `BHSeedVerbose >= 1`, each created seed also emits:

```text
[BHSEED_SEED] level=... x=... y=... z=... channel=... redshift=... patch_mass=...
patch_metal=... patch_density_peak=... kernel_complete=... host_dm_density=... accept_rank=...
```

This line reports per-seed metadata values written into particle attributes/HDF5.
`patch_density_peak` is the kernel maximum (not the candidate-cell density placeholder used in Phase 1).

## Seed Metadata Fields (Phase 3)
The BH seeding path stores 8 metadata attributes on MBH particles.
All are checkpoint/restart-safe and migrate with particles across MPI ranks.

Note on storage type:
- particle attributes are float-typed arrays, so int-semantic fields are stored as float and cast back to int when read for logging/analysis.

| Metadata field | HDF5 dataset label | Slot index (`no WINDS` / `WINDS`) | Phase 2 value | Phase status |
| --- | --- | --- | --- | --- |
| `BHSeedChannel` | `bhseed_channel` | `4 / 7` | parameter value `BHSeedChannel` | active in Phase 1+ |
| `BHSeedRedshift` | `bhseed_redshift` | `5 / 8` | current redshift (or `0` in non-cosmological runs) | active in Phase 1+ |
| `BHSeedPatchMass` | `bhseed_patch_mass` | `6 / 9` | kernel enclosed gas mass | active in Phase 2 |
| `BHSeedPatchMetallicity` | `bhseed_patch_metallicity` | `7 / 10` | kernel mass-weighted metallicity | active in Phase 2 |
| `BHSeedPatchDensityPeak` | `bhseed_patch_density_peak` | `8 / 11` | maximum density inside evaluated kernel | active in Phase 2 |
| `BHSeedKernelComplete` | `bhseed_kernel_complete` | `9 / 12` | `1` complete, `0` truncated | active in Phase 2 |
| `BHSeedHostDMDensity` | `bhseed_host_dm_density` | `10 / 13` | local DM density if available, else `-1.0` | active in Phase 2 |
| `BHSeedAcceptRank` | `bhseed_accept_rank` | `11 / 14` | global acceptance order (`1, 2, 3, ...`) | active in Phase 3 |

Quick inspection example:

```bash
h5dump -d "bhseed_channel" DD0001/data0001
```

If your output uses per-CPU files/groups, use the full dataset path:

```bash
h5dump -d "/Grid00000001/bhseed_channel" DD0001/data0001.cpu0000
```

## Test You Can Run From Scratch (Step-by-Step)
All files needed are already in:
- `run/BHSeed/TS3_wrap/bhseed_ts3wrap.enzo`
- `run/BHSeed/TS3_wrap/make_ts3wrap_ic.py`
- `run/BHSeed/TS3_wrap/assert_ts3wrap.sh`

### Step 1: Build Enzo
```bash
cd /mnt/home/boh10/bh_proj/enzo-dev/src/enzo
module load gcc openmpi
make -j8
```

What you should see:
- build completes with no errors,
- `enzo.exe` exists in this directory.

### Step 2: Initial Condition Setup (Updated for Phase 1)
You no longer need to run `make_ts3wrap_ic.py` for the pass/fail TS3_wrap check.
The test IC pattern is now planted directly by the `ProblemType = 90` initializer in code (Option C fix), so the test is self-contained when launched with `-d bhseed_ts3wrap.enzo`.

Optional (legacy utility only): you can still run the old IC generator script for manual HDF5 inspection/debugging.

```bash
cd /mnt/home/boh10/bh_proj/enzo-dev/run/BHSeed/TS3_wrap
python3 make_ts3wrap_ic.py
```

What you should see:
- script prints created file info,
- `DD0000/data0000` is generated (legacy artifact; not required by the current TS3_wrap pass/fail runner).

### What the TS3_wrap IC Pattern Contains (Simple Explanation)
The TS3_wrap setup uses a tiny controlled pattern so you can clearly see when BH seeding works.
In the current test flow, this pattern is generated internally by the initializer.
The legacy `DD0000/data0000` script output encodes the same pattern for inspection.

Grid and fields:
- Active grid is `50 x 50 x 50` (125,000 active cells).
- HDF5 arrays are `56 x 56 x 56` because of 3 ghost zones on each side.
- The file has 7 fields:
  - `Density`
  - `x-velocity`
  - `y-velocity`
  - `z-velocity`
  - `TotalEnergy`
  - `GasEnergy`
  - `Metal_Density`

Cell setup:
- Most cells are background:
  - density = `1.0`
  - temperature = `2e4 K`
  - metallicity fraction = `0.002`
- Two special cells are made to pass the BH gates:
  - Cell A at `(ix,iy,iz) = (0,24,24)` with:
    - density = `10.0`
    - temperature = `500 K`
    - metallicity fraction = `1e-5`
  - Cell B at `(ix,iy,iz) = (49,24,24)` with:
    - density = `9.0`
    - temperature = `500 K`
    - metallicity fraction = `1e-5`

Why this is useful:
- In `bhseed_ts3wrap.enzo`, thresholds are:
  - `BHSeedOverdensityThreshold = 5.0`
  - `BHSeedTemperatureThreshold = 1e4`
  - `BHSeedMetallicityThreshold = 1e-4`
- So only Cell A and Cell B are good BH candidates.
- Cell A has slightly higher density, so it seeds first.

Distance geometry (periodic-box test):
- Box size is `100 kpc`, so each active cell is `2 kpc`.
- Cell A is near `x = 0.01`, Cell B is near `x = 0.99`.
- Direct distance is `98 kpc` (far), but wrapped periodic distance is `2 kpc` (near).
- Exclusion radius is `30 kpc`, so Cell B gets blocked by wrap after Cell A seeds.

### Step 3: Run Single-Rank Smoke Test
```bash
cd /mnt/home/boh10/bh_proj/enzo-dev/run/BHSeed/TS3_wrap
mpirun -n 1 ../../../src/enzo/enzo.exe -d bhseed_ts3wrap.enzo > ts3wrap_np1.log 2>&1
grep '\[BHSEED\]' ts3wrap_np1.log
```

Expected result (key fields):
- first `[BHSEED]` line contains `ncand_global=2` and `created=1`,
- second `[BHSEED]` line contains `dist_blocked=2` and `created=0`.

Simple explanation:
- pass 1: one BH is created,
- pass 2: two candidate cells are blocked by the exclusion rule, so no new BH is created.

### Step 4: Run 4-Rank Check (Determinism)
```bash
cd /mnt/home/boh10/bh_proj/enzo-dev/run/BHSeed/TS3_wrap
mpirun -n 4 ../../../src/enzo/enzo.exe -d bhseed_ts3wrap.enzo > ts3wrap_np4.log 2>&1
grep '\[BHSEED\]' ts3wrap_np4.log
```

Expected result:
- same key behavior as single-rank run:
  - line 1 has `created=1`,
  - line 2 has `dist_blocked=2` and `created=0`.
  - appended Phase 1 counters appear at line end (`ngates_finestlevel`, `ngates_peak`),
  - appended Phase 2 counters also appear at line end (`nkernel_evaluated`, `nkernel_truncated`, `nlegacy_cellmass_would_fail`, `ngates_enclosedmass`, `seeding_wall_ms`).

### Step 5: Run Built-In Assertions
```bash
cd /mnt/home/boh10/bh_proj/enzo-dev/run/BHSeed/TS3_wrap
bash assert_ts3wrap.sh
```

Expected result:
- script prints all checks as PASS,
- final line is `TS3_WRAP: PASS`.

## Why This Test Is Useful
This one test verifies all of the following:
- BH seeding actually creates a seed (`created=1`),
- exclusion radius blocks additional seeds (`dist_blocked`),
- periodic wrap logic works,
- result is deterministic across MPI rank counts (1 rank vs 4 ranks).

## Troubleshooting (Quick)
- No `[BHSEED]` lines:
  - check `BHSeedingMethod = 1`,
  - check stdout log capture (`> ts3wrap_np1.log 2>&1`).
- `created=0` on first line:
  - gates are too strict for your setup (density/temp/metal/flow/thermal/self-bound),
  - `BHSeedMass` may be too large (`ngates_mass` will increase).
- Too many seeds:
  - exclusion radius likely too small for your box/redshift,
  - confirm `excl_phys_kpc` and `excl_com_kpch` in log line.

## Full-Scale Production Run Notes
- Always set BH parameters explicitly in your production `.enzo` file.
- Capture stdout and archive it; `[BHSEED]` lines are your primary debugging signal.
- Track `ngates_*`, `dist_blocked`, `created`, and `pre_cache_bh` over time.

## BH Accretion Phase A (Diagnostics Only)
Phase A adds a per-BH diagnostic pass for accretion physics. It is measurement-only.

What Phase A does:
- loops over MBH particles each level update,
- evaluates a spherical gas kernel around each BH,
- classifies kernel gas into hot/cold channels using `t_cool / t_dyn` (with deterministic temperature fallback),
- computes channel-averaged gas properties,
- computes diagnostic hot/cold Bondi-family rates, Eddington rate/ratio, and a diagnostic Eddington cap split,
- writes per-BH diagnostics to `[BHACCR]` log lines,
- updates BH accretion metadata attributes needed by later phases.

What Phase A does **not** do:
- no gas removal,
- no BH mass growth,
- no reservoir drainage,
- no feedback coupling.

Call ordering:
- BH seeding finalize completes first (`BHSeedFinalizeLevel`),
- then `BHAccretionDiagnosticHandler`,
- then normal star/active particle handlers.

### Accretion Parameters (Phase A)
Defaults come from `src/enzo/SetDefaultGlobalValues.C`.

| Parameter | Default | Phase A status |
| --- | --- | --- |
| `BHAccretionMethod` | `0` | active only when explicitly enabled (`0` off, `1` two-channel diagnostics/source path) |
| `BHAccretionKernelRadius` | `3.0` | active (physical kpc) |
| `BHAccretionRemovalRadius` | `1` | reserved (Phase B) |
| `BHAccretionTSplitFloor` | `5e5` | active (K, fallback split) |
| `BHAccretionColdModel` | `0` | active (AM-suppressed Bondi) |
| `BHAccretionCVisc` | `6.283` | active |
| `BHAccretionNHStar` | `0.1` | active (cm^-3) |
| `BHAccretionBeta` | `1.0` | active |
| `BHAccretionAlphaMax` | `10.0` | active |
| `BHAccretionRadiativeEfficiency` | `0.1` | active (`epsilon_r` for Eddington rate) |
| `BHAccretionSuperEddington` | `0` | reserved (Phase C) |
| `BHAccretionSuperEddFactor` | `1.0` | reserved (Phase C) |
| `BHAccretionUseReservoir` | `0` | reserved (v2) |
| `BHAccretionRemovalMode` | `0` | reserved (Phase B) |
| `BHAccretionVerbose` | `1` | active |
| `BHAccretionRunEveryTimestep` | `0` | active (`0` finest-only cadence, `1` every processed update) |

### Runtime Validation and Warnings
Hard startup errors:
- `BHAccretionMethod` must be `0` or `1`.
- `BHAccretionKernelRadius > 0`.
- `BHAccretionCVisc > 0`.
- `0 < BHAccretionRadiativeEfficiency <= 1`.
- `BHAccretionColdModel == 0` (Phase A restriction).

Runtime warning examples (`[BHACCR_WARN]`):
- kernel radius under-resolved (`kernel_radius_over_dx < 1.5`),
- kernel likely exceeds ghost support (`kernel_radius_over_dx > 3`),
- fallback fraction exceeds 50% of kernel cells,
- `alpha_boost` at configured cap (`BHAccretionAlphaMax`).

### `[BHACCR]` Log Line
One line per BH per diagnostic pass. Includes:
- context: `step`, `level`, `z`, `bh_id`, `bh_mass`,
- kernel split: `f_hot`, `f_cold`, `n_hot_cells`, `n_cold_cells`, `n_fallback_cells`,
- averaged gas properties: `rho_hot_avg`, `rho_cold_avg`, `T_hot_avg`, `T_cold_avg`, `cs_hot_avg`,
- kinematics: `V_rot_cold`, `v_rel_hot`, `v_rel_cold`,
- rates: `Mdot_hot_raw`, `Mdot_cold_raw`, `Mdot_total_raw`, `Mdot_Edd`, `f_Edd`,
- channel factors: `alpha_boost`, `f_AM`,
- diagnostic capped split: `Mdot_actual`, `Mdot_hot_actual`, `Mdot_cold_actual`, `cap_active`,
- performance: `accretion_diag_wall_ms`.

Units in logs:
- densities in cgs,
- temperatures in K,
- velocities in cm/s,
- rates in Msun/yr.

### BH Accretion Metadata Attributes
Phase A registers/checkpoints the following MBH attributes:
- `BHAccretedMass` (stays unchanged in Phase A),
- `BHReservoirMass` (stays unchanged in Phase A),
- `BHLastAccretionRedshift` (unchanged in Phase A),
- `BHLastEddingtonRatio` (updated every diagnostic pass),
- `BHFormationMass` (seed-set value, preserved).

These are mapped in `typedefs.h` and persisted via normal particle attribute grid I/O.

### Minimal Phase A Enable Block
```ini
BHAccretionMethod                = 1
BHAccretionKernelRadius          = 3.0
BHAccretionTSplitFloor           = 5.0e5
BHAccretionColdModel             = 0
BHAccretionCVisc                 = 6.283
BHAccretionNHStar                = 0.1
BHAccretionBeta                  = 1.0
BHAccretionAlphaMax              = 10.0
BHAccretionRadiativeEfficiency   = 0.1
BHAccretionVerbose               = 1
BHAccretionRunEveryTimestep      = 0
```

## BH Accretion Phase B (Active Gas Removal + BH Growth)
Phase B keeps the full Phase A diagnostic pass and then activates physical accretion.

What Phase B adds:
- converts capped diagnostic rates into requested removal mass (`dm_requested`),
- evaluates a removal kernel (`BHAccretionRemovalRadius` in cell-width units),
- removes gas from active-zone cells only (deterministic spherical k-j-i loop),
- applies proportional multi-cell removal with residual balancing,
- updates BH particle mass by the actually removed gas (`dm_removed`),
- updates cumulative BH accretion metadata (`BHAccretedMass`, `BHLastAccretionRedshift`),
- marks every cell with `dm_i > 0` in the transient SF mask,
- computes/logs momentum-omission diagnostics and optional warnings,
- skips same-pass accretion for newly seeded BHs (particles with temporary `INT_UNDEFINED` id).

What Phase B still does **not** do:
- BH velocity update from accretion momentum (diagnostics only),
- reservoir drainage,
- super-Eddington accretion mode,
- torque-limited cold channel,
- feedback coupling.

### Accretion Parameters (Phase B Active Set)
| Parameter | Default | Phase B status |
| --- | --- | --- |
| `BHAccretionMethod` | `0` | active only when explicitly enabled (`0` off, `1` two-channel) |
| `BHAccretionKernelRadius` | `3.0` | active (diagnostic kernel, physical kpc) |
| `BHAccretionRemovalRadius` | `1` | active (removal kernel radius in cell widths) |
| `BHAccretionRemovalMode` | `0` | active (`0` multi-cell, `1` single-cell debug) |
| `BHAccretionTSplitFloor` | `5e5` | active |
| `BHAccretionColdModel` | `0` | active (AM-suppressed Bondi) |
| `BHAccretionCVisc` | `6.283` | active |
| `BHAccretionNHStar` | `0.1` | active |
| `BHAccretionBeta` | `1.0` | active |
| `BHAccretionAlphaMax` | `10.0` | active |
| `BHAccretionRadiativeEfficiency` | `0.1` | active |
| `BHAccretionIgnoredDVWarn` | `1.0` | active warning threshold (km/s) |
| `BHAccretionIgnoredPFracWarn` | `0.01` | active warning threshold (dimensionless) |
| `BHAccretionVerbose` | `1` | active |
| `BHAccretionRunEveryTimestep` | `0` | active cadence control |
| `BHAccretionSuperEddington` | `0` | reserved (Phase C) |
| `BHAccretionSuperEddFactor` | `1.0` | reserved (Phase C) |
| `BHAccretionUseReservoir` | `0` | reserved (v2) |

### Phase B Bookkeeping Tiers
- Tier 1 (raw): `Mdot_hot_raw`, `Mdot_cold_raw`, `Mdot_total_raw`.
- Tier 2 (capped): Eddington downscale via `frac_cap`, producing `dm_requested`.
- Tier 3 (gas-limited): removal-kernel availability via `frac_gas`, producing `dm_removed`.

Key realized invariants:
- `Mdot_hot_realized + Mdot_cold_realized = dm_removed / dt`.
- `BH_mass_new = BH_mass_old + dm_removed`.
- `BH_mass_current = BHFormationMass + BHAccretedMass` (metadata check).

### New Phase B `[BHACCR]` Fields
Phase B appends:
- `dm_requested`, `dm_removed`, `dm_removed_msun`,
- `frac_cap`, `frac_gas`, `removal_gas_limited`,
- `Mdot_hot_realized`, `Mdot_cold_realized`,
- `bh_mass_new`,
- `removal_cells`, `n_sf_blocked_cells`,
- `acc_ignored_dv_kms`, `acc_ignored_p_frac`, `acc_momentum_warn`,
- `accretion_wall_ms`.

### Phase C `[BHACCR]` Fields
Phase C appends:
- `lambda_edd`: raw `Mdot_total_raw / Mdot_Edd`, using the physical Eddington
  rate before `BHFeedbackEddingtonFactor` is applied,
- `edd_factor`: runtime value of `BHFeedbackEddingtonFactor`.

These are appended to the end of the line. Existing Phase A/B fields are not
renamed or reordered.

### Phase B Warning Lines (`[BHACCR_WARN]`)
- gas-limited removal (`removal_gas_limited=1`),
- momentum-omission warning when
  - `acc_ignored_dv_kms > BHAccretionIgnoredDVWarn` or
  - `acc_ignored_p_frac > BHAccretionIgnoredPFracWarn`.

### Minimal Phase B Enable Block
```ini
BHAccretionMethod                = 1
BHAccretionKernelRadius          = 3.0
BHAccretionRemovalRadius         = 1
BHAccretionRemovalMode           = 0
BHAccretionTSplitFloor           = 5.0e5
BHAccretionColdModel             = 0
BHAccretionCVisc                 = 6.283
BHAccretionNHStar                = 0.1
BHAccretionBeta                  = 1.0
BHAccretionAlphaMax              = 10.0
BHAccretionRadiativeEfficiency   = 0.1
BHAccretionIgnoredDVWarn         = 1.0
BHAccretionIgnoredPFracWarn      = 0.01
BHAccretionVerbose               = 1
BHAccretionRunEveryTimestep      = 0
```

## BH Reposition Phase B (Active Movement + Diagnostics)
Phase B keeps the full Phase A diagnostic pass and then applies BH movement
toward the active-zone density peak before accretion:

`BHSeedFinalizeLevel -> BHRepositionDiagnosticHandler -> BHAccretionDiagnosticHandler`

Call-site behavior in `EvolveLevel`:
- reposition handler runs only when `BHRepositionMethod > 0` **or**
  `BHRepositionVerbose > 0`.
- `BHRepositionMethod=0` and `BHRepositionVerbose=0` means fully off
  (no handler call, no `[BHREPOS]` output).
- `BHRepositionMethod=0` with `BHRepositionVerbose>0` is diagnostic-only:
  it logs `[BHREPOS]` but never moves a particle.

What Phase B does:
- keeps deterministic BH selection/order/ownership and newly-seeded skip logic from Phase A,
- computes diagnostic and active-zone peaks exactly as in Phase A,
- moves BHs only toward the **active-zone** peak:
  - method `1`: rate-limited drift, cap = `BHRepositionMaxDisplacement * dx_local`,
  - method `2`: teleport debug mode to active-zone peak,
- applies a defensive active-zone clamp for floating-point edge cases,
- rejects active movement when the requested search sphere exceeds local
  ghost-zone/grid support; diagnostics are still logged with
  `active_reposition_rejected=1`,
- writes updated `ParticlePosition[dim][p]` directly,
- leaves `ParticleVelocity` unchanged.

What Phase B still does **not** do:
- modify baryon fields,
- modify potential fields,
- apply any velocity kick from repositioning.

### Reposition Parameters (Phase B)
| Parameter | Default | Phase B status |
| --- | --- | --- |
| `BHRepositionMethod` | `0` | active: `0` off/diagnostic-only via verbosity, `1` drift, `2` teleport |
| `BHRepositionSearchRadius` | `3.0` | active (physical kpc) |
| `BHRepositionMaxDisplacement` | `0.5` | active (cell widths per update, method 1) |
| `BHRepositionDiagnosePotential` | `0` | active optional diagnostic |
| `BHRepositionVerbose` | `1` | active warning/log control |

Runtime validation:
- hard errors:
  - `BHRepositionMethod` not in `{0,1,2}`,
  - `BHRepositionSearchRadius <= 0`,
  - `BHRepositionMaxDisplacement < 0`.
- warnings:
  - under-resolved kernel (`search_radius/dx < 1.5`),
  - kernel wider than nominal ghost support (`search_radius/dx > 3`),
  - `BHRepositionMaxDisplacement` larger than kernel radius in cell widths,
  - potential diagnostic requested but potential field unavailable (once per level pass),
  - teleport displacement larger than search radius (method 2),
  - defensive clamp triggered after movement.

### `[BHREPOS]` Log Fields
Each line reports:
- BH id, redshift, BH position (pre-move for this pass), and host-cell density,
- diagnostic peak position/density, ghost flag, and BH-to-peak offset,
- active-zone peak position/density and offset,
- `active_target_exists`,
- optional potential-min offset (`-1` when unavailable/disabled),
- `displacement_kpc`, `displacement_cells`,
- `reposition_occurred`, `reposition_clamped`,
- `newly_seeded_skip`,
- kernel cell counts (`search_cells`, `search_active_cells`),
- `reposition_wall_ms`.

## BH Feedback Phase A (Diagnostics Only)
Phase A adds a feedback diagnostics handler that runs after accretion and before
star formation:

`BHRepositionDiagnosticHandler -> BHAccretionDiagnosticHandler -> BHFeedbackHandler -> StarParticleHandler`

What Phase A does:
- reads realized accretion diagnostics from the current timestep (`f_Edd`,
  `Mdot_actual`),
- classifies mode (`THERMAL` when `f_Edd > BHFeedbackModeThreshold`,
  otherwise `KINETIC`),
- computes diagnostic thermal energy and kinetic momentum budgets,
- computes diagnostic reservoir accumulation and would-be burst trigger,
- evaluates feedback-kernel geometry/temperature diagnostics,
- writes `[BHFDBK]` log lines.

What Phase A does **not** do:
- no thermal energy injection,
- no kinetic momentum injection,
- no baryon-field writes,
- no SF mask writes from feedback,
- no mutation of `BHFeedbackEnergyReservoir`,
- no mutation of `BHLastFeedbackRedshift`.

### Feedback Parameters (Phase A)
| Parameter | Default | Phase A status |
| --- | --- | --- |
| `BHFeedbackMethod` | `0` | active (`0` off, `1` thermal framework diagnostics, `2` two-mode scaffold diagnostics) |
| `BHFeedbackModeThreshold` | `0.01` | active |
| `BHFeedbackKernelRadius` | `1.0` | active (physical kpc) |
| `BHFeedbackThermalEfficiency` | `0.02` | active |
| `BHFeedbackMinEnergyBurst` | `1e50` | active diagnostic threshold (erg) |
| `BHFeedbackKineticEfficiency` | `0.1` | parsed (Phase C active) |
| `BHFeedbackWindVelocity` | `1e4` | parsed (km/s, Phase C active) |
| `BHFeedbackKineticGeometry` | `0` | parsed (Phase C active) |
| `BHFeedbackVerbose` | `1` | active |

Validation and warnings:
- hard errors:
  - `BHFeedbackMethod` not in `{0,1,2}`,
  - `BHFeedbackKernelRadius <= 0`,
  - `BHFeedbackThermalEfficiency` outside `[0,1]`,
  - `BHFeedbackMinEnergyBurst <= 0`,
  - `BHFeedbackKineticEfficiency < 0`,
  - `BHFeedbackWindVelocity <= 0`,
  - `BHFeedbackKineticGeometry` not in `{0,1}`.
- warnings:
  - `BHFeedbackMethod=2` reports that active two-mode deposition requires Phase C,
  - under-resolved kernel (`kernel_radius/dx < 1.5`),
  - kernel exceeding nominal ghost support (`kernel_radius/dx > 3`).

### Feedback Metadata Attributes
Phase A registers/checkpoints:
- `BHFeedbackEnergyReservoir` (initialized to `0`, not mutated in Phase A),
- `BHLastFeedbackRedshift` (initialized to `-1`, not mutated in Phase A).

It also stores:
- `BHLastMdotActual` (realized accretion rate from the accretion handler, code units).

### `[BHFDBK]` Fields
Each line reports:
- `step`, `level`, `z`, `bh_id`, `bh_mass`,
- `feedback_mode`, `f_Edd`, `L_feedback`,
- `E_requested`, `reservoir_before`, `reservoir_after`, `burst_diag`,
- `E_deposited` (always `0` in Phase A),
- `p_requested`, `p_deposited` (always `0` in Phase A),
- `feedback_kernel_cells`, `feedback_kernel_active_cells`,
- `feedback_kernel_gas_msun`,
- `T_before_mean`, `T_after_mean` (equal in Phase A),
- `n_sf_blocked_feedback` (always `0` in Phase A),
- `newly_seeded_skip`,
- `feedback_wall_ms`.

## BH Feedback Phase B (Thermal Coupling)
Phase B keeps the Phase A diagnostic pass, then activates thermal reservoir
mutation and burst deposition. It remains thermal-only. There is no kinetic
momentum injection, bipolar geometry, cooling suppression, temperature floor,
or feedback timestep limiter in Phase B.

Runtime order remains:

`BHRepositionDiagnosticHandler -> BHAccretionDiagnosticHandler -> BHFeedbackHandler -> StarParticleHandler`

### Phase B Mode Behavior
- `BHFeedbackMethod = 0`: feedback handler returns immediately.
- `BHFeedbackMethod = 1`: thermal feedback is active only when
  `f_Edd > BHFeedbackModeThreshold`.
- `BHFeedbackMethod = 2`: accepted for forward compatibility, but behaves like
  method 1 in Phase B. Low-`f_Edd` BHs are logged as `KINETIC_INACTIVE`.
- `KINETIC_INACTIVE`: no reservoir accumulation, no burst decision, no
  deposition, and the existing `BHFeedbackEnergyReservoir` value is preserved.
- Newly seeded BHs (`ParticleNumber == INT_UNDEFINED`) log
  `newly_seeded_skip=1` and do not feed back in the same pass.

### Burst Mechanics
For thermal-mode BHs, the handler reads `BHLastMdotActual` from the accretion
handler, computes the thermal request, and adds it to the checkpointed
`BHFeedbackEnergyReservoir` particle attribute:

```text
E_requested = BHFeedbackThermalEfficiency
            * BHAccretionRadiativeEfficiency
            * Mdot_actual
            * c^2
            * dt
```

The reservoir is stored in code energy units in the particle attribute, with
CGS erg conversion at the feedback handler's read/write boundaries. When
`BHFeedbackEnergyReservoir >= BHFeedbackMinEnergyBurst`, the full reservoir is
selected as one burst. A successful burst deposits the full stored energy and
sets the reservoir to zero. There is no partial drain or minimum temperature
jump criterion in Phase B.

If the burst kernel has no active-zone gas deposition cells, the burst energy
is returned to `BHFeedbackEnergyReservoir`, `E_deposited=0`,
`burst_occurred=0`, `kernel_gas_zero=1`, and `BHLastFeedbackRedshift` is not
updated.

### Thermal Deposition Pattern
Phase B matches ENZO's thermal feedback convention from `Grid_AddFeedbackSphere`
for `MBH_THERMAL` feedback:

- `TotalEnergy` and `InternalEnergy`/`GasEnergy` are specific energies.
- With `DualEnergyFormalism`, the internal-energy field is increased by the
  thermal specific-energy increment. For PPM-like hydro, `TotalEnergy` is then
  reconstructed as internal specific energy plus kinetic specific energy. For
  Zeus, the energy field is treated as internal specific energy.
- Without `DualEnergyFormalism`, `TotalEnergy` receives the thermal specific
  energy increment directly.

The feedback kernel is a sphere of radius `BHFeedbackKernelRadius` in physical
kpc. Only active-zone cells are deposition targets; ghost cells are never
written. Traversal is deterministic `k-j-i`. Deposition is mass-weighted:

```text
fraction_i = (rho_i * V_cell) / sum_j(rho_j * V_cell)
dE_i       = E_burst * fraction_i
de_i       = dE_i / (rho_i * V_cell)
```

The code applies a final-cell remainder correction and hard-fails if
`abs(sum(dE_i) - E_burst) / E_burst > 1e-10`.

### SF Blocking
Every active-zone cell with `dE_i > 0` is marked in the transient SF mask:

```text
BaryonField[NumberOfBaryonFields][cell] = 1.0
```

This is the same mask used by BH seeding/accretion and by star formation to
skip cells. It is a flag, so accretion-removed cells and feedback-heated cells
naturally form a union. Ghost cells and merely inspected kernel cells are not
blocked.

### Gas-Poor Handling
`BHFeedbackKernelMassWarnThreshold` is a Phase B parameter in `Msun`, default
`1e3`. If the active-zone gas mass in the feedback kernel is below this value,
the handler logs `kernel_gas_low=1` and deposits the full burst anyway. This is
warning-only and has no physics effect.

If no active-zone gas deposition cells are available, the handler logs
`kernel_gas_zero=1` and returns the full burst energy to the reservoir.

### Metadata and Restart
- `BHFeedbackEnergyReservoir` is checkpointed as a particle scalar and is the
  Phase B reservoir state. It is stored in code energy units, not erg.
- `BHLastFeedbackRedshift` is updated only when `E_deposited > 0`.
- On restart, accumulation continues from the stored reservoir value; burst
  timing is therefore restart-consistent up to the float precision of the
  particle attribute storage.

### Multi-BH Behavior
Within a grid, BH particles are processed in ascending particle ID order.
Overlapping kernels are additive: if two BHs heat the same active-zone cell, the
cell receives both thermal increments, and later BHs see the gas state left by
earlier BHs on that grid. The SF mask remains the union of all heated cells.
There is no Phase B inter-grid communication for feedback deposition.

### Phase B Feedback Parameters
| Parameter | Default | Phase B status |
| --- | --- | --- |
| `BHFeedbackMethod` | `0` | active (`0` off, `1` thermal active, `2` thermal active with kinetic inactive) |
| `BHFeedbackModeThreshold` | `0.01` | active thermal/kinetic-inactive threshold |
| `BHFeedbackKernelRadius` | `1.0` | active, physical kpc |
| `BHFeedbackEddingtonFactor` | `1.0` | Phase C multiplier on the existing Eddington cap |
| `BHFeedbackThermalEfficiency` | `0.02` | active thermal coupling efficiency |
| `BHFeedbackMinEnergyBurst` | `1e50` | active fixed burst threshold, erg |
| `BHFeedbackKernelMassWarnThreshold` | `1e3` | active gas-poor warning threshold, Msun |
| `BHFeedbackKineticEfficiency` | `0.1` | parsed but unused until Phase C |
| `BHFeedbackWindVelocity` | `1e4` | parsed but unused until Phase C |
| `BHFeedbackKineticGeometry` | `0` | parsed but unused until Phase C |
| `BHFeedbackVerbose` | `1` | active |

Validation hard-errors if `BHFeedbackMethod` is not in `{0,1,2}`,
`BHFeedbackKernelRadius <= 0`, `BHFeedbackEddingtonFactor < 0`,
`BHFeedbackThermalEfficiency` is outside `[0,1]`,
`BHAccretionRadiativeEfficiency * BHFeedbackThermalEfficiency > 1`,
`BHFeedbackMinEnergyBurst <= 0`, or `BHFeedbackKernelMassWarnThreshold <= 0`.
`BHFeedbackEddingtonFactor = 0` warns and prevents BH accretion; values above
`100` warn that the cap is effectively disabled.

### Phase B `[BHFDBK]` Fields
Phase B preserves the Phase A fields and adds:

- `reservoir_after_accum`: reservoir after adding this timestep's thermal
  request, before burst/zero-gas handling,
- `reservoir_final`: actual stored particle reservoir after all Phase B logic,
- `burst_occurred`: `1` only when energy was actually deposited,
- `kernel_gas_low`: `1` for gas-poor warning-only kernels,
- `kernel_gas_zero`: `1` when no active-zone gas cells received feedback,
- `dT_mean`: `T_after_mean - T_before_mean`,
- `sum_deposited`, `deposit_rel_err`, `mass_weight_rel_err`: conservation and
  mass-weighting diagnostics for Phase B tests.

`E_deposited`, `T_after_mean`, and `n_sf_blocked_feedback` are no longer
diagnostic placeholders in Phase B.

## BH Feedback Phase C (Eddington Factor and Reservoir Diagnostics)

Phase C adds diagnostics and one Eddington-control parameter. It does not add
new efficiency parameters and does not multiply another efficiency into the
Phase B feedback-energy formula.

### Eddington Factor
`BHFeedbackEddingtonFactor` multiplies the existing Eddington cap in the
accretion handler:

```text
Mdot_cap = BHFeedbackEddingtonFactor * Mdot_Edd
```

The default is `1.0`. At the default, the code uses the original Phase B cap
expression path rather than relying on a `1.0 * Mdot_Edd` multiply. Values below
`1.0` tighten the cap; values above `1.0` allow proportionally
super-Eddington requested accretion.

`lambda_edd` in `[BHACCR]` remains the raw accretion estimate divided by the
physical Eddington rate, not the factor-modified cap. `edd_factor` records the
runtime factor for traceability.

### Cumulative Reservoir Counters
Phase C adds two non-WINDS particle attributes:

- `BHCumulativeReservoirIn`: cumulative stored-code-energy added to the
  feedback reservoir.
- `BHCumulativeReservoirOut`: cumulative stored-code-energy actually removed
  from the reservoir after the burst path fully resolves.

Both counters are float-backed ParticleAttributes, the same storage class as
the reservoir. They are diagnostic bug detectors for catching missing,
double-counted, or wrong-sign updates. They are not high-precision long-run
energy-conservation accounting. Once a float32 cumulative counter is large
enough, small later increments can round away.

Newly created BH particles initialize both counters to `0.0` in
`Grid_BHSeedHandler.C` (`BHSeedCreateLocalBestParticle` and
`BHSeedCreateAcceptedCandidate`). For pre-Phase-C restarts where
`BHCumulativeReservoirIn` is missing, the reader initializes it to the current
`BHFeedbackEnergyReservoir` value so the conservation residual starts from zero
at the restart boundary. Missing `BHCumulativeReservoirOut` initializes to
`0.0`; pre-Phase-C history is not reconstructed.

### Conservation Residual
Each `[BHFDBK]` line appends:

- `cumul_reservoir_in_cgs`,
- `cumul_reservoir_out_cgs`,
- `conservation_residual_cgs`.

The residual is

```text
BHCumulativeReservoirIn - BHCumulativeReservoirOut - BHFeedbackEnergyReservoir
```

converted from stored code energy to erg for logging. A warning is emitted when
the absolute residual exceeds `1e-5` of the cumulative input. This threshold is
a practical implementation check; float32 accumulation can still drift in long
runs.

### Unit Conversion Path
`BHLastMdotActual` is written by the accretion handler in code mass per code
time. The feedback handler converts it as:

```text
MassUnits        = DensityUnits * LengthUnits^3
VelocityUnits    = LengthUnits / TimeUnits
EnergyUnits      = MassUnits * VelocityUnits^2
mdot_cgs         = Mdot_actual_code * MassUnits / TimeUnits
dt_cgs           = dt_code * TimeUnits
E_requested_erg  = epsilon_f * epsilon_r * mdot_cgs * c^2 * dt_cgs
```

For each heated cell:

```text
mcell_code       = rho_code * V_cell_code
mcell_cgs        = mcell_code * MassUnits
dE_i_erg         = E_burst_erg * mcell_code / sum(mcell_code)
de_specific_code = dE_i_erg / (mcell_code * EnergyUnits)
```

Constants used by the implementation:

- `clight = 2.99792458e10 cm/s`,
- `SolarMass = 1.9891e33 g`,
- `yr_s = 3.1556952e7 s`,
- `kboltz = 1.3806504e-16 erg/K`,
- `mh = 1.67262171e-24 g`,
- code `MassUnits = DensityUnits * LengthUnits^3`,
- code `LengthUnits`, `TimeUnits`, and `VelocityUnits` from `GetUnits`,
- code `EnergyUnits = MassUnits * VelocityUnits^2`.
