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
- `src/enzo/mbh_maker2.C` (distance/mass filtering and local winner)
- `src/enzo/Grid_MBHMaker2Handler.C` (grid wrapper)
- `src/enzo/Grid_BHSeedHandler.C` (global cache, MPI winner, logging)

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
8. Choose one deterministic winner globally (MPI-safe, still one winner per pass).
9. Create at most one MBH per level pass and write per-seed metadata.
10. Mark the seeded cell so star formation skips that exact cell in the same pass.
11. Print a `[BHSEED]` diagnostic line (and `[BHSEED_SEED]` lines when verbose logging is on).

## BH Parameters (Beginner Table)
Defaults come from `src/enzo/SetDefaultGlobalValues.C`.

| Parameter | Default | What it means in plain language | Units |
| --- | --- | --- | --- |
| `BHSeedingMethod` | `1` | `int`, **Phase 1 active**. Main on/off switch (`0` = disabled, `1` = enabled). | none |
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
| `BHSeedExclusionMode` | `2` | `int`, **parsed-only placeholder (Phase 3+)**. Future exclusion algorithm selector. | mode id |
| `BHSeedExclusionRadius` | `100.0` | `float`, **Phase 1 active**. Minimum allowed distance from any existing BH/MBH. If candidate is closer, it is blocked. | **physical kpc** |
| `BHSeedExclusionCells` | `16` | `int`, **parsed-only placeholder (Phase 3+)**. Future cell-based exclusion window. | cells |
| `BHSeedMinCandidateSeparation` | `3.0` | `float`, **parsed-only placeholder (Phase 3+)**. Future candidate de-dup separation setting. | cells |
| `BHSeedMaxPerPass` | `10` | `int`, **parsed-only placeholder (Phase 3+)**. Future per-pass multi-seed cap. | count |
| `BHSeedRunEveryTimestep` | `0` | `int (bool-style)`, **Phase 1 active**. If `1`, run every sub-cycle; if `0`, run with root-grid cadence logic. | boolean-style int |
| `BHSeedRankingOrder` | `0` | `int`, **parsed-only placeholder (Phase 3+)**. Future ranking-order selector for tie resolution. | enum id |
| `BHSeedVerbose` | `1` | `int`, **Phase 2 active**. Controls BH seeding log verbosity (`[BHSEED]` + `[BHSEED_SEED]` kernel metadata line). | verbosity level |
| `BHSeedDeterministicTiebreak` | `1` | `int (bool-style)`, **parsed-only placeholder (Phase 3+)**. Reserved deterministic tie-break selector. | boolean-style int |

### Important Unit Note for `BHSeedExclusionRadius`
You set it in physical kpc.

In cosmological runs, code converts internally for comparison:
- `R_comoving_kpc/h = R_physical_kpc * h / a_phys`

The runtime log shows both values so you can verify conversion:
- `excl_phys_kpc`
- `excl_com_kpch`

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
BHSeedExclusionMode               = 2        # Phase 3+ parsed-only in Phase 1
BHSeedExclusionRadius             = 100.0
BHSeedExclusionCells              = 16       # Phase 3+ parsed-only in Phase 1
BHSeedMinCandidateSeparation      = 3.0      # Phase 3+ parsed-only in Phase 1
BHSeedMaxPerPass                  = 10       # Phase 3+ parsed-only in Phase 1
BHSeedRunEveryTimestep            = 0
BHSeedRankingOrder                = 0        # Phase 3+ parsed-only in Phase 1
BHSeedVelDivCrit                  = 1
BHSeedLegacyCellMassGate          = 0
BHSeedThermalCrit                 = 0
BHSeedSelfBoundCrit               = 0
BHSeedRequireFinestLevel          = 1
BHSeedRequireLocalPeak            = 1
BHSeedVerbose                     = 1
BHSeedDeterministicTiebreak       = 1        # Phase 3+ parsed-only in Phase 1
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

Log compatibility note:
- New fields are appended at the end of the existing `[BHSEED]` line after all legacy and Phase 1 fields.

When `BHSeedVerbose >= 1`, each created seed also emits:

```text
[BHSEED_SEED] level=... x=... y=... z=... channel=... redshift=... patch_mass=...
patch_metal=... patch_density_peak=... kernel_complete=... host_dm_density=... accept_rank=...
```

This line reports per-seed metadata values written into particle attributes/HDF5.
`patch_density_peak` is the kernel maximum (not the candidate-cell density placeholder used in Phase 1).

## Seed Metadata Fields (Phase 2)
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
| `BHSeedAcceptRank` | `bhseed_accept_rank` | `11 / 14` | `-1.0` (`-1` semantic) | placeholder (Phase 3+) |

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
