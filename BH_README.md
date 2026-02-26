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
2. Find gas cells that pass local gates:
   - dense enough,
   - cool enough,
   - metal-poor enough,
   - and optional extra checks (converging flow / thermal / self-bound).
3. Reject candidates too close to existing BHs using a linked-cell spatial hash.
4. Reject candidates that do not have enough gas mass for `BHSeedMass`.
5. Choose one deterministic winner globally (MPI-safe).
6. Create at most one MBH per level pass.
7. Print a `[BHSEED]` diagnostic line.

## BH Parameters (Beginner Table)
Defaults come from `src/enzo/SetDefaultGlobalValues.C`.

| Parameter | Default | What it means in plain language | Units |
| --- | --- | --- | --- |
| `BHSeedingMethod` | `0` | Main on/off switch (`0` = disabled, `1` = enabled). | none |
| `BHSeedOverdensityThreshold` | `100.0` | Gas density must be at least this value to be considered for a seed. | code density units |
| `BHSeedMetallicityThreshold` | `1e-4` | Gas metallicity must be below this value. Lower value = more metal-poor requirement. | absolute metal mass fraction |
| `BHSeedMetallicityThresholdInSolar` | `FLOAT_UNDEFINED` | Convenience input in solar units (`Z/Zsun`). Code converts using `Zsun = 0.02`. | solar metallicity units |
| `BHSeedMass` | `1e5` | Mass of each seeded BH particle. | Msun |
| `BHSeedTemperatureThreshold` | `1e4` | Gas temperature must be below this value. | K |
| `BHSeedExclusionRadius` | `100.0` | Minimum allowed distance from any existing BH/MBH. If candidate is closer, it is blocked. | **physical kpc** |
| `BHSeedVelDivCrit` | `1` | If `1`, require converging flow (`div(v) < 0`). | boolean-style int |
| `BHSeedThermalCrit` | `0` | If `1`, apply cooling-time criterion (`tcool` vs `tdyn`). | boolean-style int |
| `BHSeedSelfBoundCrit` | `0` | If `1`, require self-bound gas (`alpha < 1`). | boolean-style int |
| `BHSeedRunEveryTimestep` | `1` | If `1`, run every sub-cycle; if `0`, run with root-grid cadence logic. | boolean-style int |

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
BHSeedMetallicityThresholdInSolar = 1e-4
BHSeedExclusionRadius             = 50
BHSeedRunEveryTimestep            = 1
BHSeedVelDivCrit                  = 1
BHSeedThermalCrit                 = 1
BHSeedSelfBoundCrit               = 0
```

## Understanding `[BHSEED]` Log Output
You will see lines like:

```text
[BHSEED] step=... level=... z=... a_phys=... excl_phys_kpc=... excl_com_kpch=...
cell_code=... nbins=... ncand_local_min=... ncand_local_max=... ncand_global=...
ngates_density=... ngates_temp=... ngates_metal=... ngates_conv=... ngates_cool=...
ngates_bound=... ngates_mass=... dist_blocked=... created=... total_mbh=... pre_cache_bh=...
```

Key fields:
- `created`: number of BH seeds created in this level pass.
- `dist_blocked`: candidates rejected because they are too close to an existing BH.
- `ngates_mass`: candidates rejected because cell gas mass is less than `BHSeedMass`.
- `ncand_global`: total candidates after local gates.
- `total_mbh`: BH count in global cache after this pass.
- `pre_cache_bh`: BH count already in cache before this pass.

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

### Step 2: Generate Initial Conditions for Test
```bash
cd /mnt/home/boh10/bh_proj/enzo-dev/run/BHSeed/TS3_wrap
python3 make_ts3wrap_ic.py
```

What you should see:
- script prints created file info,
- `DD0000/data0000` is generated.

### What `DD0000/data0000` Contains (Simple Explanation)
This file is a tiny controlled setup so you can clearly see when BH seeding works.

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
mpirun -n 1 ../../src/enzo/enzo.exe -d bhseed_ts3wrap.enzo > ts3wrap_np1.log 2>&1
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
mpirun -n 4 ../../src/enzo/enzo.exe -d bhseed_ts3wrap.enzo > ts3wrap_np4.log 2>&1
grep '\[BHSEED\]' ts3wrap_np4.log
```

Expected result:
- same key behavior as single-rank run:
  - line 1 has `created=1`,
  - line 2 has `dist_blocked=2` and `created=0`.

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
