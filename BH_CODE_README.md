# BH_CODE_README.md

## Why This File Exists
This is the technical companion to `BH_README.md`.

Use this file if you want to know:
- exactly which source files were changed,
- what each change does,
- why each change was needed,
- where to look in code when debugging.

## Quick Summary
The BH seeding feature was added as a standalone pipeline that runs before normal star formation.

Core behavior:
1. Build a list of BH seed candidates on each grid.
2. Apply local gates (density, temperature, metallicity, optional flow/thermal/bound checks).
3. Apply global distance exclusion against existing BHs using linked-cell spatial hashing.
4. Select one deterministic winner across MPI ranks.
5. Create one MBH particle (max one per level pass) and log diagnostics.

## Code Flow (Start Here)
If you are reading code for the first time, use this order:

1. `src/enzo/EvolveLevel.C`
   - Calls `BHSeedBeginLevel(...)`
   - Calls `grid::MBHMaker2Handler(...)` on each grid
   - Calls `BHSeedFinalizeLevel()`

2. `src/enzo/Grid_MBHMaker2Handler.C`
   - Prepares required fields (`temperature`, `cooling_time`, `dmfield`, `metal_fraction`)
   - Calls `mbh_maker2(...)`

3. `src/enzo/mbh_maker2.C`
   - Calls local kernel `star_maker_bh_seed(...)`
   - Applies distance and mass gates
   - Reports local best candidate to global state

4. `src/enzo/star_maker_bh_seed.C`
   - Implements local gates and per-gate counters
   - Sorts candidates deterministically

5. `src/enzo/Grid_BHSeedHandler.C`
   - Handles global BH cache, spatial bins, MPI winner, MBH creation, and `[BHSEED]` logging

## New Source Files (What Each One Does)

### `src/enzo/star_maker_bh_seed.C`
Purpose:
- Candidate kernel (Gates 1-7).
- Produces:
  - candidate indices,
  - candidate densities,
  - diagnostic counters (`ngates_density/temp/metal/conv/cool/bound`).

Important detail:
- Candidates are sorted by density descending.
- Tie-break is linear cell index ascending (deterministic).

Gates 1-7 (local candidate gates) are:
1. Finest-level gate:
   - cell is skipped if it is covered by a finer grid (`r[idx] != 0`).
2. Density gate:
   - reject if `density < BHSeedOverdensityThreshold`.
3. Temperature gate:
   - reject if `temperature > BHSeedTemperatureThreshold`.
4. Metallicity gate:
   - reject if `metallicity > BHSeedMetallicityThreshold`.
5. Converging-flow gate (optional):
   - active only when `BHSeedVelDivCrit = 1`,
   - reject if `div(v) >= 0`.
6. Thermal/cooling gate (optional):
   - active only when `BHSeedThermalCrit = 1`,
   - reject if thermal criterion fails (`tdyn < tcool` in this implementation).
7. Self-bound gate (optional):
   - active only when `BHSeedSelfBoundCrit = 1`,
   - reject if bound criterion fails (`alpha >= 1`).

### `src/enzo/mbh_maker2.C`
Purpose:
- Takes kernel candidates and applies:
  - Gate 8: distance exclusion,
  - Gate 9: mass availability (`cell gas mass >= BHSeedMass`).
- Sends one local winner candidate to global state.

Important detail:
- Uses `BHSeedCandidateBlocked(...)` from global handler (linked-cell lookup).
- Counts `dist_blocked` and `ngates_mass`.

### `src/enzo/Grid_MBHMaker2Handler.C`
Purpose:
- Grid-level adapter layer.
- Computes fields and units needed by `mbh_maker2(...)`.

Why it exists:
- Keeps grid data extraction separate from BH core logic.

### `src/enzo/Grid_BHSeedHandler.C`
Purpose:
- Owns global BH seeding state and level lifecycle.
- Provides:
  - cache build/update,
  - linked-cell bin build,
  - deterministic MPI winner selection,
  - MBH particle creation,
  - `[BHSEED]` logging.

## Existing Files Changed (and Why)

| File | Change | Why |
| --- | --- | --- |
| `src/enzo/EvolveLevel.C` | Added begin/per-grid/finalize BH calls | Integrates BH seeding into timestep flow |
| `src/enzo/Grid.h` | Added `MBHMaker2Handler`, friend declarations, `ReturnMakeStars()` | Exposes required APIs for BH path |
| `src/enzo/Make.config.objects` | Added BH object files | Ensures new BH files compile and link |
| `src/enzo/global_data.h` | Added BH parameters + `Outfptr` extern | Global config and logging access |
| `src/enzo/SetDefaultGlobalValues.C` | Added BH defaults | Predictable default behavior |
| `src/enzo/ReadParameterFile.C` | Added BH parsing, validation, metallicity-in-solar conversion | User config + safety checks |
| `src/enzo/WriteParameterFile.C` | Added BH parameter writeout with units note | Restart continuity and provenance |
| `src/enzo/InitializeNew.C` | Particle attributes path and `Outfptr` lifecycle cleanup | Prevents missing attributes and stale pointer usage |
| `src/enzo/ReadAllData.C` | Include BH path in attribute allocation logic | Restart/read compatibility |
| `src/enzo/Group_ReadAllData.C` | Same as above for grouped IO | Restart/read compatibility |
| `src/enzo/Grid_StarParticleHandler.C` | Formatting-only differences | No functional BH logic change |

## Parameter Semantics (Current)
- `BHSeedExclusionRadius` input is **physical kpc**.
- In cosmological runs, runtime conversion is:
  - `R_comoving_kpc/h = R_physical_kpc * h / a_phys`
- `BHSeedMetallicityThresholdInSolar` is converted using `Zsun = 0.02`.

Current defaults (from `SetDefaultGlobalValues.C`):
- `BHSeedingMethod = 0`
- `BHSeedOverdensityThreshold = 100.0`
- `BHSeedMetallicityThreshold = 1e-4`
- `BHSeedMass = 1e5`
- `BHSeedTemperatureThreshold = 1e4`
- `BHSeedExclusionRadius = 100.0` (physical kpc)
- `BHSeedVelDivCrit = 1`
- `BHSeedThermalCrit = 0`
- `BHSeedSelfBoundCrit = 0`
- `BHSeedRunEveryTimestep = 1`

## Important Correctness Fixes in This Implementation

### 1) Correct scale-factor usage for exclusion conversion
Problem:
- Enzo internal `A` is not direct physical `a`.

Fix:
- Convert using `a_phys = A / (1 + InitialRedshift)` before `h/a_phys` scaling.

Impact:
- Prevents systematic exclusion-radius mis-scaling in cosmological runs.

### 2) Cache includes BHs from all AMR levels
Problem:
- Current-level-only scans can miss BHs on other levels.

Fix:
- Initial cache fill scans `LevelArray[0..MaximumRefinementLevel]`.

Impact:
- Distance exclusion now sees all existing BHs, not only same-level BHs.

### 3) Logging made diagnostic-friendly
Added fields:
- `a_phys`
- `excl_phys_kpc`
- `excl_com_kpch`
- `pre_cache_bh`

Impact:
- Unit conversion and cache state are visible in logs.

## Performance Improvements in This Implementation

### A) Linked-cell rebuild dirty-check
Change:
- Rebuild bins only if bin counts changed or BH count changed.

Why:
- Avoids repeated large memory clears in unchanged states.

### B) One-time full BH cache initialization
Change:
- Full all-level BH scan happens once, then cache is updated incrementally as seeds are created.

Why:
- Avoids repeated AMR-wide particle scans.

### C) Full-box short-circuit for very large exclusion radius
Change:
- If exclusion radius exceeds maximum periodic-box separation, return blocked immediately.

Why:
- Skips unnecessary linked-cell traversal in large-radius suppression runs.

## Determinism Guarantees
- Local candidate order is deterministic.
- Global winner tie-break is deterministic:
  1. highest density,
  2. smallest `(x, y, z)`,
  3. smallest rank ID.
- Winner creation is synchronized with `MPI_Bcast`.

Result:
- Tested behavior matches across `np=1` and `np=4` for validation cases.

## What Was Intentionally Not Changed
To keep scope tight and avoid unrelated risk:
- No changes to `src/enzo/star_maker2.F`.
- No MPI collectives added inside per-grid inner loops.
- No broad refactors in unrelated star formation code.

## User-Facing Test Assets Included
Directory:
- `run/BHSeed/TS3_wrap/`

Files:
- `bhseed_ts3wrap.enzo`
- `make_ts3wrap_ic.py`
- `assert_ts3wrap.sh`

Purpose:
- Reproducible periodic-wrap regression test with deterministic `np=1` vs `np=4` checks.

## Practical Debug Checklist
When results are unexpected, check `[BHSEED]` fields in this order:
1. `pre_cache_bh` (did cache see existing BHs?)
2. `ncand_global` (are any cells passing local gates?)
3. `ngates_density/temp/metal/conv/cool/bound` (which gate is removing candidates?)
4. `dist_blocked` (is exclusion active?)
5. `ngates_mass` (is seed mass too large for candidate cells?)
6. `created` and `total_mbh` (was a seed actually created?)
