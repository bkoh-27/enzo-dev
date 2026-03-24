# Black Hole Seeding Guide (Enzo)

## What this guide is for
This guide is for someone who is new to black hole seeding in simulations and wants to quickly understand:
- what the seeding module does,
- what each parameter means,
- how to run a simple sanity test,
- how to read the log output.

## Big-picture workflow
At each level update, the code does this:
1. Build a list of existing BH particles.
2. Scan gas cells and keep only cells that pass local gates (density, temperature, metallicity, etc.).
3. For each surviving cell, evaluate gas in a spherical kernel around it.
4. Reject cells with too little enclosed gas mass.
5. Gather candidates across MPI ranks.
6. Sort candidates deterministically (same ranking rules every run).
7. Accept candidates in sorted order, while enforcing distance rules.
8. Create BH particles (up to `BHSeedMaxPerPass`).
9. Remove gas mass from active kernel cells in a conservative way.
10. Mark seeded cells so star formation skips those exact cells in the same update.

## Main files (for code readers)
- `src/enzo/star_maker_bh_seed.C` (local gates)
- `src/enzo/mbh_maker2.C` (candidate kernel evaluation)
- `src/enzo/Grid_MBHMaker2Handler.C` (grid wrapper)
- `src/enzo/Grid_BHSeedHandler.C` (global gather, ranking, acceptance, creation, logs)

## Parameter quick reference
Defaults come from `src/enzo/SetDefaultGlobalValues.C`.

| Parameter | Default | Plain-language meaning | Units |
| --- | --- | --- | --- |
| `BHSeedingMethod` | `1` | Global on/off switch (`1` on, `0` off). | none |
| `BHSeedRunEveryTimestep` | `0` | Run cadence (`1` every substep, `0` normal cadence). | bool-style int |
| `BHSeedOverdensityThreshold` | `1000.0` | Minimum gas density to be considered. | code density |
| `BHSeedTemperatureThreshold` | `1e4` | Maximum gas temperature allowed. | K |
| `BHSeedMetallicityThreshold` | `1e-4` | Maximum metallicity fraction allowed. | fraction |
| `BHSeedMetallicityThresholdInSolar` | `FLOAT_UNDEFINED` | Optional metallicity input in `Z/Zsun` (converted internally). | solar units |
| `BHSeedRequireFinestLevel` | `1` | If on, reject cells covered by finer AMR grids. | bool-style int |
| `BHSeedRequireLocalPeak` | `1` | If on, require 26-neighbor local density peak. | bool-style int |
| `BHSeedVelDivCrit` | `1` | If on, require converging flow (`div(v) < 0`). | bool-style int |
| `BHSeedThermalCrit` | `0` | If on, apply cooling/dynamical-time criterion. | bool-style int |
| `BHSeedSelfBoundCrit` | `0` | If on, apply self-bound criterion. | bool-style int |
| `BHSeedPatchRadius` | `3.0` | Spherical kernel radius for enclosed-property evaluation. | code length |
| `BHSeedMinEnclosedMass` | `1e6` | Minimum enclosed kernel gas mass required. | Msun |
| `BHSeedMass` | `1e5` | Mass assigned to each created BH particle. | Msun |
| `BHSeedLegacyCellMassGate` | `0` | `0`: shadow diagnostic only, `1`: hard reject if cell mass < seed mass. | bool-style int |
| `BHSeedExclusionMode` | `2` | Exclusion radius mode (`0` physical kpc, `1` comoving kpc/h, `2` cell-scaled). | enum |
| `BHSeedExclusionRadius` | `100.0` | Exclusion radius value (meaning depends on mode). | kpc or comoving kpc/h |
| `BHSeedExclusionCells` | `16` | Cell multiplier for mode `2` radius (`N * dx`). | cells |
| `BHSeedMinCandidateSeparation` | `3.0` | Minimum distance between accepted candidates in one pass. | physical kpc |
| `BHSeedMaxPerPass` | `10` | Maximum seeds created per level pass. | count |
| `BHSeedRankingOrder` | `0` | Ranking priority (`0`: enclosed-mass-first, `1`: peak-density-first). | enum |
| `BHSeedDeterministicTiebreak` | `1` | If on, use deterministic position tie-break. | bool-style int |
| `BHSeedChannel` | `0` | User tag written to seed metadata. | int |
| `BHSeedVerbose` | `1` | Log verbosity (`[BHSEED]`, `[BHSEED_SEED]`, optional debug). | int |

## Parameter explanations (simple, practical)

### 1) Gates: deciding whether a cell is even a candidate
- `BHSeedOverdensityThreshold`: raise this to reduce candidate count quickly.
- `BHSeedTemperatureThreshold`: lower values favor colder gas.
- `BHSeedMetallicityThreshold`: lower values favor metal-poor gas.
- `BHSeedRequireLocalPeak`: if `1`, only local maxima pass. This strongly cuts noisy candidates.
- `BHSeedRequireFinestLevel`: helps avoid parent-grid duplicates in AMR.

### 2) Kernel controls: deciding if surrounding gas is sufficient
- `BHSeedPatchRadius`: size of environment sampled around the candidate.
- `BHSeedMinEnclosedMass`: minimum enclosed gas mass in that kernel.
- Important: `BHSeedMinEnclosedMass` must be `>= BHSeedMass`.

### 3) Exclusion controls: preventing over-crowded seeding
- `BHSeedExclusionMode=0`: `BHSeedExclusionRadius` is physical kpc.
- `BHSeedExclusionMode=1`: `BHSeedExclusionRadius` is comoving kpc/h.
- `BHSeedExclusionMode=2`: radius is `BHSeedExclusionCells * dx_candidate`.
- `BHSeedMinCandidateSeparation`: minimum distance between accepted candidates in one pass.

### 4) Ranking and multiplicity
- `BHSeedRankingOrder=0`: prioritize larger enclosed mass.
- `BHSeedRankingOrder=1`: prioritize higher peak density.
- `BHSeedMaxPerPass`: hard cap on how many accepted candidates are created.
- `BHSeedDeterministicTiebreak`: keeps ordering reproducible across runs.

### 5) Legacy cell-mass switch
- `BHSeedLegacyCellMassGate=0`: only records a diagnostic counter.
- `BHSeedLegacyCellMassGate=1`: turns that old check into a hard reject.

### 6) Logging controls
- `BHSeedVerbose=0`: minimal structured log.
- `BHSeedVerbose=1`: adds per-seed metadata lines.
- `BHSeedVerbose>=2`: adds debug conservation lines.

## Runtime safety checks
The code aborts at startup if:
- `BHSeedMinEnclosedMass < BHSeedMass`

This protects fixed-mass gas removal from impossible parameter combinations.

## Metadata written to BH particles
Each created seed stores these attributes:
- `bhseed_channel`
- `bhseed_redshift`
- `bhseed_patch_mass`
- `bhseed_patch_metallicity`
- `bhseed_patch_density_peak`
- `bhseed_kernel_complete`
- `bhseed_host_dm_density`
- `bhseed_accept_rank`

Quick check example:
```bash
h5dump -d "/Grid00000001/bhseed_accept_rank" DD0001/data0001.cpu0000
```

## Recommended explicit settings
Do not rely on defaults in production runs. Set values explicitly in your `.enzo` file.

```ini
BHSeedingMethod                   = 1
BHSeedRunEveryTimestep            = 0
BHSeedOverdensityThreshold        = 1000
BHSeedTemperatureThreshold        = 1e4
BHSeedMetallicityThreshold        = 1e-4
BHSeedMetallicityThresholdInSolar = 1e-4
BHSeedRequireFinestLevel          = 1
BHSeedRequireLocalPeak            = 1
BHSeedVelDivCrit                  = 1
BHSeedThermalCrit                 = 0
BHSeedSelfBoundCrit               = 0
BHSeedPatchRadius                 = 3.0
BHSeedMinEnclosedMass             = 1e6
BHSeedMass                        = 1e5
BHSeedLegacyCellMassGate          = 0
BHSeedExclusionMode               = 2
BHSeedExclusionRadius             = 100.0
BHSeedExclusionCells              = 16
BHSeedMinCandidateSeparation      = 3.0
BHSeedMaxPerPass                  = 10
BHSeedRankingOrder                = 0
BHSeedDeterministicTiebreak       = 1
BHSeedChannel                     = 0
BHSeedVerbose                     = 1
```

## Reading `[BHSEED]` logs
Typical line contains counters like:
- `ncand_global`: candidates after local gates.
- `dist_blocked`: candidates blocked by existing BH distance rule.
- `ngates_enclosedmass`: candidates failing enclosed-mass gate.
- `nseeds_created`: seeds actually created this pass.
- `nseeds_skipped_insufficient_gas`: accepted candidates skipped at creation because active-zone removable gas was insufficient.
- `seeding_wall_ms`: seeding wall-clock time for that pass.

When `BHSeedVerbose >= 1`, each created seed also prints `[BHSEED_SEED]` with position, kernel metrics, and `accept_rank`.

## Minimal sanity test (TS3_wrap)

### 1) Build
```bash
cd /mnt/home/boh10/bh_proj/enzo-dev/src/enzo
module load gcc openmpi
make -j8
```

### 2) Run 1-rank and 4-rank
```bash
cd /mnt/home/boh10/bh_proj/enzo-dev/run/BHSeed/TS3_wrap
mpirun -n 1 ../../../src/enzo/enzo.exe -d bhseed_ts3wrap.enzo > ts3wrap_np1.log 2>&1
mpirun -n 4 ../../../src/enzo/enzo.exe -d bhseed_ts3wrap.enzo > ts3wrap_np4.log 2>&1
bash assert_ts3wrap.sh
```

Expected behavior:
- first BH pass creates one seed,
- next pass blocks nearby candidates by exclusion,
- assertion script prints `TS3_WRAP: PASS`.

## Troubleshooting
- No `[BHSEED]` lines:
  - check `BHSeedingMethod = 1`.
- Too many candidates:
  - raise density threshold and keep local peak on.
- Too many exclusion blocks:
  - re-check exclusion mode and exclusion radius units.
- Frequent `nseeds_skipped_insufficient_gas`:
  - lower `BHSeedMass`, reduce kernel radius, or improve local active-zone gas support.
