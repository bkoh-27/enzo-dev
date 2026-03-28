# BH Seeding User Guide

## Purpose
This guide explains the BH seeding module for simulation users who are not familiar with the implementation.

Use this document to:
- understand what BH seeding does,
- set parameters with the right units,
- verify that seeding is working,
- interpret `[BHSEED]` and `[BHSEED_SEED]` logs,
- troubleshoot common setup issues.

This guide is focused only on seeding. Accretion is documented in `BH_ACCRETION_README.md`.

## Where Seeding Happens in the Code
Main files:
- `src/enzo/star_maker_bh_seed.C` (local gate checks)
- `src/enzo/mbh_maker2.C` (candidate kernel evaluation)
- `src/enzo/Grid_MBHMaker2Handler.C` (grid wrapper)
- `src/enzo/Grid_BHSeedHandler.C` (global gather/sort/accept/create)

Execution order at each level update (simplified):
1. local candidate screening on each grid,
2. global candidate gather and deterministic ranking,
3. acceptance walk with exclusion and de-duplication,
4. seed creation and kernel-distributed gas removal,
5. continue to accretion/star formation modules.

## High-Level Seeding Algorithm
For each level pass:
1. Build local candidate list from gas cells.
2. Apply local gates in fixed order:
- overdensity threshold
- local density peak requirement (optional)
- temperature threshold
- metallicity threshold
- optional flow/thermal/self-bound criteria
3. Apply legacy cell-mass gate behavior:
- `BHSeedLegacyCellMassGate = 0`: shadow-only diagnostic
- `BHSeedLegacyCellMassGate = 1`: hard reject
4. Reject cells too close to existing BHs.
5. Evaluate spherical patch properties for survivors:
- enclosed gas mass
- patch metallicity
- patch density peak
- kernel completeness
6. Reject candidates with enclosed mass below `BHSeedMinEnclosedMass`.
7. Gather all surviving candidates across MPI ranks.
8. Deterministically sort candidates (`BHSeedRankingOrder` + tiebreak).
9. Run acceptance walk:
- reject by exclusion against existing BHs
- reject by within-pass separation from already accepted candidates
- stop at `BHSeedMaxPerPass`
10. Create accepted seeds and remove gas with active-zone-only, proportional kernel removal.
11. Set SF mask on created-seed cells so star formation skips those cells in the same pass.
12. Emit structured logs.

## Parameter Quick Reference
Defaults are from `src/enzo/SetDefaultGlobalValues.C`.

### Core toggles
| Parameter | Default | Meaning | Units |
| --- | --- | --- | --- |
| `BHSeedingMethod` | `1` | Master switch (`0` off, `1` on). | none |
| `BHSeedRunEveryTimestep` | `0` | `0` root-like cadence, `1` every sub-step. | bool int |
| `BHSeedVerbose` | `1` | Verbosity for `[BHSEED]` and `[BHSEED_SEED]`. | level |

### Local gate thresholds
| Parameter | Default | Meaning | Units |
| --- | --- | --- | --- |
| `BHSeedOverdensityThreshold` | `1000.0` | Minimum density gate. | code density |
| `BHSeedTemperatureThreshold` | `1e4` | Maximum temperature gate. | K |
| `BHSeedMetallicityThreshold` | `1e-4` | Maximum metallicity gate. | mass fraction |
| `BHSeedMetallicityThresholdInSolar` | `FLOAT_UNDEFINED` | Optional metallicity input in `Z/Zsun`. | solar units |
| `BHSeedVelDivCrit` | `1` | Require converging flow if enabled. | bool int |
| `BHSeedThermalCrit` | `0` | Apply cooling-time criterion if enabled. | bool int |
| `BHSeedSelfBoundCrit` | `0` | Require self-bound gas if enabled. | bool int |
| `BHSeedRequireLocalPeak` | `1` | Require 26-neighbor local density peak. | bool int |
| `BHSeedRequireFinestLevel` | `1` | Reject parent-level cells covered by finer AMR grids. | bool int |

### Kernel and mass controls
| Parameter | Default | Meaning | Units |
| --- | --- | --- | --- |
| `BHSeedPatchRadius` | `3.0` | Spherical patch radius used for enclosed checks. | code length |
| `BHSeedMinEnclosedMass` | `1e6` | Minimum enclosed gas mass in patch. | Msun |
| `BHSeedMass` | `1e5` | Seed particle mass at creation. | Msun |
| `BHSeedLegacyCellMassGate` | `0` | Legacy cell-mass gate mode (`0` shadow, `1` hard gate). | bool int |

### Exclusion and multiplicity controls
| Parameter | Default | Meaning | Units |
| --- | --- | --- | --- |
| `BHSeedExclusionMode` | `2` | `0` physical-kpc, `1` comoving-kpc/h, `2` cell-scaled. | enum |
| `BHSeedExclusionRadius` | `100.0` | Radius input for mode 0 or 1. | mode-dependent |
| `BHSeedExclusionCells` | `16` | Radius multiplier for mode 2. | cells |
| `BHSeedMinCandidateSeparation` | `3.0` | Min distance between accepted candidates in a pass. | physical kpc |
| `BHSeedMaxPerPass` | `10` | Max seeds accepted per level pass. | count |
| `BHSeedRankingOrder` | `0` | `0` mass-first, `1` density-first ranking. | enum |
| `BHSeedDeterministicTiebreak` | `1` | Enable deterministic tie resolution. | bool int |

### Metadata field
| Parameter | Default | Meaning | Units |
| --- | --- | --- | --- |
| `BHSeedChannel` | `0` | User tag stored on each seed particle. | integer tag |

## Important Unit Notes
Exclusion radius units depend on `BHSeedExclusionMode`:
- Mode `0`: `BHSeedExclusionRadius` is physical kpc.
- Mode `1`: `BHSeedExclusionRadius` is comoving kpc/h.
- Mode `2`: exclusion radius is `BHSeedExclusionCells * dx_local`.

The `[BHSEED]` log prints both physical and comoving radius forms for clarity.

## Runtime Validation
Seeding aborts at startup if:
- `BHSeedMinEnclosedMass < BHSeedMass`

Reason: fixed seed mass cannot be guaranteed from the removal kernel in that case.

## Recommended Explicit Parameter Block
Use explicit settings in your `.enzo` file instead of relying on defaults:

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

## How to Read `[BHSEED]` Logs
Example fields you should monitor:
- candidate flow:
  - `ncand_global`
  - `ncandidates_gathered`
- gate diagnostics:
  - `ngates_density`, `ngates_temp`, `ngates_metal`
  - `ngates_peak`, `ngates_finestlevel`
  - `ngates_enclosedmass`
  - `nlegacy_cellmass_would_fail`
- acceptance and creation:
  - `created` and `nseeds_created`
  - `ncandidates_dedup_rejected`
  - `ncandidates_exclusion_rejected`
  - `walk_stopped_at_max`
  - `nseeds_skipped_insufficient_gas`
- performance:
  - `seeding_wall_ms`

Per-seed line (`[BHSEED_SEED]`) includes:
- position
- channel
- seed redshift
- patch mass/metallicity/density peak
- kernel completeness
- host DM density
- global acceptance rank

## Seed Metadata Written to Particles
Seeding stores particle attributes that are checkpoint-safe and available in HDF5 outputs.

Common datasets:
- `bhseed_channel`
- `bhseed_redshift`
- `bhseed_patch_mass`
- `bhseed_patch_metallicity`
- `bhseed_patch_density_peak`
- `bhseed_kernel_complete`
- `bhseed_host_dm_density`
- `bhseed_accept_rank`

Quick check:
```bash
h5dump -d "/Grid00000001/bhseed_accept_rank" DD0001/data0001.cpu0000
```

## Quick Validation Workflow (TS3_wrap)
1. Build Enzo.
2. Run TS3_wrap with `np=1` and `np=4`.
3. Confirm logs show:
- first relevant pass: `created=1`
- follow-up pass: exclusion blocks the second seed as expected
- deterministic counters across MPI layouts

Useful commands:
```bash
mpirun -n 1 ../../../src/enzo/enzo.exe -d bhseed_ts3wrap.enzo > ts3wrap_np1.log 2>&1
mpirun -n 4 ../../../src/enzo/enzo.exe -d bhseed_ts3wrap.enzo > ts3wrap_np4.log 2>&1
rg "\[BHSEED\]|\[BHSEED_SEED\]" ts3wrap_np1.log ts3wrap_np4.log
```

## Troubleshooting
### No `[BHSEED]` lines
- Check `BHSeedingMethod = 1`.
- Verify log capture and grep pattern.

### No seeds created
- Thresholds may be too strict.
- Check `ngates_*` counters to identify the blocking gate.
- Check enclosed-mass and legacy cell-mass diagnostics.

### Too many seeds
- Increase exclusion/separation or tighten local gates.
- Check `BHSeedMaxPerPass`.

### Many accepted candidates skipped for insufficient gas
- `nseeds_skipped_insufficient_gas` high means removal kernel lacks mass.
- Consider smaller `BHSeedMass`, smaller patch/removal demand, or better local resolution.

## Practical Tuning Advice
- If candidate count is huge, first tighten local gates.
- If seeds cluster too tightly, adjust exclusion mode/radius and candidate separation.
- If seeding is too sparse, relax local gates before changing ranking order.
- Keep `BHSeedMinEnclosedMass >= BHSeedMass` always.

## Interaction With Accretion Module
Seeding and accretion are separate modules.
- Seeding creates BH particles.
- Accretion updates BH growth and gas consumption on later passes.

Accretion setup and logs are documented in `BH_ACCRETION_README.md`.
