# BH Accretion User Guide

## Purpose
This guide explains the BH accretion module for simulation users who are not familiar with the implementation details.

Use this document to:
- understand what the accretion module does each timestep,
- configure the key parameters safely,
- interpret `[BHACCR]` and `[BHACCR_WARN]` log lines,
- verify mass-growth bookkeeping,
- troubleshoot common run-time behavior.

This guide focuses only on accretion. Seeding is documented in `BH_SEEDING_README.md`.

## Is This Separate From Legacy Enzo BH Accretion?
Yes.

This module is the new BH accretion path controlled by `BHAccretion*` parameters and `[BHACCR]` logs.
It is separate from older legacy MBH controls such as `MBHAccretion*`.

Quick way to confirm your run is using this module:
- set `BHAccretionMethod = 1`,
- check logs for `[BHACCR]` lines,
- verify parameters in your `.enzo` file use `BHAccretion*` names.

## Where Accretion Happens in the Code
Main files:
- `src/enzo/Grid_BHAccretionHandler.C` (diagnostics + active accretion)
- `src/enzo/star_maker_bh_seed.C` and `src/enzo/Grid_BHSeedHandler.C` (seeding runs first)
- `src/enzo/Grid_StarParticleHandler.C` (reads SF mask to suppress star formation in flagged cells)

Execution order at each level update (simplified):
1. seeding finalization,
2. BH accretion handler,
3. star formation and other active particle routines.

## High-Level Accretion Algorithm (Plain Language)
For each BH particle on a grid:
1. Build a diagnostic spherical kernel (`BHAccretionKernelRadius`, physical kpc).
2. Split kernel gas into hot/cold using `t_cool / t_dyn`.
- if cooling time is unavailable/unreliable for a cell, fallback uses temperature split (`BHAccretionTSplitFloor`).
3. Compute raw hot/cold rates, then apply Eddington cap.
4. Convert capped rate to requested mass this step:
- `dm_requested = Mdot_total_capped * dt`.
5. Build a removal kernel (`BHAccretionRemovalRadius`, cell widths) and remove gas proportionally from active-zone cells only.
- optional debug mode: host cell only (`BHAccretionRemovalMode = 1`).
6. Cap each cell's removal so density, internal energy, and pressure remain valid.
7. Grow BH mass by exactly the removed mass (`dm_removed`) and update metadata.

Important behavior:
- Newly seeded BHs are skipped for accretion in the same pass.
- They seed first, then accrete starting from the next pass.
- If the diagnostic kernel exceeds ghost-zone/grid support, active gas removal
  is rejected and logged (`accretion_source_rejected=1`); diagnostics are still
  emitted.

## Phase A vs Phase B (What Is Active)
Phase A behavior (still present):
- diagnostic hot/cold classification,
- raw and capped rate calculation,
- Eddington ratio diagnostics.

Phase B additions (active):
- real gas removal,
- BH mass growth,
- realized channel rates after gas limitation,
- per-cell density/energy/pressure floors and source-update diagnostics,
- SF blocking for removal cells,
- momentum-omission diagnostics and warnings,
- invariant checks and bookkeeping enforcement.

## Key Parameters (Most Users)
Defaults are from `src/enzo/SetDefaultGlobalValues.C`.

| Parameter | Default | Meaning | Units |
| --- | --- | --- | --- |
| `BHAccretionMethod` | `0` | Module switch (`0` off, `1` on). The new `[BHACCR]` source path is opt-in. | enum |
| `BHAccretionKernelRadius` | `3.0` | Diagnostic kernel radius. | physical kpc |
| `BHAccretionRemovalRadius` | `1` | Removal kernel radius. | cell widths |
| `BHAccretionRemovalMode` | `0` | `0` multi-cell, `1` host-cell debug. | enum |
| `BHAccretionTSplitFloor` | `5e5` | Temperature fallback split for hot/cold. | K |
| `BHAccretionCVisc` | `6.283` | Angular-momentum suppression scale; larger values suppress cold accretion more strongly. | dimensionless |
| `BHAccretionNHStar` | `0.1` | Density threshold for hot-channel boost. | cm^-3 |
| `BHAccretionBeta` | `1.0` | Slope of hot boost factor. | dimensionless |
| `BHAccretionAlphaMax` | `10.0` | Maximum hot boost cap. | dimensionless |
| `BHAccretionRadiativeEfficiency` | `0.1` | `epsilon_r` for Eddington rate. | dimensionless |
| `BHAccretionIgnoredDVWarn` | `1.0` | Momentum-omission warning threshold. | km/s |
| `BHAccretionIgnoredPFracWarn` | `0.01` | Momentum-omission warning threshold. | dimensionless |
| `BHAccretionRunEveryTimestep` | `0` | Cadence control. | bool int |
| `BHAccretionVerbose` | `1` | Logging verbosity. | level |

Reserved/not active for this phase:
- `BHAccretionSuperEddington`,
- `BHAccretionSuperEddFactor`,
- `BHAccretionUseReservoir`.

## Unit and Radius Notes
- Diagnostic kernel radius (`BHAccretionKernelRadius`) is in physical kpc.
- Removal kernel radius (`BHAccretionRemovalRadius`) is in local cell widths.
- `BHAccretionIgnoredDVWarn` is configured in km/s.
- BH internal mass bookkeeping is in code mass units.
- Log output includes converted rates in Msun/yr and many cgs diagnostics.

## Three-Tier Bookkeeping (Core Concept)
This module tracks accretion in three tiers:
1. Raw physics rates:
- `Mdot_hot_raw`, `Mdot_cold_raw`, `Mdot_total_raw`.
2. Eddington-capped rates:
- `frac_cap`, `Mdot_hot_capped`, `Mdot_cold_capped`, `dm_requested`.
3. Gas-limited realized rates:
- `frac_gas`, `dm_removed`, `Mdot_hot_realized`, `Mdot_cold_realized`.

BH growth is always tied to Tier 3 (`dm_removed`), not Tier 1 or Tier 2.

## Invariants You Should Monitor
The code enforces and logs values consistent with these checks:
- `Mdot_hot_realized + Mdot_cold_realized = dm_removed / dt`
- `bh_mass_new = bh_mass_old + dm_removed`
- `BH_mass = BHFormationMass + BHAccretedMass`

If these drift materially, treat as a blocking issue.

## `[BHACCR]` Log Guide
Each BH pass emits one `[BHACCR]` line (when verbose is enabled).

Fields include:
- context: `step`, `level`, `z`, `bh_id`, `bh_mass`
- split/classification: `f_hot`, `f_cold`, `n_hot_cells`, `n_cold_cells`, `n_fallback_cells`
- averaged gas properties and velocities
- rates: `Mdot_hot_raw`, `Mdot_cold_raw`, `Mdot_total_raw`, `Mdot_Edd`, `f_Edd`
- cap and realization: `cap_active`, `frac_cap`, `frac_gas`, `dm_requested`, `dm_removed`
- realized rates: `Mdot_hot_realized`, `Mdot_cold_realized`
- bookkeeping/output: `bh_mass_new`, `removal_cells`, `n_sf_blocked_cells`
- momentum monitor: `acc_ignored_dv_kms`, `acc_ignored_p_frac`, `acc_momentum_warn`
- timing: `accretion_diag_wall_ms`, `accretion_wall_ms`

### `[BHACCR_WARN]` Lines
Expected warning types:
- under-resolved or oversized kernel warnings,
- high fallback-cell fraction,
- hot boost at cap (`alpha_boost`),
- gas-limited removal (`removal_gas_limited=1`),
- momentum-omission warning (threshold exceeded).

Warnings are diagnostics, but repeated warnings should be investigated.

## Metadata Updated on Particles
Main accretion metadata fields:
- `BHAccretedMass` (increases by `dm_removed`),
- `BHLastAccretionRedshift` (updated when `dm_removed > 0`),
- `BHLastEddingtonRatio` (updated each pass),
- `BHFormationMass` (seed-set reference mass),
- `BHReservoirMass` (reserved; unchanged in current phase).

## Star Formation Coordination
All active-zone cells that actually lose gas (`dm_i > 0`) are flagged in the transient SF mask.
Star formation then skips those cells during the same pass.

Practical expectation:
- `n_sf_blocked_cells` should match cells where `dm_i > 0` in removal.

## Cooling Backend Note (Built-in vs Grackle)
Hot/cold splitting depends on `ComputeCoolingTime`.
That means cooling backend choice matters:
- built-in Enzo cooling and Grackle can produce different cooling times,
- different cooling times can change hot/cold fractions and channel rates,
- this is expected physics-path sensitivity, not automatically a bug.

What should remain true in either backend:
- no NaN rates,
- no negative removal mass,
- invariant checks still hold,
- deterministic behavior within expected MPI floating-point tolerance.

## Minimal Recommended Parameter Block
Use explicit parameters in your `.enzo` file:

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

## Quick Validation Checklist
After enabling accretion, check logs for:
1. `[BHACCR]` lines appear for BH particles (except same-pass newborn seeds).
2. `dm_removed` is non-negative and finite.
3. `bh_mass_new` increases when `dm_removed > 0`.
4. `Mdot_hot_realized + Mdot_cold_realized` matches `dm_removed / dt`.
5. `BH_mass = BHFormationMass + BHAccretedMass` remains true.
6. `removal_cells` and `n_sf_blocked_cells` are sensible.

## Troubleshooting
### No `[BHACCR]` lines
- confirm `BHAccretionMethod = 1`,
- confirm BH particles exist and have valid IDs,
- check `BHAccretionVerbose` and stdout capture.

### `dm_removed = 0` repeatedly
- removal kernel may be gas-poor,
- `BHAccretionRemovalRadius` may be too small,
- local gas density may be very low.

### Frequent `removal_gas_limited=1`
- accretion demand exceeds available local gas,
- increase removal radius or reduce demand-driving conditions.

### Frequent momentum warnings
- flow is strongly asymmetric around BH,
- warning is diagnostic (velocity update is intentionally omitted in this phase),
- tune thresholds if warnings are too noisy for your campaign.

### Different results between np=1 and np>1
- small floating-point drift can appear from ghost-zone staleness in diagnostic quantities,
- verify differences are at expected small relative levels,
- investigate if differences are large or monotonic.

## Relationship to Seeding Guide
- Read `BH_SEEDING_README.md` for seed creation, ranking, exclusion, and seeding-time gas removal.
- Read this guide for post-seeding mass growth and active accretion behavior.
