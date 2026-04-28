# V0a: Publication Validation Blueprint for the Enzo Black-Hole Module

**Status:** DRAFT
**Date:** 2026-04-28
**Scope:** Validation design — no code, no simulations
**Prerequisite:** D0a committed (685c7a9c)
**Authors:** bkoh + Claude Opus 4.6

---

## 1. Executive Summary

The Enzo BH module implements seeding, two-channel Bondi accretion with
Eddington capping, density-peak repositioning, and reservoir-based thermal
feedback. Phase B and Phase C are merged. D0a provides an engineering
smoke/regression harness for Phase C feedback diagnostics.

**The module is engineering-smoke-tested but not publication-validated.** D0a
tests run on a 50³ uniform non-cosmological box with no self-gravity, no
cooling, 2 timesteps, and 2 kpc cells. No burst event has ever occurred in any
test. No gas limitation has been triggered. The Eddington cap has never been
active (f_Edd ~ 10⁻⁶ in all runs). These tests prove code paths execute; they
do not demonstrate that the subgrid model produces physically reasonable
galaxies or BH populations.

**This document (V0a) defines the validation requirements** before D0b-lite,
D0c, D2 decision, pilot runs, or production simulations. It specifies:

- what the code implements (source-grounded inventory)
- what diagnostics exist and what is missing
- what each component must prove for publication
- how to frame controlled tests honestly
- which references anchor the model
- what calibration and validation mean for this module
- decision protocols for D1 and D2
- go/no-go criteria for production compute

**Recommended project sequence:**

1. V0a: this document (now)
2. D0b-lite: burst, conservation-through-burst, restart, MPI tolerance (after V0a)
3. D0c: formula check, Eddington cap, gas limitation, repositioning ON/OFF, mismatch quantification (after D0b-lite, with V0a-designed ICs)
4. D2 decision (using D0c data)
5. Pilot simulation (after D2 decision, with V0a-declared science question)
6. Calibration (using pilot data)
7. Production (after go/no-go criteria met)

---

## 2. Module Inventory

All claims below are source-confirmed. File references use the form
`File.C:NNN` for line numbers.

### 2A. Seeding

**Source files:** `Grid_BHSeedHandler.C`, `star_maker_bh_seed.C`,
`ReadParameterFile.C:1096-1124` (parsing), `SetDefaultGlobalValues.C:747-770`
(defaults), `ReadParameterFile.C:2240-2250` (validation).

**Gate order** (star_maker_bh_seed.C:49-166):

| Order | Gate | Parameter | Default | Notes |
|-------|------|-----------|---------|-------|
| 1 | Finest level | BHSeedRequireFinestLevel | 1 (ON) | Phase 1 |
| 2 | Density threshold | BHSeedOverdensityThreshold | 1000.0 (code) | Must be > 0 |
| 3 | Local density peak | BHSeedRequireLocalPeak | 1 (ON) | 26-neighbor stencil |
| 4 | Temperature ceiling | BHSeedTemperatureThreshold | 1e4 K | |
| 5 | Metallicity ceiling | BHSeedMetallicityThreshold | 1e-4 (abs) | Or BHSeedMetallicityThresholdInSolar via Z_sun=0.02 |
| 6 | Convergent flow | BHSeedVelDivCrit | 1 (ON) | Stencil depends on HydroMethod |
| 7 | Thermal criterion | BHSeedThermalCrit | 0 (OFF) | t_cool < t_dyn |
| 8 | Self-boundedness | BHSeedSelfBoundCrit | 0 (OFF) | Virial parameter α < 1 |

Note: No Jeans-mass gate exists. The thermal criterion is t_cool vs t_dyn, not
a Jeans comparison. The self-boundedness uses a virial parameter
α = (div²+curl²)/(2Gρ), not a binding-energy calculation.

**Exclusion zone** (Grid_BHSeedHandler.C:304-414):

| Mode | Parameter | Default | Description |
|------|-----------|---------|-------------|
| 0 | BHSeedExclusionRadius | 100.0 kpc | Fixed physical kpc |
| 1 | (same) | — | Comoving kpc/h |
| 2 | BHSeedExclusionCells | 16 | Resolution-scaled (cells × dx) |

Default mode: BHSeedExclusionMode = 2 (global_data.h:1134). Some source
comments carry historical "Phase 3+: not yet active" annotations, but source
inspection shows the exclusion implementation is active: it uses a linked-cell
spatial hash (lines 198-276) with periodic boundary wrapping (lines 313-332).

**Seed mass:** Fixed per BH type, from BHSeedMass (Msun, default 1e5),
converted to code density units (Grid_BHSeedHandler.C:647).

**Enclosed-mass kernel:** BHSeedPatchRadius (default 3.0 cell widths),
BHSeedMinEnclosedMass (default 1e6 Msun, must >= BHSeedMass). Phase 2 active.

**Legacy cell-mass gate:** BHSeedLegacyCellMassGate (default 0 = shadow-only
diagnostic; 1 = hard reject if cell mass < seed mass).

**Active parameters with historical "Phase 3+" annotations:**
BHSeedExclusionMode, BHSeedExclusionCells, BHSeedMinCandidateSeparation,
BHSeedMaxPerPass, BHSeedRankingOrder, BHSeedDeterministicTiebreak. These
parameters may carry historical "Phase 3+ not yet active" annotations, but
source inspection shows they are active in the current implementation. Their
exact behavior must be documented and validated. Do not describe them as inert
unless future source inspection proves otherwise.

**Diagnostics:**

- `[BHSEED]`: per-level gate census with gate rejection counts, candidate
  counts, creation counts, seed totals, cache diagnostics, wall time
  (Grid_BHSeedHandler.C:1717-1744)
- `[BHSEED_SEED]`: per-seed creation details (position, channel, redshift,
  patch mass/metallicity/density, kernel completeness)
- `[BHSEED_WARN]`: large candidate counts, kernel evaluation truncation
- `[BHSEED_DEBUG]`: mass/momentum conservation at BHSeedVerbose >= 2

**HDF5/restart particle attributes** (typedefs.h, Grid_WriteGrid.C):

| Index (non-WINDS) | Label | Set when |
|--------------------|-------|----------|
| 4 | bhseed_channel | Seed creation |
| 5 | bhseed_redshift | Seed creation |
| 6 | bhseed_patch_mass | Seed creation |
| 7 | bhseed_patch_metallicity | Seed creation |
| 8 | bhseed_patch_density_peak | Seed creation |
| 9 | bhseed_kernel_complete | Seed creation |
| 10 | bhseed_host_dm_density | Seed creation |
| 11 | bhseed_accept_rank | Seed creation |

### 2B. Accretion

**Source files:** `Grid_BHAccretionHandler.C`,
`ReadParameterFile.C:1151-1168` (parsing), `SetDefaultGlobalValues.C:778-795`
(defaults), `ReadParameterFile.C:2275-2300` (validation).

**Hot/cold temperature split** (Grid_BHAccretionHandler.C:355-367):

```
if cooling available AND t_cool/t_dyn < 1.0:
    cold_cell = TRUE
else if cooling unavailable:
    cold_cell = (T < BHAccretionTSplitFloor)    # default 5e5 K
```

Note: the fallback path uses BHAccretionTSplitFloor, not a physics-based
criterion. The fraction of cells classified via fallback is tracked as
n_fallback_cells in the `[BHACCR]` diagnostic.

**Hot channel — boosted Bondi** (lines 444-458):

```
mdot_bondi_hot = 4π G² M_BH² ρ_hot / (c_s,hot² + v_rel,hot²)^(3/2)
n_H,hot = ρ_hot × rho_to_nh
α = clamp( (n_H,hot / BHAccretionNHStar)^BHAccretionBeta, 1, BHAccretionAlphaMax )
mdot_hot_raw = α × mdot_bondi_hot
```

Defaults: BHAccretionNHStar = 0.1 cm⁻³, BHAccretionBeta = 1.0,
BHAccretionAlphaMax = 10.0.

**Cold channel — AM-suppressed Bondi** (lines 464-478):

```
mdot_bondi_cold = 4π G² M_BH² ρ_cold / (c_s,cold² + v_rel,cold²)^(3/2)
f_AM = min(1, BHAccretionCVisc × (c_s,cold / v_rot,cold)³)     if v_rot > 0
f_AM = 1.0                                                       if v_rot = 0
mdot_cold_raw = f_AM × mdot_bondi_cold
```

Default: BHAccretionCVisc = 2π ≈ 6.283. v_rot is computed from specific
angular momentum in the kernel (lines 388-390, 427-428).

**Eddington cap** (lines 480-506):

```
mdot_edd_cgs = 4π G M_BH m_H / (ε_r σ_T c)
lambda_edd = mdot_total_raw / mdot_edd_code
if BHFeedbackEddingtonFactor == 1.0:
    frac_cap = min(1, mdot_edd_code / mdot_total_raw)     # exact Phase B path
else:
    frac_cap = min(1, EddFactor × mdot_edd_code / mdot_total_raw)
cap_active = (frac_cap < 1)
```

Default: BHFeedbackEddingtonFactor = 1.0, BHAccretionRadiativeEfficiency = 0.1.

Note: lambda_edd always uses the base physical Eddington rate, not the
factor-modified rate. The Eddington factor only modifies the cap.

**Three-tier rate system:**

| Tier | Variable | Line | Used for |
|------|----------|------|----------|
| 1 raw | mdot_total_raw | 480 | Diagnostic only |
| 2 capped | mdot_total_capped = mdot_actual | 510, 513 | Stored as BHACCR_LAST_MDOT_ACTUAL (line 817); read by feedback for L_feedback |
| 3 gas-realized | dm_removed_total | 628-641 | BH mass growth (line 781) |

**This is the source-confirmed requested-vs-realized mismatch:**

- BH mass grows from tier 3: `ParticleMass[p] = bh_mass_code + dm_removed_total` (line 781)
- Feedback energy derives from tier 2: `L_feedback = ε_r × mdot_actual × c²` (Grid_BHFeedbackHandler.C:283)
- When frac_gas < 1: feedback energy > mass-energy consumed

**Gas removal** (lines 525-696): Two modes:

| Mode | Parameter | Description |
|------|-----------|-------------|
| 0 | BHAccretionRemovalRadius (default 1 cell) | Spherical removal |
| 1 | — | Host-cell only |

Gas removal distributes dm proportionally to cell mass, with iterative
residual backfill (lines 658-696). Star-formation blocking applied to
removed cells (lines 750-754).

**Mass bookkeeping** (lines 781-837):

- Runtime invariant: |M_BH - (M_formation + M_accreted)| < 1e-8 × M_BH
  (ENZO_FAIL if violated, line 835-836)
- M_formation set once at seed creation (constant thereafter)
- M_accreted incremented by dm_removed_total each step (line 803)

**Diagnostics — `[BHACCR]`** (lines 923-962):

All three rate tiers, channel splits (f_hot, f_cold, n_hot_cells,
n_cold_cells, n_fallback_cells), kernel-averaged ρ/T/c_s/v for each channel,
alpha_boost, f_AM, cap_active, frac_cap, frac_gas, removal_gas_limited,
dm_requested, dm_removed, momentum diagnostics (acc_ignored_dv_kms,
acc_ignored_p_frac), lambda_edd, edd_factor.

**Parameters:**

| Parameter | Default | Units | Validation |
|-----------|---------|-------|------------|
| BHAccretionMethod | 1 | — | {0,1} |
| BHAccretionKernelRadius | 3.0 | physical kpc | > 0 |
| BHAccretionRemovalRadius | 1 | cell widths | — |
| BHAccretionRemovalMode | 0 | — | — |
| BHAccretionTSplitFloor | 5e5 | K | — |
| BHAccretionColdModel | 0 | — | Must be 0 in Phase B |
| BHAccretionCVisc | 6.283 | — | > 0 |
| BHAccretionNHStar | 0.1 | cm⁻³ | — |
| BHAccretionBeta | 1.0 | — | — |
| BHAccretionAlphaMax | 10.0 | — | — |
| BHAccretionRadiativeEfficiency | 0.1 | — | (0,1]; ε_r×ε_f ≤ 1.0 |
| BHAccretionRunEveryTimestep | 0 | — | {0,1} |
| BHAccretionIgnoredDVWarn | 1.0 | km/s | ≥ 0 |
| BHAccretionIgnoredPFracWarn | 0.01 | — | ≥ 0 |

**HDF5/restart particle attributes** (non-WINDS):

| Index | Label | Updated |
|-------|-------|---------|
| 12 | bhaccr_accreted_mass | Every accretion step (cumulative) |
| 13 | bhaccr_reservoir_mass | Legacy; reset if invalid |
| 14 | bhaccr_last_accretion_redshift | When dm_removed > 0 |
| 15 | bhaccr_last_eddington_ratio | Every accretion step |
| 16 | bh_formation_mass | Set once; constant |
| 17 | bhaccr_last_mdot_actual | Every accretion step |

### 2C. Repositioning / Dynamics

**Source files:** `Grid_BHRepositionHandler.C`,
`ReadParameterFile.C:2266-2273` (validation),
`SetDefaultGlobalValues.C:772-776` (defaults).

**Methods** (Grid_BHRepositionHandler.C:441-468):

| Method | Description | Default |
|--------|-------------|---------|
| 0 | Diagnostics only — no movement | Default |
| 1 | Rate-limited drift toward active-zone density peak | — |
| 2 | Teleport to density peak (debug) | — |

**Search target:** Density peak within BHRepositionSearchRadius (default
3.0 physical kpc), converted to code units with cosmology support (lines
78-111, 141).

**Peak selection:**

- Diagnostic peak: all cells including ghosts (lines 351-362)
- Active peak: active-zone cells only, used for actual repositioning
  (lines 364-375)
- Lexicographic tiebreak for determinism (lines 358-362, 371-374)
- Optional: potential-minimum tracking if BHRepositionDiagnosePotential=1
  (lines 377-397; diagnostic only, not used for repositioning target)

**Movement limits:**

- BHRepositionMaxDisplacement (default 0.5 cell widths): maximum per-step
  movement
- Position clamped to active grid zone with edge tolerance (lines 470-505)
- Newly-seeded BHs skip repositioning on creation step (line 257)

**Diagnostics — `[BHREPOS]`** (lines 512-534):

bh_pos, diag_peak_pos/density/offset, active_peak_pos/density/offset,
potential-minimum offset (if enabled), displacement_kpc, displacement_cells,
reposition_occurred, reposition_clamped, newly_seeded_skip, search_cells,
search_active_cells, wall time.

**Missing: dynamical friction.** No dynamical friction subgrid model exists in
this code. Repositioning to density peak is the only mechanism for keeping the
BH at the halo center. This is standard in many codes (EAGLE, Illustris,
Horizon-AGN) but is known to artificially enhance accretion by placing the BH
at ρ_max. Romulus (Tremmel+ 2015, 2017) and FIRE (Ma+ 2021) use alternatives.

**Parameters:**

| Parameter | Default | Units | Validation |
|-----------|---------|-------|------------|
| BHRepositionMethod | 0 | — | {0,1,2} |
| BHRepositionSearchRadius | 3.0 | physical kpc | > 0 |
| BHRepositionMaxDisplacement | 0.5 | cell widths | ≥ 0 |
| BHRepositionDiagnosePotential | 0 | — | {0,1} |

### 2D. Feedback

**Source files:** `Grid_BHFeedbackHandler.C`,
`ReadParameterFile.C:2302-2330` (validation),
`SetDefaultGlobalValues.C:798-806` (defaults).

**Feedback energy source** (Grid_BHFeedbackHandler.C:264-283):

```
mdot_actual_code = ParticleAttribute[BHACCR_LAST_MDOT_ACTUAL][p]
mdot_cgs = mdot_actual_code × mass_rate_to_cgs
L_feedback = ε_r × mdot_cgs × c²
```

This uses the tier-2 capped rate stored by the accretion handler. It does NOT
use the tier-3 gas-realized rate. See D2 decision protocol (Section 9).

**Mode switching** (line 271):

```
thermal_mode = (f_Edd > BHFeedbackModeThreshold)    # default threshold 0.01
feedback_mode = thermal_mode ? "THERMAL" : "KINETIC_INACTIVE"
```

Kinetic mode is scaffolded but not coupled — no energy deposition occurs for
KINETIC_INACTIVE BHs. If BHFeedbackMethod=2, a runtime warning is emitted
(lines 273-280).

**Reservoir accumulation** (lines 301-318):

```
if thermal_mode AND mdot_actual > 0:
    E_requested = ε_f × L_feedback × dt
    cumul_reservoir_in += E_requested                    # non-WINDS only
    reservoir = reservoir_before + E_requested
    if reservoir >= BHFeedbackMinEnergyBurst:
        E_burst_candidate = reservoir        # tentative release candidate
        burst_diag = 1
    else:
        reservoir_final = reservoir
```

Default BHFeedbackMinEnergyBurst = 1e50 erg.

**Thermal deposition** (lines 378-449): Mass-weighted across kernel cells
within BHFeedbackKernelRadius (default 1.0 physical kpc).
`dE_cell = E_burst × (m_cell / m_kernel_total)`. Energy update respects
dual-energy formalism and HydroMethod (PPM vs Zeus). Star-formation blocking
applied to heated cells.

If the burst threshold is reached but the feedback kernel contains no gas
(`kernel_gas_zero=1`), the burst energy is returned to the reservoir. This is
not a successful deposited burst, should not satisfy D0b-lite F1, and does
not imply `reservoir_final=0`.

**Conservation tracking** (lines 479-523, Grid_BHFeedbackHandler.C):

- cumul_reservoir_in: total E_requested accumulated (code energy units)
- cumul_reservoir_out: total E_burst released (code energy units)
- conservation_residual_cgs = (cumul_in - cumul_out - reservoir) converted to CGS; logged signed
- Logged but not enforced via ENZO_FAIL (diagnostic, not invariant)

**WINDS vs non-WINDS caveat:**

The cumulative reservoir tracking is gated by `#ifndef WINDS`
(Grid_BHFeedbackHandler.C:307). In a WINDS build:

- cumul_reservoir_in and cumul_reservoir_out are NOT incremented
- The particle attributes for these fields do not exist (typedefs.h:368-386)
- The HDF5 labels bh_cumul_reservoir_in/out are not written
  (Grid_WriteGrid.C:93-104 vs 106-119)

**Production build must be confirmed as non-WINDS if cumulative tracking and
the D0a T6 HDF5 label regression guard are to remain valid.**

**Diagnostics — `[BHFDBK]`** (lines 558-584):

feedback_mode, f_Edd, L_feedback, E_requested, reservoir_before/after/final,
burst_diag, burst_occurred, E_deposited, p_requested/deposited,
kernel_cells/active_cells/gas_msun, T_before/after_mean, dT_mean,
n_sf_blocked_feedback, deposit_rel_err, mass_weight_rel_err,
newly_seeded_skip, feedback_wall_ms, cumul_reservoir_in/out_cgs,
conservation_residual_cgs.

**Parameters:**

| Parameter | Default | Units | Validation |
|-----------|---------|-------|------------|
| BHFeedbackMethod | 0 | — | {0,1,2} |
| BHFeedbackModeThreshold | 0.01 | — | — |
| BHFeedbackKernelRadius | 1.0 | physical kpc | > 0 |
| BHFeedbackEddingtonFactor | 1.0 | — | ≥ 0 |
| BHFeedbackThermalEfficiency | 0.02 | — | [0,1]; ε_r×ε_f ≤ 1.0 |
| BHFeedbackMinEnergyBurst | 1e50 | erg | > 0 |
| BHFeedbackKernelMassWarnThreshold | 1e3 | Msun | > 0 |
| BHFeedbackKineticEfficiency | 0.1 | — | ≥ 0 (parsed, Phase C) |
| BHFeedbackWindVelocity | 1e4 | km/s | > 0 (parsed, Phase C) |
| BHFeedbackKineticGeometry | 0 | — | {0,1} (parsed, Phase C) |

Source default is `BHFeedbackMethod = 0`; feedback is disabled unless
explicitly enabled. D0a/D0b/D0c feedback tests and production parameter files
must explicitly set `BHFeedbackMethod = 1` when feedback is intended. No
roadmap item or go/no-go criterion should imply default feedback-on behavior.

### 2E. Restart / HDF5 State

**Particle attribute layout** (typedefs.h:364-408, Grid_WriteGrid.C:88-120):

Two compile-time layouts exist:

| Build | Standard attrs | BH attrs start | Cumul reservoir | Total BH attrs |
|-------|----------------|-----------------|-----------------|----------------|
| non-WINDS | 0-3 (4 attrs) | Index 4 | Indices 20-21 | 22 |
| WINDS | 0-6 (7 attrs, includes jet xyz) | Index 7 | NOT PRESENT | 23 |

**Non-WINDS HDF5 labels** (complete BH-relevant list from Grid_WriteGrid.C):

```
 4: bhseed_channel                 12: bhaccr_accreted_mass
 5: bhseed_redshift                13: bhaccr_reservoir_mass
 6: bhseed_patch_mass              14: bhaccr_last_accretion_redshift
 7: bhseed_patch_metallicity       15: bhaccr_last_eddington_ratio
 8: bhseed_patch_density_peak      16: bh_formation_mass
 9: bhseed_kernel_complete         17: bhaccr_last_mdot_actual
10: bhseed_host_dm_density         18: bhfdbk_energy_reservoir
11: bhseed_accept_rank             19: bhfdbk_last_feedback_redshift
                                   20: bh_cumul_reservoir_in
                                   21: bh_cumul_reservoir_out
```

**Restart-critical fields** (fields whose corruption would invalidate a
resumed simulation):

- ParticleMass (BH mass)
- ParticlePosition
- ParticleVelocity, if available and meaningful in the output format
- bh_formation_mass (seed mass, constant)
- bhaccr_accreted_mass (cumulative accretion)
- bhfdbk_energy_reservoir (thermal reservoir, determines next burst)
- bhaccr_last_mdot_actual (feedback energy source for next step)
- bhaccr_last_eddington_ratio (mode switching for next step)
- bh_cumul_reservoir_in / out (conservation tracking, non-WINDS only)

---

## 3. Missing Diagnostics Needed for Publication

The following information is not currently summarized by the code or its
diagnostic logs. All would need to be extracted via post-processing of
`[BHSEED]`, `[BHACCR]`, `[BHFDBK]`, `[BHREPOS]` log lines or HDF5 snapshots.

| Missing diagnostic | What it would provide | Needed for |
|--------------------|----------------------|------------|
| Per-BH lifetime mass history | M_BH(t) from birth to present | BH growth tracks figure |
| Per-snapshot BH-halo association | Which BH is in which halo | M_BH-M_star, M_BH-σ |
| Burst event catalog | When/where/how-much for each burst | Duty cycle measurement, D1 decision |
| Population Eddington-ratio distribution | f_Edd histogram at each snapshot | Comparison to Aird+ 2018 |
| Hot/cold channel fraction summary | Fraction of accretion via hot vs cold, per BH or per snapshot | Model characterization figure |
| Requested-vs-realized energy summary | Cumulative E(mdot_actual) vs E(dm_removed) per BH | D2 decision metric |
| Repositioning displacement history | Time series of displacement_kpc per BH | R4 analysis |
| Bolometric luminosity in HDF5 | L_bol per BH per snapshot, queryable without log parsing | AGN luminosity function |
| Gas-limitation frequency | Fraction of steps with frac_gas < 1, per BH | D2 relevance assessment |

None of these require code changes to the Enzo source — they are
post-processing tasks on existing log output and HDF5 snapshots. An analysis
script suite should be designed alongside the figure pipeline.

---

## 4. Component Validation Matrix

### 4.1 Seeding

**What must be proven:**

- Seeds form in physically plausible environments (overdense, cool,
  metal-poor gas) and not elsewhere
- Exclusion zones prevent spurious over-seeding, including near periodic
  boundaries
- Seed mass is correctly assigned in code units
- Results are stable (not necessarily convergent) with resolution

**Controlled tests needed:**

| Test | Description | Acceptance criterion | Phase |
|------|-------------|---------------------|-------|
| S1 | Exclusion zone: two qualifying cells within active exclusion radius/mode | Exactly 1 seed; documented behavior for BHSeedExclusionMode and BHSeedExclusionCells | D0c |
| S2 | Exclusion zone wrap: qualifying cell near periodic boundary, existing BH across boundary within exclusion distance | Correctly blocked | D0c |
| S3 | Gate sensitivity: enable each gate independently, document rejection fractions | No crash; fractions documented | D0c |
| S4 | Resolution: same physical setup at 2x resolution | Seed count stable ±1; location converges | D0c |
| S5 | Active candidate controls: MinCandidateSeparation, MaxPerPass, RankingOrder, DeterministicTiebreak | Exact behavior documented; deterministic ordering/tiebreak verified where enabled | D0c |

S3 and S4 require ICs with self-gravity and cooling (for thermal and
self-bound gates). S1 and S2 may be feasible on TS3_wrap variants.

**Comparison data to collect:**

- Habouzit+ 2017: seed occupation fraction vs halo mass across codes
  (methods anchor)
- Greene+ 2020: dwarf BH occupation fraction (observational constraint)

**Diagnostics required:** `[BHSEED]` gate statistics (exist). Per-halo
seed census requires halo finder (missing; post-processing).

**Publication risks:**

- Low-mass BH outcomes are dominated by seed mass and exclusion radius
- Paper must separate seeding-dominated (M_BH < few × M_seed) from
  accretion-dominated BHs
- Parameters with historical "Phase 3+" annotations
  (BHSeedExclusionMode, BHSeedExclusionCells, BHSeedMinCandidateSeparation,
  BHSeedMaxPerPass, BHSeedRankingOrder, BHSeedDeterministicTiebreak) are
  active in current source and must have documented, validated behavior

**Before-production requirement:** S1 (basic exclusion) is
REQUIRED_BEFORE_PRODUCTION. S3/S4 are REQUIRED_BEFORE_PUBLICATION. Full
comparison with Habouzit+ 2017 is OPTIONAL_APPENDIX or FUTURE_WORK depending
on paper scope.

### 4.2 Accretion

**What must be proven:**

- Bondi formula is correctly evaluated (units, constants, kernel averaging)
- Hot/cold channel split works as documented
- Eddington cap activates at the correct luminosity
- Angular momentum suppression reduces cold accretion by expected factor
- Alpha boost enhances hot accretion in dense environments
- Mass bookkeeping invariant holds over long runs
- Gas-limited path (frac_gas < 1) correctly limits mass growth

**Controlled tests needed:**

| Test | Description | Acceptance criterion | Phase |
|------|-------------|---------------------|-------|
| A1 | Uniform hot medium, SelfGravity=1. Compare code mdot to Bondi formula evaluated on kernel-averaged ρ, c_s. | Agreement to floating-point precision for uniform medium | D0c |
| A2 | High-density medium forcing mdot_raw > mdot_Edd. | cap_active=1; Mdot_actual = EddFactor × M_dot_Edd ± roundoff | D0c |
| A5 | Low-density or high-M_BH forcing frac_gas < 1. | frac_gas < 1; removal_gas_limited=1; M_BH grows by dm_removed only | D0c |
| A6 | Long run (100+ steps). Verify mass bookkeeping invariant. | Runtime ENZO_FAIL never triggers; post-hoc check confirms | D0c |

**Correct framing for A1:** See Section 5A. This is a formula/unit-regime
check, not hydrodynamic Bondi-flow validation.

**Comparison data to collect:**

- Booth & Schaye 2009: alpha boost motivation and calibration
- Rosas-Guevara+ 2015: EAGLE angular momentum implementation (closest
  comparison to this module's cold-channel f_AM)

**Diagnostics required:** All `[BHACCR]` fields (exist).

**Publication risks:**

- Boosted Bondi is resolution-dependent by construction; paper must state
  target resolution and acknowledge calibration at that resolution
- f_AM depends on kernel-scale velocity structure (resolution-dependent)
- BHAccretionTSplitFloor fallback (5e5 K) is ad hoc; if many cells use
  fallback path, the hot/cold split is not physically grounded

**Before-production requirement:** A1 and A2 are REQUIRED_BEFORE_PRODUCTION.
A5 is REQUIRED_BEFORE_PRODUCTION (also feeds D2 decision). A6 is
REQUIRED_BEFORE_PUBLICATION.

### 4.3 Repositioning / Dynamics

**What must be proven:**

- Repositioning to density peak is deterministic and reproducible
- Displacement bounded by MaxDisplacement
- The enhancement to BH growth from repositioning is quantified

**Controlled tests needed:**

| Test | Description | Acceptance criterion | Phase |
|------|-------------|---------------------|-------|
| R1 | BH at density peak in isolated halo; verify no movement | displacement_kpc ≈ 0 at all steps | D0c |
| R2 | Offset BH; verify drift converges toward peak | Monotonic convergence, bounded by MaxDisplacement per step | D0c |
| R4 | Repositioning ON vs OFF; same IC; compare BH mass growth | Enhancement factor documented | D0c |

R1, R2, R4 all require non-uniform density (isolated halo IC). Cannot be
done on TS3_wrap (uniform density; repositioning has no effect).

**Comparison data to collect:**

- Tremmel+ 2015, 2017: dynamical friction alternative (key comparison)
- Ma+ 2021: FIRE without repositioning
- Sijacki+ 2007: repositioning effects

**Diagnostics required:** `[BHREPOS]` fields (exist).

**Publication risks:**

- Repositioning is the most contentious modeling choice. Referees will
  question it.
- If R4 shows >2x enhancement, repositioning dominates BH growth and the
  paper is effectively about repositioning, not accretion physics
- Paper must discuss this honestly

**Before-production requirement:** R4 is REQUIRED_BEFORE_PRODUCTION (informs
whether repositioning settings need adjustment). R1/R2 are
REQUIRED_BEFORE_PUBLICATION.

### 4.4 Feedback

**What must be proven:**

- Burst events occur and deposit measurable thermal energy
- Reservoir accumulation and release conserve energy
- Conservation tracking holds through burst events
- The requested-vs-realized mismatch is quantified

**Controlled tests needed:**

| Test | Description | Acceptance criterion | Phase |
|------|-------------|---------------------|-------|
| F1 | Low BHFeedbackMinEnergyBurst; verify gas-present burst fires | kernel_gas_zero=0; burst_diag=1; burst_occurred=1; E_deposited > 0; reservoir_final=0 or tiny tolerance; cumul_reservoir_out_cgs increases; abs(conservation_residual_cgs) passes criterion | D0b-lite |
| F3 | Conservation audit through burst: cumul_in - cumul_out = reservoir at every step | If cumul_reservoir_in_cgs > 0: abs(conservation_residual_cgs) <= 1e-5 × cumul_reservoir_in_cgs. If cumul_reservoir_in_cgs == 0: abs(conservation_residual_cgs) <= tiny absolute tolerance | D0b-lite |
| F6 | Gas-limited run; compute mismatch metric (see Section 9) | Metric value documented | D0c |

**Comparison data to collect:**

- Booth & Schaye 2009: reservoir model (closest comparison)
- Dalla Vecchia & Schaye 2012: minimum dT requirements for overcooling
  prevention

**Diagnostics required:** All `[BHFDBK]` fields (exist). Burst event catalog
requires post-processing (missing).

**Publication risks:**

- No burst has ever occurred in any test run. The burst code path is
  untested under load. D0b-lite F1 addresses this.
- `kernel_gas_zero=1` is a failed F1 outcome, not a success-like burst. It
  may be kept as a separate negative or edge-case test outside D0b-lite.
- The mdot_actual vs dm_removed mismatch is the biggest known accounting
  issue. D0c F6 + D2 decision protocol addresses this.
- WINDS build disables conservation tracking. Must confirm build.

**Before-production requirement:** F1 and F3 are REQUIRED_BEFORE_PRODUCTION.
F6 is REQUIRED_BEFORE_PRODUCTION (feeds D2 decision).

### 4.5 End-to-End Observables

**What must be proven:**

- BH mass–stellar mass relation has reasonable normalization
- Stellar mass function is not destroyed by feedback
- BH feedback quenches massive galaxies

**Tests needed:** These require cosmological simulations (pilot or
production). No controlled tests possible on isolated ICs for emergent
population statistics.

| Test | Description | Acceptance criterion | Phase |
|------|-------------|---------------------|-------|
| E2 | Cosmological box to z=0 | M_BH-M_star within 0.5 dex of observations at z=0 | Pilot/prod |
| E3 | Resolution study | Key relations stable or non-convergence documented | Production |

**Comparison data to collect:**

- Kormendy & Ho 2013; Reines & Volonteri 2015 (M_BH-M_star)
- Schaye+ 2015; Pillepich+ 2018; Davé+ 2019 (EAGLE/TNG/SIMBA comparison)
- Habouzit+ 2021 (cross-simulation BH comparison)

**Diagnostics required:** Halo finder + BH-halo association + galaxy
property extraction. None of this exists in the Enzo BH module; requires
external analysis pipeline.

**Observational comparison pipeline:** This is a required pre-production
infrastructure item, not part of D0c. It includes halo finder setup,
BH-halo association, galaxy property extraction, M_BH-M_star extraction,
stellar mass/SFR catalogs, observational comparison tables, and
figure-generation scripts. D0c may define diagnostic requirements consumed by
this pipeline, but the pipeline itself is a separate pre-production
infrastructure task.

**Publication risks:**

- Without a target science question, end-to-end validation is unbounded
- Paper must decide: methods paper vs science paper (see Section 13)

**Before-production requirement:** Science question and calibration targets
must be declared (V0a decision; see Section 8). The observational comparison
pipeline must be prototyped as pre-production infrastructure outside D0c.
Both are REQUIRED_BEFORE_PRODUCTION.

---

## 5. Correct Test Framing

### 5A. Bondi Test

A uniform static-medium test with SelfGravity=1 validates that the code
correctly evaluates:

```
M_dot = 4π G² M² ρ / (c_s² + v²)^(3/2)
```

using kernel-averaged ρ, c_s, v, with correct unit conversions (CGS
constants, code-unit scaling, cosmology factors).

For a truly uniform medium (no gradients, no bulk velocity), the kernel
average equals the analytic input exactly. The acceptance criterion is:
code M_dot agrees with the analytic formula to floating-point precision
(~1e-12 relative). **This is a formula and unit-regime check.**

This test should not be described as "Bondi-flow validation" or "within 10%
after 3+ Bondi times." Those phrases imply a resolved steady-state inflow
solution, which requires resolving the Bondi radius (r_B = GM/c_s²). For a
10⁵ M_sun BH in 10⁷ K gas, r_B ~ 0.3 pc — unresolved by 2-3 orders of
magnitude at any production resolution. The code does not resolve the
accretion flow; it evaluates the Bondi formula at the grid scale.

A future non-uniform variant (density gradient across the kernel) would test
kernel-averaging accuracy. That test's tolerance would reflect kernel
averaging, not a physics discrepancy.

### 5B. MPI Test

np=1 vs np=4 comparison must be tolerance-based, not exact-equality.

MPI domain decomposition changes the order of floating-point reductions.
Kernel-averaged quantities (ρ, c_s, v) depend on which cells each rank
contributes, and finite-precision sum ordering.

Acceptance criteria for fields included in the hard comparison (see Section
5E):

- Single-step (step 2 of TS3_wrap): all diagnostic fields agree to ~1e-10
  relative tolerance (roundoff-level for one reduction pass)
- Multi-step (if tested): tolerance grows with step count. Exact criteria
  should be determined empirically during D0b-lite implementation. A
  practical starting point: ~1e-6 relative after 10 steps in a uniform
  medium.

### 5C. D0b-lite F1 Burst Configuration

D0b-lite F1 base: TS3_wrap copied into the output directory, with feedback
explicitly enabled. Base configuration:

```
BHAccretionMethod = 1
BHAccretionRunEveryTimestep = 1
BHAccretionVerbose = 1
BHFeedbackMethod = 1
BHFeedbackModeThreshold = 1e-10
BHFeedbackVerbose = 1
```

Burst override:

```
BHFeedbackMinEnergyBurst = 1e40   # erg
```

The implementation must verify from logs that `E_requested` exceeds
`BHFeedbackMinEnergyBurst`. If no gas-present burst occurs, F1 fails with an
actionable message rather than silently retuning the threshold or accepting a
no-gas branch.

Required gas-present burst acceptance:

- `kernel_gas_zero = 0`
- `burst_diag = 1`
- `burst_occurred = 1`
- `E_deposited > 0`
- `reservoir_final = 0` or within tiny numerical tolerance
- `cumul_reservoir_out_cgs` increases
- `abs(conservation_residual_cgs)` passes the absolute residual criterion

`kernel_gas_zero=1` is a failed F1 result. It may be tested separately as a
negative or edge-case outside D0b-lite.

### 5D. D0b-lite T7 Restart Contract

T7 compares continuous vs restart output at the same physical time, matched
by BH identity.

Fields to compare:

- ParticleMass
- ParticlePosition
- ParticleVelocity, if available and meaningful
- bhseed_channel
- bhseed_redshift
- bhseed_patch_mass
- bhseed_patch_metallicity
- bhseed_patch_density_peak
- bhseed_kernel_complete
- bhseed_host_dm_density
- bhseed_accept_rank
- bhaccr_accreted_mass
- bhaccr_reservoir_mass
- bhaccr_last_accretion_redshift
- bhaccr_last_eddington_ratio
- bh_formation_mass
- bhaccr_last_mdot_actual
- bhfdbk_energy_reservoir
- bhfdbk_last_feedback_redshift
- bh_cumul_reservoir_in
- bh_cumul_reservoir_out

Matching and tolerances:

- Match by BH particle ID if available; otherwise match by nearest position
  plus mass/formation-mass consistency.
- Categorical or integer-like fields compare exactly.
- Floating HDF5/restart fields use `rtol = 1e-6` unless source precision
  supports a tighter tolerance.
- Expected-zero energy fields use a tiny absolute tolerance.

### 5E. D0b-lite T10 MPI Diagnostic Contract

T10 compares np=1 vs np=4 diagnostics after normalizing rank-local ordering.
Matching keys:

- diagnostic tag
- step
- BH id
- feedback mode, where applicable

Compare BHACCR fields:

- lambda_edd
- edd_factor
- mdot_total_raw
- mdot_actual or mdot_total_capped
- mdot_edd
- frac_cap
- frac_gas
- f_hot
- f_cold
- alpha_boost
- f_AM
- dm_requested
- dm_removed

Compare BHFDBK fields:

- f_Edd
- E_requested
- reservoir_before
- reservoir_after
- reservoir_final
- E_deposited
- cumul_reservoir_in_cgs
- cumul_reservoir_out_cgs
- conservation_residual_cgs

Compare exactly/categorically:

- cap_active
- removal_gas_limited
- feedback_mode
- burst_diag
- burst_occurred
- kernel_gas_zero

Exclude from hard comparison:

- wall-clock timing fields
- absolute output paths
- rank-local ordering artifacts
- raw log ordering
- kernel/cell-count fields unless D0b-lite implementation confirms they are
  deterministic

Tolerances:

- Short TS3_wrap-style run: `rtol = 1e-10` for core scalar diagnostics.
- Expected-zero values: tiny absolute tolerance.
- Multi-step run: `rtol = 1e-6` unless D0b-lite empirically justifies tighter or
  looser tolerance.

### 5F. D2 Threshold

V0a defines the mismatch metric:

```
M = |Σ(ε_f × ε_r × mdot_actual × c² × dt) − Σ(ε_f × ε_r × (dm_removed/dt) × c² × dt)|
    / Σ(ε_f × ε_r × mdot_actual × c² × dt)
```

summed over all BHs and all timesteps in a representative run where gas
limitation occurs (frac_gas < 1 for a meaningful fraction of steps).

V0a does **not** set a hard implementation threshold. The D2 decision
compares M against other dominant modeling uncertainties (see Section 9).

### 5G. Pilot Run

Do not choose box size or resolution until the paper's target science
question is declared (see Section 13, Open Question 3). Box size depends on
the science target:

- Early-universe BH seeding → small box (10-25 Mpc/h), high resolution
- AGN feedback and quenching → larger box (50-100 Mpc/h)
- BH mass function statistics → volume matters
- Zoom of a single galaxy → different IC entirely

---

## 6. Reference Matrix

Priority: MUST_READ (~10 papers to read before D0c design), SHOULD_READ
(~15 papers for context), OPTIONAL (background).

### Seeding

| Reference | Category | Priority | Placement | Why |
|-----------|----------|----------|-----------|-----|
| Bromm & Loeb 2003 | methods anchor | SHOULD_READ | methods | Direct collapse seed physics |
| Begelman+ 2006 | methods anchor | SHOULD_READ | methods | Heavy seed channel |
| Habouzit+ 2017 | comparison sim | MUST_READ | validation | Cross-code seeding comparison; closest analog to our comparison needs |
| Bhowmick+ 2022 (BRAHMA) | comparison sim | SHOULD_READ | validation | Modern seeding-focused simulation |
| Ni+ 2022 (ASTRID) | comparison sim | SHOULD_READ | validation | Large-volume seeding study |
| Reines & Volonteri 2015 | obs constraint | SHOULD_READ | validation | Dwarf BH occupation; low-mass BH anchor |
| Greene+ 2020 | obs constraint | SHOULD_READ | validation | Low-mass BH demographics review |

### Accretion

| Reference | Category | Priority | Placement | Why |
|-----------|----------|----------|-----------|-----|
| Bondi 1952 | methods anchor | MUST_READ | methods | Foundational formula |
| Booth & Schaye 2009 | methods anchor | MUST_READ | methods | Alpha boost; reservoir model; closest comparison |
| Rosas-Guevara+ 2015 | methods anchor | MUST_READ | methods | EAGLE angular momentum suppression; closest f_AM comparison |
| Angles-Alcazar+ 2017 | comparison sim | SHOULD_READ | methods | Torque-limited accretion; alternative to boosted Bondi |
| Hopkins & Quataert 2011 | methods anchor | SHOULD_READ | methods | Gravitational torque theory |
| Kauffmann & Heckman 2009 | obs constraint | SHOULD_READ | validation | Eddington ratio distributions |
| Aird+ 2018 | obs constraint | SHOULD_READ | validation | Specific accretion rate distributions |

### Repositioning / Dynamics

| Reference | Category | Priority | Placement | Why |
|-----------|----------|----------|-----------|-----|
| Tremmel+ 2015 | methods anchor | MUST_READ | methods/discussion | Dynamical friction alternative; key comparison point |
| Tremmel+ 2017 (Romulus) | comparison sim | MUST_READ | discussion | Full code without repositioning |
| Sijacki+ 2007 | methods anchor | SHOULD_READ | methods | Repositioning effects quantified |
| Ma+ 2021 (FIRE) | comparison sim | SHOULD_READ | discussion | No repositioning; BH dynamics |
| Pfister+ 2019 | methods anchor | OPTIONAL | discussion | DF in cosmological sims |

### Feedback

| Reference | Category | Priority | Placement | Why |
|-----------|----------|----------|-----------|-----|
| Booth & Schaye 2009 | methods anchor | (already listed) | methods | Reservoir thermal model |
| Dalla Vecchia & Schaye 2012 | methods anchor | MUST_READ | methods | Minimum dT for overcooling prevention; stochastic injection |
| Weinberger+ 2017 (TNG) | comparison sim | SHOULD_READ | discussion | Dual-mode feedback (thermal + kinetic); contrast with our thermal-only |
| Davé+ 2019 (SIMBA) | comparison sim | SHOULD_READ | discussion | Bipolar jet feedback; alternative approach |

### Comparison Simulations

| Reference | Category | Priority | Placement | Why |
|-----------|----------|----------|-----------|-----|
| Springel+ 2005 | comparison sim | MUST_READ | methods | Foundational thermal feedback; original repositioning |
| Schaye+ 2015 (EAGLE) | comparison sim | MUST_READ | validation | Stochastic thermal + AM suppression; most similar to this module |
| Pillepich+ 2018 (TNG) | comparison sim | SHOULD_READ | validation | Major comparison target for BH populations |
| Vogelsberger+ 2014 (Illustris) | comparison sim | OPTIONAL | background | TNG predecessor |
| Dubois+ 2014 (Horizon-AGN) | comparison sim | OPTIONAL | background | Dual-mode with jet axis |
| Habouzit+ 2021 | comparison sim | MUST_READ | validation | Multi-code BH comparison at fixed resolution |

### Observational Constraints

| Reference | Category | Priority | Placement | Why |
|-----------|----------|----------|-----------|-----|
| Kormendy & Ho 2013 | obs constraint | MUST_READ | validation | M_BH-M_star and M_BH-σ (primary calibration candidate) |
| McConnell & Ma 2013 | obs constraint | SHOULD_READ | validation | M_BH-σ |
| Shankar+ 2009 | obs constraint | SHOULD_READ | validation | BH mass function |
| Hopkins+ 2007 | obs constraint | SHOULD_READ | validation | AGN luminosity function |
| Wetzel+ 2012 | obs constraint | SHOULD_READ | validation | Quenched fractions |
| Baldry+ 2012 | obs constraint | SHOULD_READ | validation | Stellar mass function |
| Donnari+ 2021 | obs constraint | OPTIONAL | validation | Quenched fractions (TNG comparison) |
| Sun+ 2009 | obs constraint | OPTIONAL | appendix | Hot gas fractions in groups |

**Total: ~30 references. MUST_READ: ~10. SHOULD_READ: ~15. OPTIONAL: ~5.**

---

## 7. Paper Figure/Table Skeleton

### Main Text Figures

| # | Purpose | Data source | Phase | Calibration? |
|---|---------|-------------|-------|--------------|
| 1 | Module schematic: seeding → accretion → repositioning → feedback flow diagram | Diagram (no simulation) | V0a/paper | N/A |
| 2 | Bondi formula check: code mdot vs analytic, demonstrating unit correctness | D0c A1 run | D0c | No (validation) |
| 3 | Eddington cap activation: mdot vs time showing cap engage | D0c A2 run | D0c | No (validation) |
| 4 | Reservoir sawtooth: reservoir energy vs time showing accumulation/burst cycles | D0c or pilot run with bursts | D0c or pilot | No (characterization) |
| 5 | BH mass growth tracks: representative BHs showing M_BH(t) from pilot/production | Pilot/production | Pilot+ | No (prediction) |
| 6 | M_BH vs M_star at z=0 with observational data | Production + Kormendy & Ho 2013 | Production | Yes (calibration target, if chosen) |
| 7 | Eddington ratio distribution vs Aird+ 2018 | Production | Production | No (validation) |
| 8 | Stellar mass function at z=0 vs Baldry+ 2012 | Production | Production | Possibly (calibration target) |
| 9 | Quenched fraction vs M_star | Production + Wetzel+ 2012 | Production | No (validation) |
| 10 | Hot/cold accretion channel fractions vs z or M_BH | Production | Production | No (diagnostic prediction) |

### Main Text Tables

| # | Purpose | Data source | Phase |
|---|---------|-------------|-------|
| T1 | Parameter table: all BH module parameters, defaults, values used | Source code + production config | Paper |
| T2 | Controlled test summary: A1, A2, F1, etc. with pass/fail | D0b-lite + D0c results | D0c |

### Appendix Figures/Tables

| # | Purpose | Data source | Phase |
|---|---------|-------------|-------|
| A1 | MPI reproducibility: np=1 vs np=4 diagnostic comparison | D0b-lite T10 | D0b-lite |
| A2 | Restart preservation: attribute comparison pre/post restart | D0b-lite T7 | D0b-lite |
| A3 | Conservation audit: signed cumul_in - cumul_out - reservoir residual and abs residual criterion vs time | D0b-lite F3 + D0c | D0b-lite+ |
| A4 | D2 mismatch: E(mdot_actual) vs E(dm_removed) over representative run | D0c F6 | D0c |
| A5 | Repositioning enhancement: M_BH growth ON vs OFF | D0c R4 | D0c |
| A6 | Seeding gate rejection breakdown | D0c S3 | D0c |
| A7 | Resolution sensitivity for key relations | Production E3 | Production |
| AT1 | Full parameter validation matrix (all gates, all checks) | Source code | Paper |

This skeleton is provisional. It will be revised after pilot runs reveal
which figures are informative and which are redundant.

---

## 8. Calibration-vs-Validation Declaration

### Proposed policy

**Free parameters available for tuning:**

| Parameter | Default | What it controls | Tuning range |
|-----------|---------|------------------|--------------|
| BHFeedbackThermalEfficiency (ε_f) | 0.02 | Feedback energy per unit accreted mass | 0.005-0.1 |
| BHFeedbackMinEnergyBurst | 1e50 erg | Burst frequency / duty cycle | 1e48-1e52 |
| BHSeedMass | 1e5 Msun | Seed mass | 1e3-1e6 |
| BHSeedExclusionRadius | 100 kpc | Seed density | 30-300 |
| BHAccretionAlphaMax | 10 | Hot boost ceiling | 1-100 |

**Parameters fixed by physics (not tuned):**

| Parameter | Value | Justification |
|-----------|-------|---------------|
| BHAccretionRadiativeEfficiency (ε_r) | 0.1 | Standard thin-disk value |
| BHAccretionCVisc | 2π | Rosas-Guevara+ 2015 / Booth & Schaye 2009 |
| BHAccretionBeta | 1.0 | Booth & Schaye 2009 |
| BHAccretionNHStar | 0.1 cm⁻³ | Star-formation threshold |
| BHFeedbackEddingtonFactor | 1.0 | Physical Eddington limit |
| BHFeedbackModeThreshold | 0.01 | Conventional high/low state boundary |

**Calibration targets (at most 2):**

1. M_BH-M_star normalization at z=0 (via ε_f and/or BHSeedMass)
2. Stellar mass function at z=0 (cross-check with galaxy formation params)

**Validation targets (predicted, not tuned):**

- BH mass function shape
- Eddington ratio distribution
- Quenched fraction vs M_star
- Hot/cold accretion channel fractions (unique to this model)
- Burst frequency / duty cycle (measured, not tuned; informs D1 decision)

**Diagnostic predictions (unique model outputs):**

- Angular momentum suppression factor f_AM distribution
- Alpha boost distribution
- Gas-limitation frequency (frac_gas < 1 incidence)
- Reservoir burst energy distribution

**Rule: no more than 2 calibration targets.** M_BH-M_star and M_BH-σ are
not independent (both driven by central potential). Using both as calibration
targets double-counts the same physics. Choose one.

---

## 9. D2 Decision Protocol

### The issue

`Grid_BHAccretionHandler.C:513` defines `mdot_actual = mdot_total_capped`
(tier 2, Eddington-capped). This is stored at line 817 and consumed by
`Grid_BHFeedbackHandler.C:282-283` as `L_feedback = ε_r × mdot_cgs × c²`.
BH mass grows from `dm_removed_total` (tier 3, gas-realized, line 781).

When frac_gas < 1, the BH receives feedback energy for mass it did not
actually consume. This is an energy accounting inconsistency.

### Mismatch metric

```
M = |Σ_i Σ_t [ ε_f ε_r mdot_actual(i,t) c² dt  −  ε_f ε_r (dm_removed(i,t)/dt) c² dt ]|
    / Σ_i Σ_t [ ε_f ε_r mdot_actual(i,t) c² dt ]
```

where i indexes BHs, t indexes timesteps, and the sum is over all BHs and
all steps in a representative run. This simplifies to:

```
M = Σ_i Σ_t [ mdot_actual(i,t) − dm_removed(i,t)/dt ] × dt
    / Σ_i Σ_t [ mdot_actual(i,t) × dt ]
```

since ε_f, ε_r, c cancel. Both mdot_actual and dm_removed are available
in the `[BHACCR]` log line (fields Mdot_actual and dm_removed).

M = 0 when frac_gas = 1 for all steps. M > 0 only when gas limitation
occurs.

### Data required

- D0c test F6: a run where frac_gas < 1 for a meaningful fraction of
  steps. This requires an IC with either low gas density or high BH mass
  relative to the local gas supply.
- Pilot run: M should also be measured in the pilot cosmological run to
  assess frequency under realistic conditions.

### Comparison uncertainties

The D2 decision compares M against:

| Uncertainty | Symbol | How to measure | Phase |
|-------------|--------|----------------|-------|
| Repositioning sensitivity | R_reposition | Total feedback energy change when BHRepositionMethod=0 vs 1 (R4 test) | D0c |
| Kernel sensitivity | R_kernel | Total feedback energy change at 2× kernel radius | D0c or pilot |
| Resolution sensitivity | R_resolution | Total feedback energy change at 2× resolution | Pilot or production |
| Calibration uncertainty | R_calibration | Range of ε_f producing acceptable M_BH-M_star | Calibration (P1) |
| Stochastic variability | R_stochastic | Variance across different random seeds or ICs | Pilot |

### Decision rule

- If M ≪ min(R_reposition, R_kernel, R_resolution): mismatch is
  subdominant. **Document and defer D2.** The paper states: "feedback
  energy is computed from the Eddington-capped rate; in gas-limited
  episodes, this overestimates feedback relative to actual mass growth
  by [M]. This is smaller than other modeling uncertainties
  (repositioning: [R_reposition], kernel averaging: [R_kernel])."

- If M is comparable to or larger than the smallest R_*: mismatch is a
  leading-order uncertainty. **Implement D2** (tie feedback energy to
  dm_removed-based rate).

- If M is large but gas limitation is rare (<1% of BH-steps): mismatch
  is large per-event but negligible in aggregate. **Document with
  per-event analysis.** No D2 implementation needed.

### Possible D2 implementation (if needed)

Change Grid_BHAccretionHandler.C to store dm_removed_total/dt as the
rate for feedback instead of mdot_actual. Requires:

- New or repurposed particle attribute for gas-realized rate
- Grid_BHFeedbackHandler.C reads the gas-realized rate instead of
  BHACCR_LAST_MDOT_ACTUAL
- Full D0a + D0b-lite + D0c regression

---

## 10. D1 Decision Protocol

### What D1 is

D1 (duty-cycle gating) would modulate feedback activation frequency to
match observed AGN duty cycles, beyond what the reservoir burst model
already provides.

### Why D1 is deferred

The reservoir model already produces duty-cycle-like behavior: energy
accumulates until the burst threshold is reached, then releases. The
natural burst frequency is determined by:

```
burst_interval ≈ BHFeedbackMinEnergyBurst / (ε_f × ε_r × mdot_actual × c²)
```

This has never been measured in a realistic run because no burst has
occurred in any test (burst_occurred=0 in all D0a logs). Adding D1
before measuring the natural duty cycle is premature optimization.

### Required measurements before D1 decision

1. D0b-lite F1: confirm bursts occur at all (engineering prerequisite)
2. Pilot run: measure natural burst interval distribution across the BH
   population
3. Compare natural burst interval to observed AGN duty cycles at
   comparable Eddington ratios (Kauffmann & Heckman 2009)

### Evidence that would justify D1

- Natural burst frequency is >10× higher or lower than observed duty
  cycles at the same f_Edd range
- AND the paper makes a specific claim about AGN duty cycles or
  variability

### D1 consequences if implemented

- Adds at least one new free parameter (duty-cycle modulation timescale
  or probability)
- Requires recalibration of ε_f and BHFeedbackMinEnergyBurst
- Changes the model's character enough that all controlled tests must
  be rerun

### Current recommendation

Do not implement D1. Measure natural burst frequency in pilot runs.
Revisit only if the evidence threshold above is met.

---

## 11. Production Go/No-Go Checklist

Production simulations should not begin until all items below are
satisfied.

### Engineering prerequisites

- [x] D0a: Phase C feedback smoke tests pass (685c7a9c)
- [ ] V0a: this document reviewed and approved
- [ ] D0b-lite F1: gas-present burst event exercised (kernel_gas_zero=0,
      burst_occurred=1)
- [ ] D0b-lite F3: conservation audit passes through burst using
      abs(conservation_residual_cgs)
- [ ] D0b-lite T7: restart preserves all BH particle attributes per
      Section 5D
- [ ] D0b-lite T10: np=1 vs np=4 agree within Section 5E tolerances

### Controlled physics prerequisites

- [ ] D0c A1: Bondi formula check passes (see Section 5A framing)
- [ ] D0c A2: Eddington cap activates correctly
- [ ] D0c A5: gas-limited accretion path works (frac_gas < 1)
- [ ] D0c F6: D2 mismatch metric M measured
- [ ] D0c R4: repositioning enhancement factor documented

### Decision prerequisites

- [ ] D2 decision made (implement, defer, or document-as-limitation)
- [ ] D1 decision made (almost certainly: defer to post-pilot)
- [ ] Science question declared (see Section 13)
- [ ] Calibration-vs-validation targets approved (Section 8)
- [ ] WINDS vs non-WINDS build confirmed for production

### Infrastructure prerequisites

- [ ] Reference matrix: at least MUST_READ papers (~10) read and annotated
- [ ] Observational comparison pipeline prototyped outside D0c: halo finder
      setup, BH-halo association, galaxy property extraction,
      M_BH-M_star extraction, stellar mass/SFR catalogs, observational
      comparison tables, and figure-generation scripts
- [ ] BH diagnostic log parser validated against D0a/D0b outputs
- [ ] Figure pipeline prototyped for at least Figures 2-4 (controlled tests)

### Pilot prerequisites

- [ ] Pilot run completes to target redshift without crash
- [ ] >50 BHs form; no runaway masses (>10¹¹ Msun)
- [ ] Calibration target(s) within ~0.5 dex of observations (pre-tuning)
- [ ] No pathological behavior (universal quenching, empty BH population,
      burst-never-fires-in-production)

### Reproducibility

- [ ] Exact parameter file for production checked into repo
- [ ] IC generation script or archive documented
- [ ] Enzo build configuration documented (Make.mach, configure flags,
      WINDS status)
- [ ] Random seed (if any) recorded

---

## 12. Proposed Phase Roadmap

### Phase D0b-lite

- **Goal:** Exercise untested engineering paths: burst, conservation
  through burst, restart, MPI consistency.
- **Non-goal:** No new physics ICs. No gas-limitation test. No mismatch
  quantification.
- **Files touched:** `run/BHSeed/feedback_phaseC_tests/run_feedback_phaseC_matrix.py`
  (extend with F1, F3, T7, T10). Update EXPECTED_HEAD to 685c7a9c.
- **Compute cost:** ~1-2 CPU-hours (TS3_wrap variants).
- **Acceptance criteria:** F1 gas-present burst fires using the Section 5C
  configuration and criteria. F3 conservation holds using
  `abs(conservation_residual_cgs)`. T7 restart preserves attributes per
  Section 5D. T10 diagnostics agree within Section 5E tolerance.
- **Required before production:** Yes.

### Phase D0c

- **Goal:** Controlled physics validation only, using purpose-built ICs.
  D0c includes the Bondi formula/unit check, Eddington cap activation,
  gas-limited accretion, requested-vs-realized mismatch quantification,
  seeding controlled ICs, repositioning ON/OFF, and feedback controlled
  physics tests as specified.
- **Non-goal:** No cosmological runs. No parameter calibration. No
  observational comparison pipeline implementation.
- **Files touched:** New IC generator(s); new test script
  `run/BHSeed/physics_tests/run_physics_validation.py`.
- **Compute cost:** ~10-50 CPU-hours (isolated halos, SelfGravity=1, moderate
  resolution).
- **Acceptance criteria:** A1 formula check passes. A2 cap activates. A5
  gas limitation handled. S1/S2/S3/S4/S5 seeding behavior documented as
  applicable. R4 enhancement factor quantified. F6 mismatch metric M
  measured. Any feedback controlled physics tests assigned to D0c pass their
  specified criteria.
- **Required before production:** Yes.

### Phase D2-Decision

- **Goal:** Using D0c F6 data and comparison uncertainties from R4, decide
  whether to implement D2.
- **Non-goal:** No code changes unless decision is "implement."
- **Compute cost:** 0 (analysis of existing data).
- **Required before production:** Yes (the decision, not necessarily the
  implementation).

### Phase D2-Implementation (conditional)

- **Goal:** If D2 decided, tie feedback energy to gas-realized rate.
- **Files touched:** `Grid_BHAccretionHandler.C`,
  `Grid_BHFeedbackHandler.C`, possibly `typedefs.h`.
- **Compute cost:** ~5-10 CPU-hours (rerun D0a+D0b-lite+D0c regression).
- **Required before production:** Only if D2-Decision says yes.

### Phase Pilot

- **Goal:** First cosmological run at target resolution. Measure BH
  population, calibration proximity, burst frequency, gas-limitation
  frequency, mismatch metric under realistic conditions.
- **Non-goal:** Not the production run. Not publication-quality.
- **Compute cost:** ~10,000-50,000 CPU-hours (depends on box/resolution,
  which depends on science question).
- **Required before production:** Yes.

### Phase P1 (Calibration)

- **Goal:** Tune ε_f and possibly BHSeedMass/BHFeedbackMinEnergyBurst to
  match calibration targets. 3-5 runs.
- **Compute cost:** ~50,000-150,000 CPU-hours.
- **Required before production:** Yes.

### Phase Production

- **Goal:** Publication-quality runs at calibrated parameters.
- **Compute cost:** ~100,000-500,000 CPU-hours.
- **Required before production:** N/A (this IS production).

### Phase FIG

- **Goal:** Automated analysis and observational comparison pipeline
  producing paper figures and comparison tables from simulation output. This
  includes halo finder setup, BH-halo association, galaxy property extraction,
  M_BH-M_star extraction, stellar mass/SFR catalogs, observational comparison
  tables, and figure-generation scripts.
- **Relationship to D0c:** Separate pre-production infrastructure task. D0c
  may define diagnostic requirements consumed by this pipeline, but does not
  implement it.
- **Required before production:** Observational comparison pipeline
  prototyped before production; finalized after.

### Phase PAPER

- **Goal:** Write the paper.

---

## 13. Open Questions Requiring Human Decision

These cannot be resolved by code inspection or automated testing. They
require scientific judgment.

| # | Question | Options | Impact | When to decide |
|---|----------|---------|--------|----------------|
| 1 | **Paper type** | (a) Methods paper: "we present a new BH module and demonstrate it produces reasonable galaxy populations." (b) Science paper: "we use our BH module to study [specific question]." (c) Hybrid. | Determines figure list, comparison depth, validation scope. Methods paper has lower bar for predictive accuracy but higher bar for controlled validation. | Before pilot design |
| 2 | **Primary comparison simulation** | EAGLE (closest model: stochastic thermal + AM suppression), TNG (most cited), SIMBA, or Habouzit+ 2021 cross-code | Determines which parameter choices to justify and which figures to include | Before pilot design |
| 3 | **Target science question** (if science paper) | BH seeding at high-z, AGN feedback and quenching, BH mass function, BH-galaxy co-evolution, or other | Determines box size, resolution, redshift range, analysis focus | Before pilot design |
| 4 | **Target resolution** | ~0.5-2 kpc physical at z=2 (typical for cosmological galaxy formation) | Affects all parameter calibration; boosted Bondi is resolution-dependent | Before pilot design |
| 5 | **Calibration observables** | Proposed: M_BH-M_star + SMF at z=0. Alternatives: M_BH-σ, or single-target calibration | Determines which parameters to tune and what to predict | Before calibration |
| 6 | **Simulation geometry** | Cosmological box (uniform volume), zoom (single halo at high resolution), or isolated suite (controlled but limited) | Cost/science tradeoff; zoom is cheapest for single-galaxy questions but not representative | Before pilot design |
| 7 | **Repositioning acceptability** | (a) Accept density-peak repositioning as modeling choice, quantify R4 enhancement, discuss limitations. (b) Implement dynamical friction subgrid model (major new development; out of scope for first paper). (c) Run without repositioning and document consequences. | Affects BH growth rates; referee sensitivity high | Before production, after R4 data |
| 8 | **Is D2 likely required?** | Unknown until F6 data exists. Prior: probably not for first paper (Booth & Schaye 2009 used similar approximation), but must be quantified. | Determines whether additional code development is needed before production | After D0c |

---

## Appendix: Glossary

| Term | Definition |
|------|------------|
| D0a | Phase C feedback engineering smoke test (committed: 685c7a9c) |
| D0b-lite | Engineering depth: burst, conservation, restart, MPI |
| D0c | Controlled physics tests with purpose-built ICs |
| D1 | Duty-cycle gating (deferred) |
| D2 | Tie feedback energy to gas-realized accretion rate |
| V0a | This document |
| ε_r | Radiative efficiency (BHAccretionRadiativeEfficiency) |
| ε_f | Thermal coupling efficiency (BHFeedbackThermalEfficiency) |
| f_AM | Angular momentum suppression factor |
| frac_gas | Gas-realized fraction (dm_removed / dm_requested) |
| frac_cap | Eddington cap fraction (mdot_edd_eff / mdot_raw) |
| lambda_edd | Eddington ratio (mdot_raw / mdot_edd) |
| WINDS | Compile-time macro affecting particle attribute layout |
