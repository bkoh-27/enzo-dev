# Black-hole dynamical-friction diagnostics

Enzo currently includes black-hole dynamical-friction diagnostics and an active
logging scaffold. The implemented method states are:

- Method `0`: disabled/no-output.
- Method `1`: diagnostics-only Rev1 collisionless estimator.
- Method `2`: active-schema dry-run logging scaffold only.

Active black-hole dynamical-friction force, active kicks, `ActiveApply=1`, and
active DF validation remain blocked. The Method 2 path emits active-format log
rows for parser and runtime plumbing evidence, but it does not apply forces,
move particles, change particle masses, update grid fields, modify timesteps,
or change HDF5 outputs.

`BLOCKING_FOR_FUTURE_ACTIVE_DF=true` remains.

## Status Summary

This module is useful when you want to inspect the local collisionless
background around black-hole particles and see either the Rev1 diagnostic
estimate or the newer active-schema dry-run rows.

The current implementation:

- uses dark-matter particles as the collisionless background for Method 1;
- excludes stars and black holes from the Method 1 background estimate;
- emits structured Rev1 `[BHDF]`, `[BHDF_SUMMARY]`, and `[BHDF_WARN]` rows for
  Method 1 when `BHDynamicalFrictionVerbose` permits;
- emits structured `[BHDF_ACTIVE]` and `[BHDF_ACTIVE_SUMMARY]` dry-run rows for
  Method 2 when `BHDynamicalFrictionActiveVerbose` permits;
- emits no `[BHDF_ACTIVE_WARN]` rows by default;
- keeps Method 2 separated from Rev1 `[BHDF]` tags;
- leaves `active_applied=0`, `kick_applied=0`, and
  `applied_delta_v_code=0` in Method 2 detail rows;
- uses `suppression_reason_name=diagnostics_only` and
  `active_validation_flags=dry_run_no_physics` for Method 2 detail rows;
- does not apply a velocity kick or acceleration correction;
- does not modify particle state, grid state, timestep state, or HDF5 schema.

## Quick-Start Parameter Blocks

Use this block to enable detailed Rev1 diagnostics:

```text
BHDynamicalFrictionMethod = 1
BHDynamicalFrictionVerbose = 2
BHDynamicalFrictionKernelRadius = 1.0
BHDynamicalFrictionMinParticles = 8
BHDynamicalFrictionMinSlowParticles = 2
```

Use this block only for Method 2 dry-run active logging:

```text
BHDynamicalFrictionMethod = 2
BHDynamicalFrictionActiveApply = 0
BHDynamicalFrictionActiveSchemaVersion = 2
BHDynamicalFrictionActiveVerbose = 2
BHDynamicalFrictionKernelRadius = 1.0
```

`BHDynamicalFrictionMethod = 1` enables diagnostics only.
`BHDynamicalFrictionMethod = 2` enables dry-run active-schema logging only.
Neither mode enables an active force.

## Method Values

| Method | State | Log families | State mutation | Notes |
| --- | --- | --- | --- | --- |
| `0` | off | none from this module | none | Disabled/no-output control. |
| `1` | diagnostics-only Rev1 estimator | `[BHDF]`, `[BHDF_SUMMARY]`, `[BHDF_WARN]` when verbose | none | Computes and logs diagnostic quantities only. |
| `2` | active dry-run logging scaffold | `[BHDF_ACTIVE]`, `[BHDF_ACTIVE_SUMMARY]` when active verbose | none | Parser-valid active-schema rows, no force, no kick. |
| `>= 3` | invalid/reserved | none | none | Rejected by parameter validation unless future source changes define it. |

Validation currently allows `BHDynamicalFrictionMethod = 0`, `1`, or `2`.
Values below `0` or above `2` are invalid in the pushed SKELETON-0B source.

## Parameter Reference

Defaults are set in `src/enzo/SetDefaultGlobalValues.C`. The parameters are
read from parameter files in `src/enzo/ReadParameterFile.C` and written to
restart/output parameter files in `src/enzo/WriteParameterFile.C`.

### Base Parameters

| Parameter | Type | Default | Allowed values | Meaning | Notes |
| --- | --- | --- | --- | --- | --- |
| `BHDynamicalFrictionMethod` | integer | `0` | `0`, `1`, or `2` | Main switch. | `0` disables output. `1` runs Rev1 diagnostics. `2` runs active dry-run logging only. |
| `BHDynamicalFrictionVerbose` | integer | `0` | `0`, `1`, or `2` | Method 1 logging level. | `0` is quiet. `1` emits Rev1 summary/warning rows. `2` also emits Rev1 per-BH detail rows. |
| `BHDynamicalFrictionKernelRadius` | float | `1.0` | positive when Method `> 0` | Spherical diagnostic/logging kernel radius. | Units are physical kpc. The handler converts this to code length at runtime. |
| `BHDynamicalFrictionMinParticles` | integer | `8` | `>= 1` | Minimum dark-matter particles required in the Method 1 kernel. | Rows with fewer DM particles are marked under-resolved. |
| `BHDynamicalFrictionMinSlowParticles` | integer | `2` | `>= 1` | Minimum slow dark-matter particles required by Method 1. | Slow particles are selected relative to the BH/background motion. |

### Active Dry-Run Parameters

The following parameters are plumbed, defaulted, read, validated, and written
for SKELETON-0A/0B. Most are future active-force controls; in the current
pushed source they do not activate force calculation, subcycling, local-gravity
guards, momentum exchange, or kicks.

| Parameter | Default | Current status |
| --- | ---: | --- |
| `BHDynamicalFrictionActiveApply` | `0` | Must remain `0`; any nonzero value is fatal/blocked. |
| `BHDynamicalFrictionActiveSchemaVersion` | `2` | Must be `2`; used in active dry-run rows. |
| `BHDynamicalFrictionActiveVerbose` | `2` | Method 2 active logging level: `0` quiet, `1` summary, `2` detail plus summary. |
| `BHDynamicalFrictionActiveMaxKickFraction` | `0.01` | Future kick cap parameter; currently logged as dry-run metadata only. |
| `BHDynamicalFrictionActiveMaxLevelKickFraction` | `0.05` | Future per-level cap parameter; currently logged as dry-run metadata only. |
| `BHDynamicalFrictionActiveUseSubcycling` | `1` | Future subcycling control; no subcycling is performed. |
| `BHDynamicalFrictionActiveMaxSubcycles` | `64` | Future subcycle limit; currently logged as dry-run metadata only. |
| `BHDynamicalFrictionActiveMaxDtOverDtDF` | `10.0` | Future timestep/DF guard; no timestep coupling is performed. |
| `BHDynamicalFrictionActiveRequireKernelComplete` | `1` | Future guard; Method 2 dry-run rows use neutral parser-valid placeholders. |
| `BHDynamicalFrictionActiveMinParticles` | `64` | Future active minimum; no active candidate acceleration is computed. |
| `BHDynamicalFrictionActiveMinSlowParticles` | `16` | Future active slow-particle minimum; no active candidate acceleration is computed. |
| `BHDynamicalFrictionActiveVrelFloorKmS` | `0.0` | Future active guard; no active force is computed. |
| `BHDynamicalFrictionActiveSigmaFloorKmS` | `0.0` | Future active guard; no active force is computed. |
| `BHDynamicalFrictionActiveLnLambdaMax` | `-1.0` | Future active limiter; no active force is computed. |
| `BHDynamicalFrictionActiveRequireLocalGravity` | `1` | Future local-gravity guard; dry-run rows use neutral local-gravity placeholders. |
| `BHDynamicalFrictionActiveMaxADFOverAGrav` | `0.1` | Future local-gravity ratio guard; no active acceleration is computed. |
| `BHDynamicalFrictionActiveMomentumPolicy` | `0` | Must be `0` (`bh_only`) in the current skeleton. |
| `BHDynamicalFrictionActiveSuppressNewSeeds` | `1` | Future active suppression guard; no active kick path exists. |
| `BHDynamicalFrictionActiveSuppressRepositionConflict` | `1` | Future active suppression guard; no active kick path exists. |

`BHDynamicalFrictionActiveKernelMode` is deferred and rejected if present in a
parameter file.

## Verbosity Values

`BHDynamicalFrictionVerbose` controls Method 1 Rev1 diagnostics:

- `0`: quiet; no Rev1 BHDF log rows.
- `1`: Rev1 summary and warning rows.
- `2`: Rev1 summary, warning, and per-BH detail rows.

`BHDynamicalFrictionActiveVerbose` controls Method 2 active dry-run logging:

- `0`: quiet; no active tags.
- `1`: `[BHDF_ACTIVE_SUMMARY]` rows only.
- `2`: `[BHDF_ACTIVE]` detail rows plus matching
  `[BHDF_ACTIVE_SUMMARY]` rows.

## Method 1 Rev1 Diagnostics

For each local black-hole particle, Method 1:

1. converts `BHDynamicalFrictionKernelRadius` from physical kpc to code units;
2. sorts black holes deterministically by particle ID;
3. skips non-authoritative or newly seeded particles when appropriate;
4. counts particles inside the spherical kernel;
5. excludes star particles and other black holes from the background;
6. uses dark-matter particles to compute a mass-weighted center-of-mass
   velocity;
7. computes a mass-weighted one-dimensional velocity dispersion, `sigma_1d`;
8. identifies slow dark-matter particles whose speed relative to the local
   background is below the BH/background relative speed;
9. estimates `b_min`, `b_max`, the Coulomb logarithm, candidate acceleration,
   diagnostic timestep, and velocity-kick cap fraction;
10. records whether the kernel is complete, under-resolved, slow-particle
    deficient, or otherwise rejected.

The candidate acceleration is a logged Method 1 diagnostic. It is not applied
to the black hole or to any other particle.

### `[BHDF]`

Per-BH Method 1 detail row emitted when `BHDynamicalFrictionVerbose >= 2`.

Useful fields include:

- `bh_id`, `particle_index`, `bh_mass_msun`;
- `kernel_radius_code`, `kernel_radius_phys_kpc`, `kernel_complete`;
- `neighbor_particles_total`, `n_dm_in_kernel`, `n_star_excluded`,
  `n_bh_excluded`, `n_slow`;
- `v_CoM_*_kms`, `sigma_1d_kms`, `ln_lambda`;
- `a_DF_*_cgs`, `a_DF_mag_cgs`, `dt_DF_Myr`;
- `velocity_kick_cap_fraction`, `cap_triggered`, `under_resolved`,
  `rejection_reason`;
- `mode=diagnostic`, `applied=0`.

### `[BHDF_SUMMARY]`

Per-grid/per-invocation Method 1 accounting row emitted when
`BHDynamicalFrictionVerbose >= 1`.

Important summary fields include:

- `rows` and `eligible_bh`;
- `valid_or_conditioned_rows`;
- `diagnostics_only_rows`;
- `active_rows=0`;
- `applied_kicks=0`;
- skip and warning counts;
- `detail_rows_emitted` and `detail_rows_complete`.

### `[BHDF_WARN]`

Method 1 warning or diagnostic-condition row emitted when
`BHDynamicalFrictionVerbose >= 1`.

Examples include unsupported grid rank, active repositioning also enabled,
guard skips, under-resolved kernels, no slow particles, kernel truncation, zero
or large Coulomb logarithm, and hypothetical velocity-kick cap triggers.

Warnings are diagnostic signals. They do not mean an active force was applied.

## Method 2 Active Dry-Run Logging

Method 2 is an active-schema logging scaffold. It does not run the Method 1
Rev1 estimator and does not emit Rev1 `[BHDF]`, `[BHDF_SUMMARY]`, or
`[BHDF_WARN]` rows.

With `BHDynamicalFrictionActiveVerbose >= 2`, Method 2 emits one
`[BHDF_ACTIVE]` detail row per eligible local BH-type particle. These rows use
dry-run identity fields:

- `schema_version=2`;
- `mode=active_dry_run`;
- `active_method=2`;
- `active_policy=dry_run`;
- `suppression_reason_code=1`;
- `suppression_reason_name=diagnostics_only`;
- `active_applied=0`;
- `kick_applied=0`;
- `applied_delta_v_code=0`;
- `active_validation_flags=dry_run_no_physics`.

All candidate acceleration, limited acceleration, kick, application,
bookkeeping, and momentum-transfer fields are zero or neutral parser-valid
placeholders. Method 2 does not compute an active force formula, candidate
acceleration, cap-limited acceleration, subcycle application, or local-gravity
force ratio.

With `BHDynamicalFrictionActiveVerbose >= 1`, Method 2 emits
`[BHDF_ACTIVE_SUMMARY]` accounting rows. Required dry-run accounting is:

- `active_rows = detail_rows_emitted`;
- `applied_kicks=0`;
- `suppressed_rows = active_rows`;
- `diagnostics_only_rows = suppressed_rows`;
- `summary_accounting_passed=1`.

No `[BHDF_ACTIVE_WARN]` rows are emitted by default.

## Common Checks For Users

- Method `0`: logs should contain no `[BHDF]`, `[BHDF_SUMMARY]`,
  `[BHDF_WARN]`, `[BHDF_ACTIVE]`, `[BHDF_ACTIVE_SUMMARY]`, or
  `[BHDF_ACTIVE_WARN]` rows from this module.
- Method `1`: Rev1 `[BHDF]`, `[BHDF_SUMMARY]`, and `[BHDF_WARN]` rows may
  appear when `BHDynamicalFrictionVerbose` permits, but no active tags should
  appear.
- Method `2`: `[BHDF_ACTIVE]` and `[BHDF_ACTIVE_SUMMARY]` rows may appear when
  `BHDynamicalFrictionActiveVerbose` permits, but Rev1 `[BHDF]`,
  `[BHDF_SUMMARY]`, and `[BHDF_WARN]` rows should not appear.
- Method `2` detail rows must retain `active_applied=0`, `kick_applied=0`,
  `applied_delta_v_code=0`, `suppression_reason_name=diagnostics_only`, and
  `active_validation_flags=dry_run_no_physics`.
- Under-resolved Method 1 warnings mean the kernel did not meet the configured
  dark-matter particle-count requirement.
- `no_slow_particles` Method 1 warnings mean the kernel did not meet the
  configured slow-particle requirement.
- `kernel_truncated` Method 1 warnings mean the diagnostic kernel crossed the
  local grid boundary.
- Method 1 `cap_triggered=1` reports that the hypothetical diagnostic velocity
  kick would exceed the internal cap monitor. It still does not apply a kick.

## Validation Summary

DF-0 Method 1 diagnostics were previously checked with parser fixtures, build
evidence, production-like diagnostic parsing, and deterministic one-rank
collisionless invariance evidence.

SKELETON-0B Method 2 dry-run active logging has the following accepted
evidence:

- compile/link verification passed;
- accepted parser package checksum passed;
- synthetic source-format parser evidence passed;
- runtime Method 2 dry-run parser evidence passed;
- runtime job `322573` completed with Slurm state `COMPLETED` and exit `0:0`;
- Method 0 control emitted no DF tags;
- Method 2 dry-run emitted parser-valid active rows;
- Method 2 dry-run emitted `detail_rows=3`;
- Method 2 dry-run emitted `summary_rows=3`;
- `active_applied_count=0`;
- `kick_applied_count=0`;
- `applied_delta_v_code_all_zero=true`;
- `suppression_reason_only_diagnostics_only=true`;
- `summary_accounting_PASS=true`;
- exact Method 0 versus Method 2 physical-state equality passed across
  `DD0001`, `DD0002`, and `DD0003`;
- `exact_physical_state_equal=true`;
- `mismatch_count=0`;
- `particle_state_mismatch=false`;
- protected input integrity passed.

This evidence validates dry-run log emission, accepted parser compatibility on
runtime logs, zero active application, and exact state equality in one
deterministic single-rank collisionless test. It does not validate active force
physics.

## Important Limitations

This module does not currently provide:

- active dynamical-friction force;
- active velocity kicks;
- `ActiveApply=1` authorization;
- particle-position, particle-velocity, particle-mass, or particle-attribute
  updates;
- grid-field updates;
- timestep coupling;
- HDF5 output/schema changes;
- gas dynamical friction;
- inspiral validation;
- sinking-time validation;
- active dynamical-friction parameter calibration;
- active force validation;
- active DF physics validation;
- physical correctness validation for active DF;
- active kick validation;
- a guarantee beyond the deterministic single-rank collisionless test;
- a bitwise guarantee for parallel or cosmological runs;
- production science readiness.

Active dynamical friction remains blocked pending a new design and validation
phase. Do not interpret Method 2 dry-run logging as authorization to enable
`BHDynamicalFrictionActiveApply=1` or to use active kicks.

`BLOCKING_FOR_FUTURE_ACTIVE_DF=true` remains.

## Future Work

Future active dynamical friction requires a separate design authority and
validation ladder. Any future active kick requires new design review,
parser-valid active evidence, isolated validation, runtime validation, and
separate user authorization before `ActiveApply=1`, active force, or active
kick work can proceed.

Potential future work includes:

- active-force design;
- a validation ladder for active kicks and timestep coupling;
- active-force parser and runtime evidence;
- parameter calibration;
- gas dynamical-friction design;
- in-repository documentation or test harness updates for any future phase.
