# Black-hole dynamical-friction diagnostics

Enzo currently includes a diagnostics-only black-hole dynamical-friction estimator.
It can compute and log collisionless dynamical-friction quantities around black
holes, but it does not apply forces, move particles, change particle masses,
update grid fields, or modify simulation state.

Active black-hole dynamical friction is not implemented in this phase. Method
values greater than 1 are invalid and remain future work.

## Status summary

This module is useful when you want to inspect the local collisionless
background around black-hole particles and see what a Tremmel-style
dynamical-friction estimate would report.

The current implementation:

- uses dark-matter particles as the collisionless background;
- excludes stars and black holes from the background estimate;
- emits structured `[BHDF]`, `[BHDF_SUMMARY]`, and `[BHDF_WARN]` log rows;
- leaves `applied=0`, `active_rows=0`, and `applied_kicks=0`;
- does not apply a velocity kick or acceleration correction;
- does not modify particle state, grid state, timestep state, or HDF5 schema.

## Quick-start parameter block

Use this block to enable detailed diagnostics:

```text
BHDynamicalFrictionMethod = 1
BHDynamicalFrictionVerbose = 2
BHDynamicalFrictionKernelRadius = 1.0
BHDynamicalFrictionMinParticles = 8
BHDynamicalFrictionMinSlowParticles = 2
```

`BHDynamicalFrictionMethod = 1` enables diagnostics only. It does not enable an
active force.

## Parameter reference

Defaults are set in `src/enzo/SetDefaultGlobalValues.C`. The parameters are
read from parameter files in `src/enzo/ReadParameterFile.C` and written to
restart/output parameter files in `src/enzo/WriteParameterFile.C`.

| Parameter | Type | Default | Allowed values | Meaning | Notes |
| --- | --- | --- | --- | --- | --- |
| `BHDynamicalFrictionMethod` | integer | `0` | `0` or `1` | Main switch. | `0` disables the handler. `1` runs diagnostics only. Values `< 0` or `> 1` are invalid. |
| `BHDynamicalFrictionVerbose` | integer | `0` | `0`, `1`, or `2` | Logging level. | `0` is quiet. `1` emits summary and warning rows when the handler reaches a BH-containing grid. `2` also emits per-BH detail rows. |
| `BHDynamicalFrictionKernelRadius` | float | `1.0` | positive when Method `> 0` | Spherical diagnostic kernel radius. | Units are physical kpc. The handler converts this to code length at runtime. |
| `BHDynamicalFrictionMinParticles` | integer | `8` | `>= 1` | Minimum dark-matter particles required in the kernel. | Rows with fewer DM particles are marked under-resolved. |
| `BHDynamicalFrictionMinSlowParticles` | integer | `2` | `>= 1` | Minimum slow dark-matter particles required. | Slow particles are selected relative to the BH/background motion. |

## Method values

- `BHDynamicalFrictionMethod = 0`: disabled. The handler is not called and no
  `[BHDF]`, `[BHDF_SUMMARY]`, or `[BHDF_WARN]` rows should be produced by this
  module.
- `BHDynamicalFrictionMethod = 1`: diagnostics-only collisionless estimator.
  The handler computes and logs diagnostic quantities but does not apply them.
- `BHDynamicalFrictionMethod >= 2`: not implemented. These values are rejected
  by parameter validation and are reserved for future work.

## Verbosity values

- `BHDynamicalFrictionVerbose = 0`: quiet. The handler can compute diagnostics
  internally, but it does not emit BHDF log rows.
- `BHDynamicalFrictionVerbose = 1`: summary and warning logging. This includes
  `[BHDF_SUMMARY]` accounting rows and `[BHDF_WARN]` condition rows when
  applicable.
- `BHDynamicalFrictionVerbose = 2`: detailed diagnostics. This includes the
  level-1 rows plus per-BH `[BHDF]` detail rows.

The validation runs used `BHDynamicalFrictionVerbose = 2`.

## What the estimator uses

For each local black-hole particle, the handler:

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

The candidate acceleration is a logged diagnostic. It is not applied to the
black hole or to any other particle.

## Log output

The module writes structured rows to the normal Enzo output stream when
verbosity is enabled.

### `[BHDF]`

Per-BH detail row emitted when `BHDynamicalFrictionVerbose >= 2`.

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

Per-grid/per-invocation accounting row emitted when
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

Warning or diagnostic-condition row emitted when
`BHDynamicalFrictionVerbose >= 1`.

Examples include:

- unsupported grid rank;
- active repositioning also enabled;
- guard skips;
- under-resolved kernel;
- no slow particles;
- kernel truncation;
- zero or large Coulomb logarithm;
- hypothetical velocity-kick cap trigger.

Warnings are diagnostic signals. They do not mean an active force was applied.

## Common checks for users

- With `BHDynamicalFrictionMethod = 0`, logs should contain no `[BHDF]`,
  `[BHDF_SUMMARY]`, or `[BHDF_WARN]` rows from this module.
- With `BHDynamicalFrictionMethod = 1` and
  `BHDynamicalFrictionVerbose = 2`, eligible black holes on local grids should
  produce `[BHDF]` detail rows and `[BHDF_SUMMARY]` accounting rows.
- In all current DF diagnostics, `[BHDF]` rows should report `applied=0`.
- `[BHDF_SUMMARY]` rows should report `active_rows=0` and `applied_kicks=0`.
- Under-resolved warnings mean the kernel did not meet the configured
  dark-matter particle-count requirement.
- `no_slow_particles` warnings mean the kernel did not meet the configured
  slow-particle requirement.
- `kernel_truncated` warnings mean the diagnostic kernel crossed the local grid
  boundary.
- `cap_triggered=1` reports that the hypothetical diagnostic velocity kick
  would exceed the internal cap monitor. It still does not apply a kick.

## Important limitations

This is a diagnostics-only module.

It does not provide:

- active dynamical-friction force;
- velocity kicks;
- particle-position, particle-velocity, particle-mass, or particle-attribute
  updates;
- grid-field updates;
- timestep coupling;
- gas dynamical friction;
- inspiral validation;
- sinking-time validation;
- dynamical-friction parameter calibration;
- a guarantee of bitwise-identical parallel cosmological evolution;
- production science certification.

Active dynamical friction remains blocked pending a new design and validation
phase. Do not interpret this documentation as authorization to enable or use
Method `>= 2`.

## Validation summary

The DF-0 diagnostics-only implementation was checked outside this repository
page with parser fixtures, build evidence, production-like diagnostic parsing,
and a deterministic one-rank collisionless invariance test.

The deterministic test showed exact physical-state equality between control
and diagnostics runs across `DD0001`, `DD0002`, and `DD0003`. That result
supports the claim that Method 1 diagnostics do not directly modify simulation
state in that controlled setup.

Production-like cosmological diagnostics parsed successfully. A separate
cosmological exact-HDF5 mismatch in the paired run found broad differences that
are retained as a numerical-divergence caveat, not as evidence of direct
dynamical-friction state writes.

## Future work

Future active dynamical friction requires a separate design authority and
validation ladder. Existing diagnostics found cap-triggered estimates in
production-like runs, so Method `>= 2` must not be enabled by this
documentation or treated as validated.

Potential future work includes:

- active-force design;
- a validation ladder for active kicks and timestep coupling;
- parameter calibration;
- gas dynamical-friction design;
- in-repository documentation or test harness updates for any future phase.
