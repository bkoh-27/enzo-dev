# MBH Thermal Feedback Fix

**Commit:** `38140bb6` — *MBH feedback: persist lar in attr[3], enforce attr slots, clean logging*

## Problem

MBH (massive black hole) thermal feedback was silently injecting zero energy.
The accretion rate (`last_accretion_rate`) was computed correctly during the
accretion step, but by the time the feedback step read it, the value was always
zero. As a result, MBH thermal heating had no effect on the gas.

## Root Cause

Enzo rebuilds `Star` objects from raw particle data every time the AMR hierarchy
changes. The `Star` constructor reads particle attributes 0-2 (birth time,
lifetime, metallicity) but had no way to restore `last_accretion_rate`. Any
value written to that field during accretion was lost the next time grids were
reconstructed, which happens frequently in AMR simulations.

## Fix

Store `last_accretion_rate` in `ParticleAttribute[3]` so it persists across
Star object reconstruction, restarts, and load balancing.

### Changes by file

**Core persistence (read/write cycle):**

- `StarRoutines.C` — Grid-based `Star` constructor and `CopyFromParticle` now
  read `ParticleAttribute[3]` for MBH particles.
- `Star_MirrorToParticle.C` — `MirrorToParticle` now writes
  `last_accretion_rate` to `ParticleAttribute[3]` for MBH particles.
- `Star_CalculateMassAccretion.C` — Direct write of `last_accretion_rate = mdot`
  after the Bondi calculation (defense-in-depth).

**Attribute slot enforcement:**

- `InitializeNew.C` — Auto-promotes `NumberOfParticleAttributes` to 4 when
  `MBHFeedback > 0` (new simulations).
- `ReadAllData.C` — Same auto-promotion on the restart path, before grid data
  is read.

**Backward-compatible I/O:**

- `Grid_ReadGrid.C` — Both HDF4 and HDF5 read paths now handle a missing 4th
  attribute gracefully: zero-initialize the slot and print a one-time info
  message instead of crashing. This allows restarting from checkpoints written
  with only 3 attributes.

**Minor / cleanup:**

- `Star.h` — `operator=` signature changed from `Star a` (by-value) to
  `const Star &a` (by-reference) to avoid unnecessary copies.
- `StarParticleAddFeedback.C` — Added a debug-gated `[MBHFB]` log line that
  prints cycle, BH ID, mdot, radius, and energy when `debug=1`.
- `Star_CalculateFeedbackParameters.C` — Whitespace cleanup only.

### Particle attribute convention

| Index | Field | Used by |
|-------|-------|---------|
| 0 | BirthTime | All star types |
| 1 | LifeTime | All star types |
| 2 | Metallicity | All star types |
| 3 | last_accretion_rate | MBH only |

## Backward Compatibility

- Simulations without MBH feedback are unaffected.
- Restarts from 3-attribute checkpoints work: the missing slot is
  zero-initialized, so feedback begins accumulating from the next accretion
  step (one cycle of zero energy, then normal operation).
- `NumberOfParticleAttributes` is auto-promoted, so existing parameter files
  do not need to be updated unless the user wants to be explicit.

## How to Verify

Run with `debug = 1` and grep for `[MBHFB]` in stdout. A working setup shows
nonzero `mdot` and `E` values after the first accretion step:

```
[MBHFB] cycle=23 BH=100663297 mdot=6.489e-10 Msun/yr lar=6.489e-10 Radius_pc=1646.2 E=1.957e-11
```

If `mdot=0` persists after multiple cycles, check that
`NumberOfParticleAttributes >= 4` in your parameter file or restart dump.
