# BH Module Collaborator Presentation Design

Date: 2026-05-12

## Goal

Create a 20-25 minute Markdown-first collaborator presentation about the Enzo black hole physics module, covering seeding, accretion, repositioning, and feedback. Each component should follow the same flow:

1. Scientific motivation or literature context.
2. What was coded and how it was coded.
3. Results, validation evidence, or explicit status caveats from the existing artifacts.

The Markdown deck will be the working artifact for collaborative editing. PowerPoint generation comes after the content is stable.

## Source Material

Primary project documentation:

- `/gpfs/bkoh/bh_dev/docs/START_HERE_BH_MODULE.md`
- `/gpfs/bkoh/bh_dev/docs/BH_MODULE_STATUS.md`
- `/gpfs/bkoh/bh_dev/docs/FEEDBACK_PHASEB_LOCKED_SPEC.md`
- `/gpfs/bkoh/bh_dev/docs/TEST_RUNBOOK_BH_MODULE.md`
- `/gpfs/bkoh/bh_dev/docs/KNOWN_ISSUES_AND_GUARDRAILS.md`

Existing deck package:

- `/gpfs/bkoh/bh_dev/presentation_deck_20260510/draft_slide_deck.md`
- `/gpfs/bkoh/bh_dev/presentation_deck_20260510/slide_outline.md`
- `/gpfs/bkoh/bh_dev/presentation_deck_20260510/speaker_notes.md`
- `/gpfs/bkoh/bh_dev/presentation_deck_20260510/figure_captions.md`
- `/gpfs/bkoh/bh_dev/presentation_deck_20260510/caveats_and_safe_language.md`
- `/gpfs/bkoh/bh_dev/presentation_deck_20260510/collaborator_QA.md`
- `/gpfs/bkoh/bh_dev/presentation_deck_20260510/rehearsal_review/`

Results and figures:

- `/gpfs/bkoh/bh_dev/presentation_plots_20260510/reports/presentation_plots_report.md`
- `/gpfs/bkoh/bh_dev/presentation_plots_20260510/figures/`
- `/gpfs/bkoh/bh_dev/presentation_plots_20260510/data/`

Implementation anchors in source:

- `src/enzo/EvolveLevel.C`: lifecycle order.
- `src/enzo/typedefs.h`: BH particle attribute slots.
- `src/enzo/Grid_BHSeedHandler.C`: seeding candidate gather, ranking, exclusion, mass removal, metadata.
- `src/enzo/Grid_BHAccretionHandler.C`: hot/cold accretion, Eddington cap, gas-limited realized bookkeeping, mass invariant checks.
- `src/enzo/Grid_BHRepositionHandler.C`: density peak diagnostics, active-zone target, rate-limited/teleport movement modes.
- `src/enzo/Grid_BHFeedbackHandler.C`: realized-basis feedback diagnostics, thermal reservoir, burst/deposition path, SF blocking.
- `src/enzo/Grid_MoveAllParticles.C`, `src/enzo/Grid_MoveSubgridParticles.C`, `src/enzo/Grid_MoveSubgridParticlesFast.C`, `src/enzo/Grid_TransferSubgridParticles.C`: absolute BH `ParticleMass` preservation across AMR transfers.

## Deck Architecture

Use an 18-slide component-block structure:

1. Title and boundary: engineering validation, not production science.
2. Scientific need for subgrid BH physics.
3. Module architecture and lifecycle order: seed -> reposition -> accrete -> feedback -> star formation.
4. Validation artifacts and figure provenance: RD0004, DD0001, plot report.
5. Seeding science motivation.
6. Seeding implementation.
7. Seeding result.
8. Accretion science motivation.
9. Accretion implementation.
10. Accretion result.
11. AMR mass safety fix and validation.
12. Repositioning science motivation.
13. Repositioning implementation.
14. Repositioning status and caveat.
15. Feedback science motivation.
16. Feedback implementation.
17. Feedback diagnostic result.
18. Evidence summary, caveats, and next validation gates.

This structure should fit a 20-25 minute walk-through if the science/context slides stay concise and detailed code/log material goes into speaker notes or backup sections.

## Component Content

### Seeding

Scientific context:

- BH seeds are subgrid stand-ins for unresolved early massive-black-hole formation channels.
- Literature anchor: Volonteri 2010 for seed formation motivation and the link between massive black holes and galaxy evolution.

Implementation content:

- Deterministic candidate collection and MPI all-gather.
- Candidate ranking with deterministic tie-breaks.
- Multi-seed acceptance walk with exclusion and de-duplication.
- Kernel-distributed active-zone gas mass removal.
- Metadata initialization for seed channel, redshift, patch properties, accept rank, formation mass, accretion fields, and feedback reservoir fields.

Result:

- Native Enzo seeding produces an increasing BH particle count in validated artifacts.
- RD0004 has 8464 BH particles at z approximately 7.999995.
- Use `bh_population_vs_redshift.png`; optionally include `RD0004_slice_z_density_bh_overlay.png` as qualitative visual context.

Required caveat:

- Not an occupation-fraction comparison.
- Not seed calibration.

### Accretion

Scientific context:

- Bondi-like accretion is the standard unresolved baseline for gas capture.
- Rotating/cold gas makes pure Bondi accretion incomplete; angular momentum can suppress inflow.
- Literature anchors: Bondi 1952; Rosas-Guevara et al. 2015/EAGLE angular momentum model.

Implementation content:

- Kernel samples gas around each BH in deterministic k-j-i order.
- Cells are classified into hot/cold channels using `t_cool / t_dyn`, with a temperature fallback.
- Hot channel uses conservative boosted Bondi-Hoyle.
- Cold channel uses Bondi with angular-momentum suppression from kernel angular momentum.
- Eddington cap is applied before active removal.
- Active-zone-only removal kernel computes gas-limited realized growth.
- Particle attributes track actual and realized rates, accreted mass, Eddington ratio, and cumulative mass invariant.
- A pre-accretion mass-sync guard rejects unsafe `ParticleMass` mismatches before gas removal.

Result:

- DD0001 has `[BHACCR] = 636457` and `[BHACCR_SRC] = 636457`.
- Positive finite actual/realized pairs align with y=x within numerical precision.
- `BHACCR_INVARIANT_FAIL = 0`.
- `BHACCR_MASS_MISMATCH_PRE = 0`.
- Use `bh_accretion_actual_vs_realized_DD0001.png`.

Required caveat:

- Code-unit bookkeeping diagnostic, not accretion-model calibration.

### AMR Mass Safety

Scientific context:

- AMR particle transfer can create bookkeeping errors if particle masses are scaled like ordinary grid-level particles.
- For BHs, absolute particle mass must remain consistent with `BHFormationMass + BHAccretedMass`.

Implementation content:

- AMR transfer paths detect `PARTICLE_TYPE_MBH` and `PARTICLE_TYPE_BLACK_HOLE`.
- Ordinary particles retain existing `MassIncrease`/`MassDecrease` behavior.
- BH particles preserve absolute `ParticleMass` across parent/child movement and transfer paths.
- Diagnostic tag: `[BH_PARTICLE_AMR_MASS_TRANSFER]`.
- Accretion precheck and postcheck enforce the mass invariant.

Result:

- DD0001 HDF5 validation job 209132 passed with failures/warnings 0/0.
- 64 CPU files, 11542 real grid groups, 16787718 particles, 8834 BH particles.
- Maximum relative residual below `1.0e-4` tolerance.
- Use `bh_mass_conservation_residuals_DD0001.png`.

Required caveat:

- Validates the requested mass-bookkeeping gate, not long-run growth calibration.

### Repositioning

Scientific context:

- Underresolved BH orbits can wander numerically away from dense gas or halo centers.
- Repositioning is a pragmatic coupling aid; it is not a full unresolved dynamical-friction model.
- Literature anchor: Tremmel et al. 2015 for why simple advection/repositioning can be artificial and why orbital treatment matters.

Implementation content:

- Deterministic density-peak search around each BH.
- Diagnostics distinguish active+ghost diagnostic peak from active-zone movement target.
- Optional potential-minimum diagnostics degrade gracefully if the potential field is unavailable.
- Production movement is rate-limited drift to the active-zone density peak.
- Teleport mode is debug-only.
- Newly seeded BHs are skipped for one pass.
- Repositioning runs before accretion in `EvolveLevel.C`.

Status:

- Implemented and instrumented.
- Stage 4 intentionally kept repositioning off to isolate accretion bookkeeping, feedback-basis diagnostics, and AMR mass safety.
- Current evidence includes `reposition_occurred=1 = 0`.

Required caveat:

- Dynamic repositioning is not validated in the current deck.
- Next step is a reposition-only controlled movement validation.

### Feedback

Scientific context:

- AGN feedback closes the self-regulation loop between BH growth and the surrounding gas.
- Common subgrid models use thermal coupling at high accretion states and kinetic/mechanical feedback at low accretion states.
- Literature anchors: Springel, Di Matteo, and Hernquist 2005; Booth and Schaye 2009; Sijacki et al. 2007; Weinberger et al. 2017.

Implementation content:

- Feedback source term is tied to realized accretion, not requested/capped accretion.
- `Grid_BHFeedbackHandler.C` reads `BHACCR_LAST_MDOT_REALIZED`.
- Missing or anomalous realized basis falls back to requested_actual with explicit warning fields.
- Phase B path is thermal-only; low-f_Edd kinetic mode is logged as `KINETIC_INACTIVE`.
- Thermal energy is accumulated in `BHFeedbackEnergyReservoir` using code/CGS conversion helpers.
- Burst threshold controls when reservoir energy is deposited.
- Deposition kernel is active-zone-only and mass-weighted.
- Cells receiving feedback are marked in the transient SF mask.
- Cumulative reservoir in/out fields track conservation in non-WINDS builds.

Result:

- DD0001 has real `[BHFDBK] = 636457`.
- All real rows use `feedback_mdot_basis=realized`.
- Requested_actual basis count is zero.
- No feedback burst or deposition occurred in the damped one-cycle Stage 4 gate.
- `E_requested = 0` and `E_deposited = 0`.
- Use `feedback_basis_summary_DD0001.png`.

Required caveat:

- Validates realized-basis feedback diagnostics only.
- Burst/deposition remains Stage 4b.
- Not feedback calibration.

## Output Files To Create

Create a new Markdown deck under:

- `/gpfs/bkoh/bh_dev/presentation_deck_20260510/bh_module_technical_walkthrough_20min.md`

Optionally create companion files after the deck draft is stable:

- `/gpfs/bkoh/bh_dev/presentation_deck_20260510/bh_module_technical_walkthrough_speaker_notes.md`
- `/gpfs/bkoh/bh_dev/presentation_deck_20260510/bh_module_technical_walkthrough_references.md`

## Required Corrections From Existing Draft

- Avoid saying the patches are on `origin/main` unless the remote branch is actually updated. Local `main` currently points at `6327b267`, while fetched `origin/main` points at `0c13b29c`.
- Replace "growing BH population" with "increasing BH particle count" where the goal is avoiding confusion with mass growth.
- Keep on-slide caveat footers for seeding, accretion, feedback, and repositioning.
- Do not imply feedback deposition, dynamic repositioning, physical accretion calibration, z=0 occupation fraction, or production science results.

## Success Criteria

- The deck follows the approved 18-slide component-block architecture.
- Every component has scientific context, code implementation detail, and result/status.
- Figures are linked to the existing `/gpfs/bkoh/bh_dev/presentation_plots_20260510/figures/` artifacts.
- Code specifics reference actual source files and stable implementation details.
- Literature references are included as short slide citations or notes.
- Caveats are visible on-slide, not only in speaker notes.
- The deck remains honest about current validation scope.

## References

- Bondi 1952, "On Spherically Symmetrical Accretion": https://academic.oup.com/mnras/article/112/2/195/2601964
- Volonteri 2010, "Formation of Supermassive Black Holes": https://arxiv.org/abs/1003.4404
- Rosas-Guevara et al. 2015, "The impact of angular momentum on black hole accretion rates in simulations of galaxy formation": https://arxiv.org/abs/1312.0598
- Tremmel et al. 2015, "Off the Beaten Path": https://arxiv.org/abs/1501.07609
- Springel, Di Matteo, and Hernquist 2005, "Modeling feedback from stars and black holes in galaxy mergers": https://arxiv.org/abs/astro-ph/0411108
- Booth and Schaye 2009, "Cosmological simulations of the growth of supermassive black holes and feedback from active galactic nuclei": https://arxiv.org/abs/0904.2572
- Sijacki et al. 2007, "A unified model for AGN feedback in cosmological simulations of structure formation": https://arxiv.org/abs/0705.2238
- Weinberger et al. 2017, "Supermassive black holes and their feedback effects in the IllustrisTNG simulation": https://arxiv.org/abs/1710.04659
