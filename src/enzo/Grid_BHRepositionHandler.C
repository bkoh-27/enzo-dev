/***********************************************************************
/
/  GRID CLASS (BH REPOSITION DIAGNOSTICS, PHASE A)
/
/  PURPOSE:
/    Measure BH offsets to local gas density peaks (and optional potential
/    minimum) with deterministic traversal. Phase A is diagnostics-only:
/    no BH movement, no velocity changes, and no grid-field writes.
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <algorithm>
#include <vector>

#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
double ReturnWallTime();

static int BHRepositionIsBHType(int type)
{
  return (type == PARTICLE_TYPE_MBH || type == PARTICLE_TYPE_BLACK_HOLE);
}

static int BHRepositionIsNewlySeededThisPass(PINT particle_id)
{
  return (particle_id == INT_UNDEFINED);
}

struct BHRepositionParticleOrder {
  PINT *ParticleIDs;
  BHRepositionParticleOrder(PINT *ids) : ParticleIDs(ids) { }
  bool operator()(int a, int b) const {
    if (ParticleIDs[a] < ParticleIDs[b])
      return true;
    if (ParticleIDs[a] > ParticleIDs[b])
      return false;
    return (a < b);
  }
};

static int BHRepositionPositionCoveredByChild(HierarchyEntry *SubgridPointer,
                                              FLOAT *pos)
{
  for (HierarchyEntry *sub = SubgridPointer; sub != NULL;
       sub = sub->NextGridThisLevel)
    if (sub->GridData != NULL && sub->GridData->PointInGrid(pos))
      return TRUE;

  return FALSE;
}

static int BHRepositionNonAuthoritativeMassRatio(
  int num_particle_attributes,
  float *particle_attribute[],
  int p,
  double particle_mass,
  double *expected_full_mass,
  double *ratio)
{
  if (expected_full_mass != NULL)
    *expected_full_mass = -1.0;
  if (ratio != NULL)
    *ratio = -1.0;

  if (num_particle_attributes <= PARTICLE_ATTRIBUTE_BH_FORMATION_MASS ||
      num_particle_attributes <= PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS ||
      particle_attribute[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS] == NULL ||
      particle_attribute[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS] == NULL)
    return FALSE;

  if (!isfinite(particle_mass) || particle_mass <= 0.0)
    return FALSE;

  const double formation_mass =
    particle_attribute[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS][p];
  const double accreted_mass =
    particle_attribute[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS][p];

  if (!isfinite(formation_mass) || formation_mass <= 0.0 ||
      !isfinite(accreted_mass) || accreted_mass < 0.0)
    return FALSE;

  const double full_mass = formation_mass + accreted_mass;
  if (!isfinite(full_mass) || full_mass <= 0.0)
    return FALSE;

  const double mass_ratio = full_mass / particle_mass;
  if (expected_full_mass != NULL)
    *expected_full_mass = full_mass;
  if (ratio != NULL)
    *ratio = mass_ratio;

  return (isfinite(mass_ratio) && mass_ratio > 1.5) ? TRUE : FALSE;
}

static int BHRepositionLexicographicLess(const FLOAT *lhs, const FLOAT *rhs,
                                         int rank)
{
  for (int dim = 0; dim < rank; dim++) {
    if (lhs[dim] < rhs[dim])
      return TRUE;
    if (lhs[dim] > rhs[dim])
      return FALSE;
  }
  return FALSE;
}

static int BHRepositionPositionEqual(const FLOAT *lhs, const FLOAT *rhs, int rank)
{
  for (int dim = 0; dim < rank; dim++)
    if (lhs[dim] != rhs[dim])
      return FALSE;
  return TRUE;
}

static double BHRepositionSearchRadiusCode(float SearchRadiusPhysKpc,
                                           FLOAT time,
                                           float LengthUnits)
{
  float h_param = 1.0f;
  float a_phys = 1.0f;
  float box_kpch = 0.0f;
  double search_radius_code = 0.0;

  if (ComovingCoordinates) {
    FLOAT a = 1.0, dadt = 0.0;
    CosmologyComputeExpansionFactor(time, &a, &dadt);
    a_phys = float(a / (1.0 + InitialRedshift));

    h_param = HubbleConstantNow;
    if (h_param > 10.0f)
      h_param *= 0.01f;
    if (h_param <= 0.0f)
      h_param = 1.0f;
    if (a_phys <= 0.0f)
      a_phys = 1.0f;

    box_kpch = ComovingBoxSize * 1000.0f;
    const float search_com_kpch = SearchRadiusPhysKpc * h_param / a_phys;
    if (box_kpch > 0.0f)
      search_radius_code = search_com_kpch / box_kpch;
  } else {
    box_kpch = LengthUnits / kpc_cm;
    if (box_kpch > 0.0f)
      search_radius_code = SearchRadiusPhysKpc / box_kpch;
  }

  return search_radius_code;
}

int grid::BHRepositionDiagnosticHandler(HierarchyEntry* SubgridPointer,
                                        int level, int cycle_number,
                                        float dtLevelAbove)
{
  (void) dtLevelAbove;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0 || NumberOfParticles <= 0)
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL)
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");

  float DensityUnits = 1.0f, LengthUnits = 1.0f, TemperatureUnits = 1.0f;
  float TimeUnits = 1.0f, VelocityUnits = 1.0f;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL)
    ENZO_FAIL("Error in GetUnits.");

  const float cell_width = float(CellWidth[0][0]);
  if (cell_width <= 0.0f)
    return SUCCESS;

  const double search_radius_code =
    BHRepositionSearchRadiusCode(BHRepositionSearchRadius, Time, LengthUnits);
  if (search_radius_code <= 0.0)
    return SUCCESS;

  const double search_radius_over_dx = search_radius_code / cell_width;
  const double kpc_per_code = BHRepositionSearchRadius / search_radius_code;
  const double effective_active_radius_code =
    min(search_radius_code, double(NumberOfGhostZones) * double(cell_width));
  const double effective_active_radius_cells =
    effective_active_radius_code / double(cell_width);
  const int active_radius_was_capped =
    (effective_active_radius_code < search_radius_code) ? 1 : 0;

  FLOAT a = 1.0, dadt = 0.0;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  const double zred = ComovingCoordinates ?
    double((1.0 + InitialRedshift) / a - 1.0) : 0.0;

  const int nx = GridDimension[0];
  const int ny = GridDimension[1];
  const int nz = GridDimension[2];

  const int rcell = max(0, int(ceil(search_radius_code / cell_width)));
  const int rcell_effective = min(rcell, NumberOfGhostZones);
  const double r2 = search_radius_code * search_radius_code;
  const double r2_effective =
    effective_active_radius_code * effective_active_radius_code;
  const int search_exceeds_ghost =
    (search_radius_over_dx > double(NumberOfGhostZones)) ? 1 : 0;

  const int isx = this->GetGridStartIndex(0);
  const int isy = this->GetGridStartIndex(1);
  const int isz = this->GetGridStartIndex(2);
  const int iex = this->GetGridEndIndex(0);
  const int iey = this->GetGridEndIndex(1);
  const int iez = this->GetGridEndIndex(2);

  FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;

  std::vector<int> bh_particles;
  bh_particles.reserve(NumberOfParticles);
  for (int p = 0; p < NumberOfParticles; p++)
    if (BHRepositionIsBHType(ParticleType[p]))
      bh_particles.push_back(p);

  if (bh_particles.empty())
    return SUCCESS;

  static int bhrepos_invocation_seq_counter = 0;
  const int bhrepos_invocation_seq = ++bhrepos_invocation_seq_counter;

  const int bhrepos_processor = MyProcessorNumber;
  const int bhrepos_grid_id = this->GetGridID();

  int bhrepos_summary_rows = 0;
  int bhrepos_summary_moved = 0;
  int bhrepos_summary_method0_rows = 0;
  int bhrepos_summary_diagnostic_only_rows = 0;
  int bhrepos_summary_active_radius_was_capped = 0;
  int bhrepos_summary_active_target_exists = 0;
  int bhrepos_summary_effective_active_target_exists = 0;
  int bhrepos_summary_effective_reposition_rejected = 0;
  int bhrepos_summary_legacy_active_reposition_rejected = 0;
  int bhrepos_summary_search_exceeds_ghost = 0;
  int bhrepos_summary_search_truncated_by_grid = 0;
  int bhrepos_summary_effective_search_truncated_by_grid = 0;
  int bhrepos_summary_under_resolved = 0;
  int bhrepos_summary_no_effective_target = 0;
  int bhrepos_summary_newly_seeded = 0;
  int bhrepos_detail_rows_emitted = 0;

  double bhrepos_summary_max_displacement_cells = 0.0;
  double bhrepos_summary_max_effective_active_peak_offset_kpc = 0.0;
  double bhrepos_summary_max_search_radius_over_dx = search_radius_over_dx;

  const int bhrepos_detail_rows_complete =
    (BHRepositionVerbose >= 2) ? 1 : 0;

  int bhrepos_legacy_reject_warn_count = 0;
  const int bhrepos_legacy_reject_warn_first_n = 10;

  const int potential_requested = (BHRepositionDiagnosePotential == 1);
  int potential_available = (PotentialField != NULL);
  int grav_dim[MAX_DIMENSION] = {1, 1, 1};
  int pot_off[MAX_DIMENSION] = {0, 0, 0};
  long long potential_size = 1;
  if (potential_available) {
    for (int dim = 0; dim < GridRank; dim++) {
      grav_dim[dim] = GravitatingMassFieldDimension[dim];
      if (grav_dim[dim] <= 0)
        potential_available = FALSE;
    }
    if (potential_available) {
      for (int dim = 0; dim < GridRank; dim++) {
        pot_off[dim] = (grav_dim[dim] - GridDimension[dim]) / 2;
        potential_size *= (long long) grav_dim[dim];
      }
    }
  }

  if (BHRepositionVerbose >= 1) {
    if (search_radius_over_dx < 1.5) {
      fprintf(logptr,
              "[BHREPOS_WARN] step=%d level=%d search_radius_over_dx=%.6g "
              "kernel under-resolved (<1.5 cells).\n",
              cycle_number, level, search_radius_over_dx);
    }
	    if (search_radius_over_dx > 3.0) {
	      fprintf(logptr,
	              "[BHREPOS_WARN] step=%d level=%d search_radius_over_dx=%.6g "
	              "kernel may exceed ghost-zone support (>3 cells).\n",
	              cycle_number, level, search_radius_over_dx);
	    }
	    if (BHRepositionMethod == 0) {
	      fprintf(logptr,
	              "[BHREPOS_WARN] step=%d level=%d BHRepositionMethod=0 "
	              "runs diagnostics only because BHRepositionVerbose=%d; "
	              "no particle movement will be applied.\n",
	              cycle_number, level, BHRepositionVerbose);
	    }
    if (BHRepositionMaxDisplacement > search_radius_over_dx) {
      fprintf(logptr,
              "[BHREPOS_WARN] step=%d level=%d max_displacement_cells=%.6g "
              "exceeds search radius in cell widths=%.6g.\n",
              cycle_number, level, double(BHRepositionMaxDisplacement),
              search_radius_over_dx);
    }
  }

  int warned_potential_unavailable = FALSE;
  if (potential_requested && !potential_available) {
    fprintf(logptr,
            "[BHREPOS_WARN] step=%d level=%d "
            "Potential field not available for repositioning diagnostics.\n",
            cycle_number, level);
    warned_potential_unavailable = TRUE;
  }

  std::sort(bh_particles.begin(), bh_particles.end(),
            BHRepositionParticleOrder(ParticleNumber));

  for (int bp = 0; bp < int(bh_particles.size()); bp++) {
    const int p = bh_particles[bp];
    const double t0 = ReturnWallTime();

    FLOAT bh_pos[MAX_DIMENSION] = {0.0, 0.0, 0.0};
    bh_pos[0] = ParticlePosition[0][p];
    bh_pos[1] = (GridRank > 1) ? ParticlePosition[1][p] : 0.0;
    bh_pos[2] = (GridRank > 2) ? ParticlePosition[2][p] : 0.0;

    const double bh_mass_code = ParticleMass[p];

    if (SubgridPointer != NULL &&
        BHRepositionPositionCoveredByChild(SubgridPointer, bh_pos)) {
      if (BHRepositionVerbose >= 1) {
        double expected_full_mass = -1.0;
        double ratio = -1.0;
        BHRepositionNonAuthoritativeMassRatio(
          NumberOfParticleAttributes, ParticleAttribute, p, bh_mass_code,
          &expected_full_mass, &ratio);
        fprintf(logptr,
                "[BHREPOS_SKIP_NON_AUTHORITATIVE] step=%d level=%d grid_id=%d "
                "bh_id=%lld particle_index=%d ParticleMass=%.15e "
                "expected_full_mass=%.15e ratio=%.15e reason=child_covered\n",
                cycle_number, level, this->GetGridID(),
                (long long) ParticleNumber[p], p, bh_mass_code,
                expected_full_mass, ratio);
      }
      continue;
    }

    if (!isfinite(bh_mass_code) || bh_mass_code <= 0.0)
      continue;

    if (!this->PointInGrid(bh_pos))
      continue;

    double expected_full_mass = -1.0;
    double mass_attribute_ratio = -1.0;
    if (BHRepositionNonAuthoritativeMassRatio(
          NumberOfParticleAttributes, ParticleAttribute, p, bh_mass_code,
          &expected_full_mass, &mass_attribute_ratio)) {
      if (BHRepositionVerbose >= 1) {
        fprintf(logptr,
                "[BHREPOS_SKIP_NON_AUTHORITATIVE] step=%d level=%d grid_id=%d "
                "bh_id=%lld particle_index=%d ParticleMass=%.15e "
                "expected_full_mass=%.15e ratio=%.15e reason=mass_attribute_ratio\n",
                cycle_number, level, this->GetGridID(),
                (long long) ParticleNumber[p], p, bh_mass_code,
                expected_full_mass, mass_attribute_ratio);
      }
      continue;
    }

    if (BHRepositionIsNewlySeededThisPass(ParticleNumber[p])) {
      const double reposition_wall_ms = 1000.0 * (ReturnWallTime() - t0);
      bhrepos_summary_rows++;
      bhrepos_summary_method0_rows +=
        (BHRepositionMethod == 0) ? 1 : 0;
      bhrepos_summary_diagnostic_only_rows +=
        (BHRepositionMethod == 0) ? 1 : 0;
      bhrepos_summary_active_radius_was_capped +=
        active_radius_was_capped ? 1 : 0;
      bhrepos_summary_newly_seeded++;
      bhrepos_summary_under_resolved +=
        (search_radius_over_dx < 1.5) ? 1 : 0;
      bhrepos_summary_search_exceeds_ghost +=
        search_exceeds_ghost ? 1 : 0;
      bhrepos_summary_no_effective_target++;

      if (BHRepositionVerbose >= 2) {
        fprintf(logptr,
                "[BHREPOS] processor=%d invocation_seq=%d grid_id=%d "
                "step=%d level=%d z=%.8g bh_id=%lld "
                "bh_pos_x=%.15g bh_pos_y=%.15g bh_pos_z=%.15g "
                "bh_density=%.8e "
                "diag_peak_pos_x=%.15g diag_peak_pos_y=%.15g diag_peak_pos_z=%.15g "
                "diag_peak_density=%.8e diag_peak_in_ghost=%d diag_peak_offset_kpc=%.8g "
                "active_peak_pos_x=%.15g active_peak_pos_y=%.15g active_peak_pos_z=%.15g "
                "active_peak_density=%.8e active_peak_offset_kpc=%.8g active_target_exists=%d "
                "offset_from_potential_kpc=%.8g "
                "displacement_kpc=%.8g displacement_cells=%.8g "
                "reposition_occurred=%d reposition_clamped=%d newly_seeded_skip=%d "
                "method=%d diagnostic_only=%d search_kernel_truncated=%d "
                "active_reposition_rejected=%d "
                "search_exceeds_ghost=%d search_truncated_by_grid=%d "
                "search_radius_over_dx=%.8g under_resolved=%d "
                "search_cells=%d search_active_cells=%d reposition_wall_ms=%.4f "
                "effective_active_radius_cells=%.8g active_radius_was_capped=%d "
                "effective_active_target_exists=%d "
                "effective_active_peak_density=%.8e "
                "effective_active_peak_offset_kpc=%.8g "
                "effective_search_truncated_by_grid=%d "
                "effective_reposition_rejected=%d rejection_reason=%d\n",
                bhrepos_processor, bhrepos_invocation_seq, bhrepos_grid_id,
                cycle_number, level, zred, (long long) ParticleNumber[p],
                bh_pos[0], bh_pos[1], bh_pos[2],
                0.0,
                bh_pos[0], bh_pos[1], bh_pos[2],
                0.0, 0, -1.0,
                bh_pos[0], bh_pos[1], bh_pos[2],
                0.0, -1.0, 0,
                -1.0,
                0.0, 0.0,
                0, 0, 1,
                BHRepositionMethod, (BHRepositionMethod == 0) ? 1 : 0,
                0, 0,
                search_exceeds_ghost, 0,
                search_radius_over_dx,
                (search_radius_over_dx < 1.5) ? 1 : 0,
                0, 0, reposition_wall_ms,
                effective_active_radius_cells, active_radius_was_capped,
                0, 0.0, -1.0, 0, 0, 6);
        bhrepos_detail_rows_emitted++;
      }
      continue;
    }

    int i0 = int(floor((ParticlePosition[0][p] - CellLeftEdge[0][0]) / cell_width));
    int j0 = int(floor((ParticlePosition[1][p] - CellLeftEdge[1][0]) / cell_width));
    int k0 = int(floor((ParticlePosition[2][p] - CellLeftEdge[2][0]) / cell_width));
    i0 = max(0, min(nx - 1, i0));
    j0 = max(0, min(ny - 1, j0));
    k0 = max(0, min(nz - 1, k0));

    const int host_index = (k0*ny + j0)*nx + i0;
    double bh_density = BaryonField[DensNum][host_index];
    if (!isfinite(bh_density) || bh_density < 0.0)
      bh_density = 0.0;

    double diag_peak_density = bh_density;
    double active_peak_density = bh_density;
    double effective_active_peak_density = bh_density;

    FLOAT diag_peak_pos[MAX_DIMENSION] = {bh_pos[0], bh_pos[1], bh_pos[2]};
    FLOAT active_peak_pos[MAX_DIMENSION] = {bh_pos[0], bh_pos[1], bh_pos[2]};
    FLOAT effective_active_peak_pos[MAX_DIMENSION] =
      {bh_pos[0], bh_pos[1], bh_pos[2]};

    int diag_peak_in_ghost = 0;
    int diag_peak_from_cell = FALSE;
    int active_peak_from_cell = FALSE;
    int effective_active_peak_from_cell = FALSE;

    int search_cells = 0;
    int search_active_cells = 0;

    int potential_min_found = FALSE;
    double potential_min_value = DBL_MAX;
    FLOAT potential_min_pos[MAX_DIMENSION] = {bh_pos[0], bh_pos[1], bh_pos[2]};

    const int ilo = max(0, i0 - rcell);
    const int ihi = min(nx - 1, i0 + rcell);
    const int jlo = max(0, j0 - rcell);
    const int jhi = min(ny - 1, j0 + rcell);
    const int klo = max(0, k0 - rcell);
    const int khi = min(nz - 1, k0 + rcell);
    const int search_truncated_by_grid =
      (i0 - rcell < 0 || i0 + rcell >= nx ||
       j0 - rcell < 0 || j0 + rcell >= ny ||
       k0 - rcell < 0 || k0 + rcell >= nz) ? 1 : 0;
    const int search_kernel_truncated =
      (search_exceeds_ghost || search_truncated_by_grid) ? 1 : 0;
    const int active_reposition_rejected =
      (BHRepositionMethod > 0 && search_kernel_truncated) ? 1 : 0;
    const int effective_search_truncated_by_grid =
      (i0 - rcell_effective < 0 || i0 + rcell_effective >= nx ||
       j0 - rcell_effective < 0 || j0 + rcell_effective >= ny ||
       k0 - rcell_effective < 0 || k0 + rcell_effective >= nz) ? 1 : 0;
    const int effective_search_kernel_truncated =
      effective_search_truncated_by_grid ? 1 : 0;
    const int effective_reposition_rejected =
      (BHRepositionMethod > 0 &&
       effective_search_kernel_truncated) ? 1 : 0;

    /* Deterministic k-j-i traversal, matching BH seed/accretion kernels. */
    for (int k = klo; k <= khi; k++) {
      const double dz = (double(k - k0)) * double(cell_width);
      for (int j = jlo; j <= jhi; j++) {
        const double dy = (double(j - j0)) * double(cell_width);
        for (int i = ilo; i <= ihi; i++) {
          const double dx = (double(i - i0)) * double(cell_width);
          const double dist2 = dx*dx + dy*dy + dz*dz;
          if (dist2 > r2)
            continue;

          const int nindex = (k*ny + j)*nx + i;
          const int cell_is_active =
            (i >= isx && i <= iex &&
             j >= isy && j <= iey &&
             k >= isz && k <= iez) ? TRUE : FALSE;

          search_cells++;
          if (cell_is_active)
            search_active_cells++;

          const double rho_i = BaryonField[DensNum][nindex];
          if (!isfinite(rho_i))
            continue;

          FLOAT cell_pos[MAX_DIMENSION] = {0.0, 0.0, 0.0};
          cell_pos[0] = CellLeftEdge[0][0] + (FLOAT(i) + 0.5f)*cell_width;
          cell_pos[1] = CellLeftEdge[1][0] + (FLOAT(j) + 0.5f)*cell_width;
          cell_pos[2] = CellLeftEdge[2][0] + (FLOAT(k) + 0.5f)*cell_width;

          if (rho_i > diag_peak_density) {
            diag_peak_density = rho_i;
            for (int dim = 0; dim < GridRank; dim++)
              diag_peak_pos[dim] = cell_pos[dim];
            diag_peak_in_ghost = cell_is_active ? 0 : 1;
            diag_peak_from_cell = TRUE;
          } else if (rho_i == diag_peak_density && diag_peak_from_cell &&
                     BHRepositionLexicographicLess(cell_pos, diag_peak_pos, GridRank)) {
            for (int dim = 0; dim < GridRank; dim++)
              diag_peak_pos[dim] = cell_pos[dim];
            diag_peak_in_ghost = cell_is_active ? 0 : 1;
          }

          if (cell_is_active) {
            if (rho_i > active_peak_density) {
              active_peak_density = rho_i;
              for (int dim = 0; dim < GridRank; dim++)
                active_peak_pos[dim] = cell_pos[dim];
              active_peak_from_cell = TRUE;
            } else if (rho_i == active_peak_density && active_peak_from_cell &&
                       BHRepositionLexicographicLess(cell_pos, active_peak_pos, GridRank)) {
              for (int dim = 0; dim < GridRank; dim++)
                active_peak_pos[dim] = cell_pos[dim];
            }

            if (dist2 <= r2_effective) {
              if (rho_i > effective_active_peak_density) {
                effective_active_peak_density = rho_i;
                for (int dim = 0; dim < GridRank; dim++)
                  effective_active_peak_pos[dim] = cell_pos[dim];
                effective_active_peak_from_cell = TRUE;
              } else if (rho_i == effective_active_peak_density &&
                         effective_active_peak_from_cell &&
                         BHRepositionLexicographicLess(
                           cell_pos, effective_active_peak_pos, GridRank)) {
                for (int dim = 0; dim < GridRank; dim++)
                  effective_active_peak_pos[dim] = cell_pos[dim];
              }
            }
          }

          if (potential_requested && potential_available) {
            const long long pidx =
              (((long long) (k + pot_off[2]) * grav_dim[1]) +
                (long long) (j + pot_off[1])) * grav_dim[0] +
              (long long) (i + pot_off[0]);
            if (pidx >= 0 && pidx < potential_size) {
              const double phi = PotentialField[pidx];
              if (isfinite(phi)) {
                if (!potential_min_found || phi < potential_min_value) {
                  potential_min_found = TRUE;
                  potential_min_value = phi;
                  for (int dim = 0; dim < GridRank; dim++)
                    potential_min_pos[dim] = cell_pos[dim];
                } else if (phi == potential_min_value &&
                           BHRepositionLexicographicLess(cell_pos, potential_min_pos, GridRank)) {
                  for (int dim = 0; dim < GridRank; dim++)
                    potential_min_pos[dim] = cell_pos[dim];
                }
              }
            }
          }
        }
      }
    }

    const double diag_offset_code =
      sqrt((double(diag_peak_pos[0] - bh_pos[0]))*(double(diag_peak_pos[0] - bh_pos[0])) +
           (double(diag_peak_pos[1] - bh_pos[1]))*(double(diag_peak_pos[1] - bh_pos[1])) +
           (double(diag_peak_pos[2] - bh_pos[2]))*(double(diag_peak_pos[2] - bh_pos[2])));
    const double active_offset_code =
      sqrt((double(active_peak_pos[0] - bh_pos[0]))*(double(active_peak_pos[0] - bh_pos[0])) +
           (double(active_peak_pos[1] - bh_pos[1]))*(double(active_peak_pos[1] - bh_pos[1])) +
           (double(active_peak_pos[2] - bh_pos[2]))*(double(active_peak_pos[2] - bh_pos[2])));
    const double effective_active_offset_code =
      sqrt((double(effective_active_peak_pos[0] - bh_pos[0]))*(double(effective_active_peak_pos[0] - bh_pos[0])) +
           (double(effective_active_peak_pos[1] - bh_pos[1]))*(double(effective_active_peak_pos[1] - bh_pos[1])) +
           (double(effective_active_peak_pos[2] - bh_pos[2]))*(double(effective_active_peak_pos[2] - bh_pos[2])));

    const double diag_offset_kpc = diag_offset_code * kpc_per_code;
    const double active_offset_kpc = active_offset_code * kpc_per_code;
    const double effective_active_offset_kpc =
      effective_active_offset_code * kpc_per_code;

    int active_target_exists = 0;
    if (active_peak_density > bh_density)
      active_target_exists = 1;
    else if (active_peak_density == bh_density &&
             !BHRepositionPositionEqual(active_peak_pos, bh_pos, GridRank))
      active_target_exists = 1;

    int effective_active_target_exists = 0;
    if (effective_active_peak_density > bh_density)
      effective_active_target_exists = 1;
    else if (effective_active_peak_density == bh_density &&
             !BHRepositionPositionEqual(
               effective_active_peak_pos, bh_pos, GridRank))
      effective_active_target_exists = 1;

    int rejection_reason = 0;
    if (BHRepositionMethod == 0)
      rejection_reason = 5;
    else if (effective_reposition_rejected)
      rejection_reason = 2;
    else if (!effective_active_target_exists)
      rejection_reason = 4;

    double offset_from_potential_kpc = -1.0;
    if (potential_requested && potential_available && potential_min_found) {
      const double potential_offset_code =
        sqrt((double(potential_min_pos[0] - bh_pos[0]))*(double(potential_min_pos[0] - bh_pos[0])) +
             (double(potential_min_pos[1] - bh_pos[1]))*(double(potential_min_pos[1] - bh_pos[1])) +
             (double(potential_min_pos[2] - bh_pos[2]))*(double(potential_min_pos[2] - bh_pos[2])));
      offset_from_potential_kpc = potential_offset_code * kpc_per_code;
    } else if (potential_requested && !potential_available && !warned_potential_unavailable) {
      fprintf(logptr,
              "[BHREPOS_WARN] step=%d level=%d "
              "Potential field not available for repositioning diagnostics.\n",
              cycle_number, level);
      warned_potential_unavailable = TRUE;
    }

    double displacement_code = 0.0;
    int reposition_occurred = 0;
    int reposition_clamped = 0;

    FLOAT bh_pos_new[MAX_DIMENSION] = {bh_pos[0], bh_pos[1], bh_pos[2]};
    if (active_reposition_rejected) {
      bhrepos_legacy_reject_warn_count++;
      if (BHRepositionVerbose >= 2 ||
          (BHRepositionVerbose >= 1 &&
           bhrepos_legacy_reject_warn_count <=
           bhrepos_legacy_reject_warn_first_n)) {
        fprintf(logptr,
                "[BHREPOS_WARN] step=%d level=%d bh_id=%lld "
                "active_reposition_rejected=1 search_radius_over_dx=%.6g "
                "NumberOfGhostZones=%d search_truncated_by_grid=%d "
                "legacy_full_radius_diagnostic=1 effective_reposition_rejected=%d "
                "rejection_reason=%d; effective fields control particle movement.\n",
                cycle_number, level, (long long) ParticleNumber[p],
                search_radius_over_dx, NumberOfGhostZones,
                search_truncated_by_grid, effective_reposition_rejected,
                rejection_reason);
      }
    }

    if (BHRepositionMethod > 0 && effective_active_target_exists &&
        !effective_reposition_rejected) {
      if (BHRepositionMethod == 1) {
        const double max_displacement_code =
          double(BHRepositionMaxDisplacement) * double(cell_width);
        if (effective_active_offset_code <= max_displacement_code) {
          for (int dim = 0; dim < GridRank; dim++)
            bh_pos_new[dim] = effective_active_peak_pos[dim];
          displacement_code = effective_active_offset_code;
        } else if (effective_active_offset_code > 0.0 &&
                   max_displacement_code > 0.0) {
          const double frac =
            max_displacement_code / effective_active_offset_code;
          for (int dim = 0; dim < GridRank; dim++)
            bh_pos_new[dim] =
              bh_pos[dim] + FLOAT(frac * double(effective_active_peak_pos[dim] - bh_pos[dim]));
          displacement_code = max_displacement_code;
        }
      } else if (BHRepositionMethod == 2) {
        for (int dim = 0; dim < GridRank; dim++)
          bh_pos_new[dim] = effective_active_peak_pos[dim];
        displacement_code = effective_active_offset_code;
        if (BHRepositionVerbose >= 1 &&
            displacement_code > effective_active_radius_code) {
          fprintf(logptr,
                  "[BHREPOS_WARN] step=%d level=%d bh_id=%lld "
                  "teleport displacement_cells=%.8g exceeds "
                  "effective_active_radius_cells=%.8g.\n",
                  cycle_number, level, (long long) ParticleNumber[p],
                  displacement_code / cell_width,
                  effective_active_radius_cells);
        }
      }
    }

    if (displacement_code > 0.0) {
      for (int dim = 0; dim < GridRank; dim++) {
        const double left = GridLeftEdge[dim];
        const double right = GridRightEdge[dim];
        const double edge_eps =
          max(1.0e-12 * max(fabs(left), fabs(right)),
              1.0e-9 * double(cell_width));
        if (bh_pos_new[dim] < left) {
          bh_pos_new[dim] = FLOAT(left);
          reposition_clamped = 1;
        }
        if (bh_pos_new[dim] >= right) {
          bh_pos_new[dim] = FLOAT(max(left, right - edge_eps));
          reposition_clamped = 1;
        }
      }

      if (reposition_clamped && BHRepositionVerbose >= 1) {
        fprintf(logptr,
                "[BHREPOS_WARN] step=%d level=%d bh_id=%lld "
                "position clamped to active-zone bounds after reposition.\n",
                cycle_number, level, (long long) ParticleNumber[p]);
      }

      ParticlePosition[0][p] = bh_pos_new[0];
      if (GridRank > 1)
        ParticlePosition[1][p] = bh_pos_new[1];
      if (GridRank > 2)
        ParticlePosition[2][p] = bh_pos_new[2];

      displacement_code =
        sqrt((double(bh_pos_new[0] - bh_pos[0]))*(double(bh_pos_new[0] - bh_pos[0])) +
             (double(bh_pos_new[1] - bh_pos[1]))*(double(bh_pos_new[1] - bh_pos[1])) +
             (double(bh_pos_new[2] - bh_pos[2]))*(double(bh_pos_new[2] - bh_pos[2])));
      reposition_occurred = (displacement_code > 0.0) ? 1 : 0;
    }

    const double displacement_kpc = displacement_code * kpc_per_code;
    const double displacement_cells = displacement_code / double(cell_width);

    const double reposition_wall_ms = 1000.0 * (ReturnWallTime() - t0);

    bhrepos_summary_rows++;
    bhrepos_summary_moved += reposition_occurred ? 1 : 0;
    bhrepos_summary_method0_rows += (BHRepositionMethod == 0) ? 1 : 0;
    bhrepos_summary_diagnostic_only_rows +=
      (BHRepositionMethod == 0) ? 1 : 0;
    bhrepos_summary_active_radius_was_capped +=
      active_radius_was_capped ? 1 : 0;
    bhrepos_summary_active_target_exists +=
      active_target_exists ? 1 : 0;
    bhrepos_summary_effective_active_target_exists +=
      effective_active_target_exists ? 1 : 0;
    bhrepos_summary_effective_reposition_rejected +=
      effective_reposition_rejected ? 1 : 0;
    bhrepos_summary_legacy_active_reposition_rejected +=
      active_reposition_rejected ? 1 : 0;
    bhrepos_summary_search_exceeds_ghost +=
      search_exceeds_ghost ? 1 : 0;
    bhrepos_summary_search_truncated_by_grid +=
      search_truncated_by_grid ? 1 : 0;
    bhrepos_summary_effective_search_truncated_by_grid +=
      effective_search_truncated_by_grid ? 1 : 0;
    bhrepos_summary_under_resolved +=
      (search_radius_over_dx < 1.5) ? 1 : 0;
    bhrepos_summary_no_effective_target +=
      (!effective_active_target_exists) ? 1 : 0;
    bhrepos_summary_max_displacement_cells =
      max(bhrepos_summary_max_displacement_cells, displacement_cells);
    bhrepos_summary_max_effective_active_peak_offset_kpc =
      max(bhrepos_summary_max_effective_active_peak_offset_kpc,
          effective_active_offset_kpc);

    if (BHRepositionVerbose >= 2) {
      fprintf(logptr,
              "[BHREPOS] processor=%d invocation_seq=%d grid_id=%d "
              "step=%d level=%d z=%.8g bh_id=%lld "
              "bh_pos_x=%.15g bh_pos_y=%.15g bh_pos_z=%.15g "
              "bh_density=%.8e "
              "diag_peak_pos_x=%.15g diag_peak_pos_y=%.15g diag_peak_pos_z=%.15g "
              "diag_peak_density=%.8e diag_peak_in_ghost=%d diag_peak_offset_kpc=%.8g "
              "active_peak_pos_x=%.15g active_peak_pos_y=%.15g active_peak_pos_z=%.15g "
              "active_peak_density=%.8e active_peak_offset_kpc=%.8g active_target_exists=%d "
              "offset_from_potential_kpc=%.8g "
              "displacement_kpc=%.8g displacement_cells=%.8g "
              "reposition_occurred=%d reposition_clamped=%d newly_seeded_skip=%d "
              "method=%d diagnostic_only=%d "
              "search_exceeds_ghost=%d search_truncated_by_grid=%d "
              "active_reposition_rejected=%d search_radius_over_dx=%.8g "
              "under_resolved=%d "
              "search_cells=%d search_active_cells=%d reposition_wall_ms=%.4f "
              "effective_active_radius_cells=%.8g active_radius_was_capped=%d "
              "effective_active_target_exists=%d "
              "effective_active_peak_density=%.8e "
              "effective_active_peak_offset_kpc=%.8g "
              "effective_search_truncated_by_grid=%d "
              "effective_reposition_rejected=%d rejection_reason=%d\n",
              bhrepos_processor, bhrepos_invocation_seq, bhrepos_grid_id,
              cycle_number, level, zred, (long long) ParticleNumber[p],
              bh_pos[0], bh_pos[1], bh_pos[2],
              bh_density,
              diag_peak_pos[0], diag_peak_pos[1], diag_peak_pos[2],
              diag_peak_density, diag_peak_in_ghost, diag_offset_kpc,
              active_peak_pos[0], active_peak_pos[1], active_peak_pos[2],
              active_peak_density, active_offset_kpc, active_target_exists,
              offset_from_potential_kpc,
              displacement_kpc, displacement_cells,
              reposition_occurred, reposition_clamped, 0,
              BHRepositionMethod, (BHRepositionMethod == 0) ? 1 : 0,
              search_exceeds_ghost, search_truncated_by_grid,
              active_reposition_rejected, search_radius_over_dx,
              (search_radius_over_dx < 1.5) ? 1 : 0,
              search_cells, search_active_cells, reposition_wall_ms,
              effective_active_radius_cells, active_radius_was_capped,
              effective_active_target_exists,
              effective_active_peak_density,
              effective_active_offset_kpc,
              effective_search_truncated_by_grid,
              effective_reposition_rejected, rejection_reason);
      bhrepos_detail_rows_emitted++;
    }
  }

  if (BHRepositionVerbose >= 1) {
    fprintf(logptr,
            "[BHREPOS_SUMMARY] schema_version=2 "
            "processor=%d invocation_seq=%d grid_id=%d step=%d level=%d "
            "rows=%d moved=%d method0_rows=%d diagnostic_only_rows=%d "
            "active_radius_was_capped_count=%d "
            "active_target_exists_count=%d "
            "effective_active_target_exists_count=%d "
            "effective_reposition_rejected_count=%d "
            "legacy_active_reposition_rejected_count=%d "
            "search_exceeds_ghost_count=%d "
            "search_truncated_by_grid_count=%d "
            "effective_search_truncated_by_grid_count=%d "
            "under_resolved_count=%d no_effective_target_count=%d "
            "newly_seeded_count=%d "
            "max_displacement_cells=%.8g "
            "max_effective_active_peak_offset_kpc=%.8g "
            "max_search_radius_over_dx=%.8g "
            "detail_rows_emitted=%d detail_rows_complete=%d\n",
            bhrepos_processor, bhrepos_invocation_seq, bhrepos_grid_id,
            cycle_number, level,
            bhrepos_summary_rows, bhrepos_summary_moved,
            bhrepos_summary_method0_rows,
            bhrepos_summary_diagnostic_only_rows,
            bhrepos_summary_active_radius_was_capped,
            bhrepos_summary_active_target_exists,
            bhrepos_summary_effective_active_target_exists,
            bhrepos_summary_effective_reposition_rejected,
            bhrepos_summary_legacy_active_reposition_rejected,
            bhrepos_summary_search_exceeds_ghost,
            bhrepos_summary_search_truncated_by_grid,
            bhrepos_summary_effective_search_truncated_by_grid,
            bhrepos_summary_under_resolved,
            bhrepos_summary_no_effective_target,
            bhrepos_summary_newly_seeded,
            bhrepos_summary_max_displacement_cells,
            bhrepos_summary_max_effective_active_peak_offset_kpc,
            bhrepos_summary_max_search_radius_over_dx,
            bhrepos_detail_rows_emitted,
            bhrepos_detail_rows_complete);
  }

  return SUCCESS;
}
