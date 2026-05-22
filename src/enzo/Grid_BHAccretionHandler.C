/***********************************************************************
/
/  GRID CLASS (BH ACCRETION DIAGNOSTICS + ACTIVE ACCRETION, PHASE B)
/
/  PURPOSE:
/    Evaluate per-BH local gas diagnostics (Phase A) and, in Phase B,
/    convert capped rates into actual gas removal and BH mass growth.
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
#include "BHAccretionSafety.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
double ReturnWallTime();

struct BHAccretionChannelAccumulator {
  double mass;
  double rho_sum;
  double cs_sum;
  double vx_sum;
  double vy_sum;
  double vz_sum;
  double t_sum;
  double z_sum;
  double lx;
  double ly;
  double lz;
  int ncells;
};

static int BHAccretionIsBHType(int type)
{
  return (type == PARTICLE_TYPE_MBH || type == PARTICLE_TYPE_BLACK_HOLE);
}

static int BHAccretionIsNewlySeededThisPass(PINT particle_id)
{
  /*
   * Newly seeded BH particles are created with INT_UNDEFINED ids and receive
   * final unique ids later in the star particle finalize stage. Skipping
   * INT_UNDEFINED here gives a robust same-pass exclusion without float-time
   * comparisons.
   */
  return (particle_id == INT_UNDEFINED);
}

static void BHAccretionResetRealizedMdot(int NumberOfParticles,
                                         int *ParticleType,
                                         float *ParticleAttribute[],
                                         int NumberOfParticleAttributes)
{
  if (NumberOfParticleAttributes <= PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED ||
      ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED] == NULL)
    return;

  for (int p = 0; p < NumberOfParticles; p++)
    if (BHAccretionIsBHType(ParticleType[p]))
      ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED][p] = 0.0f;
}

struct BHAccretionParticleOrder {
  PINT *ParticleIDs;
  BHAccretionParticleOrder(PINT *ids) : ParticleIDs(ids) { }
  bool operator()(int a, int b) const {
    if (ParticleIDs[a] < ParticleIDs[b])
      return true;
    if (ParticleIDs[a] > ParticleIDs[b])
      return false;
    return (a < b);
  }
};

static int BHAccretionPositionCoveredByChild(HierarchyEntry *SubgridPointer,
                                             FLOAT *pos)
{
  for (HierarchyEntry *sub = SubgridPointer; sub != NULL;
       sub = sub->NextGridThisLevel)
    if (sub->GridData != NULL && sub->GridData->PointInGrid(pos))
      return TRUE;

  return FALSE;
}

static int BHAccretionNonAuthoritativeMassRatio(
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

static double BHAccretionKernelRadiusCode(float KernelRadiusPhysKpc,
                                          FLOAT time,
                                          float LengthUnits,
                                          float *a_phys_out,
                                          float *h_out,
                                          float *box_kpch_out)
{
  float h_param = 1.0f;
  float a_phys = 1.0f;
  float box_kpch = 0.0f;
  double kernel_radius_code = 0.0;

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
    const float kernel_com_kpch = KernelRadiusPhysKpc * h_param / a_phys;
    if (box_kpch > 0.0f)
      kernel_radius_code = kernel_com_kpch / box_kpch;
  } else {
    h_param = 1.0f;
    a_phys = 1.0f;
    box_kpch = LengthUnits / kpc_cm;
    if (box_kpch > 0.0f)
      kernel_radius_code = KernelRadiusPhysKpc / box_kpch;
  }

  if (a_phys_out != NULL)
    *a_phys_out = a_phys;
  if (h_out != NULL)
    *h_out = h_param;
  if (box_kpch_out != NULL)
    *box_kpch_out = box_kpch;

  return kernel_radius_code;
}

int grid::BHAccretionDiagnosticHandler(HierarchyEntry* SubgridPointer,
                                       int level, int cycle_number,
                                       float dtLevelAbove)
{
  (void) dtLevelAbove;

  if (!BHAccretionMethod)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0 || NumberOfParticles <= 0)
    return SUCCESS;

  BHAccretionResetRealizedMdot(NumberOfParticles, ParticleType,
                               ParticleAttribute, NumberOfParticleAttributes);

  /* Ensure the transient SF mask exists for accretion cell blocking. */
  if (BaryonField[NumberOfBaryonFields] == NULL) {
    this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
    for (HierarchyEntry *sg = SubgridPointer; sg != NULL; sg = sg->NextGridThisLevel)
      this->ZeroSolutionUnderSubgrid(sg->GridData, ZERO_UNDER_SUBGRID_FIELD);
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL)
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  std::vector<float> temperature(size, 0.0f);
  if (this->ComputeTemperatureField(&temperature[0]) == FAIL)
    ENZO_FAIL("Error in grid->ComputeTemperatureField for BH accretion.");

  const int cooling_available = (RadiativeCooling != 0);
  std::vector<float> cooling_time(size, 0.0f);
  if (cooling_available)
    if (this->ComputeCoolingTime(&cooling_time[0]) == FAIL)
      ENZO_FAIL("Error in grid->ComputeCoolingTime for BH accretion.");

  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum;
  int MetalIaNum, MetalIINum;
  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
                                 MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.");

  const int metal_field_present = (MetalNum != -1 || SNColourNum != -1);
  std::vector<float> metal_fraction(size, 0.0f);
  if (metal_field_present) {
    for (int n = 0; n < size; n++) {
      const float den = BaryonField[DensNum][n];
      if (den <= 0.0f)
        continue;
      if (MetalNum != -1 && SNColourNum != -1)
        metal_fraction[n] = (BaryonField[MetalNum][n] + BaryonField[SNColourNum][n]) / den;
      else if (MetalNum != -1)
        metal_fraction[n] = BaryonField[MetalNum][n] / den;
      else
        metal_fraction[n] = BaryonField[SNColourNum][n] / den;
    }
  }

  float DensityUnits = 1.0f, LengthUnits = 1.0f, TemperatureUnits = 1.0f;
  float TimeUnits = 1.0f, VelocityUnits = 1.0f;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL)
    ENZO_FAIL("Error in GetUnits.");

  const float cell_width = float(CellWidth[0][0]);
  if (cell_width <= 0.0f)
    return SUCCESS;

  const double kernel_radius_code =
    BHAccretionKernelRadiusCode(BHAccretionKernelRadius, Time, LengthUnits,
                                NULL, NULL, NULL);
  if (kernel_radius_code <= 0.0)
    return SUCCESS;

  FLOAT a = 1.0, dadt = 0.0;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  const double zred = ComovingCoordinates ?
    double((1.0 + InitialRedshift) / a - 1.0) : 0.0;

  const double cell_volume_code = pow(cell_width, 3.0);
  const double mass_units = DensityUnits * pow(LengthUnits, 3.0);
  const double mass_to_msun = mass_units / SolarMass;
  const double mdot_to_msunyr = mass_units * yr_s / (SolarMass * TimeUnits);
  const double vel_to_cgs = VelocityUnits;
  const double rho_to_cgs = DensityUnits;

  const float x_h = (CoolData.HydrogenFractionByMass > 0.0f) ?
    CoolData.HydrogenFractionByMass : 0.76f;
  const double rho_to_nh = rho_to_cgs * x_h / mh;
  const double g_code = GravConst * DensityUnits * TimeUnits * TimeUnits;

  const int nx = GridDimension[0];
  const int ny = GridDimension[1];
  const int nz = GridDimension[2];
  const int xo = 1;
  const int yo = nx;
  const int zo = nx*ny;

  const int rcell = max(0, int(ceil(kernel_radius_code / cell_width)));
  const double r2 = kernel_radius_code * kernel_radius_code;

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
    if (BHAccretionIsBHType(ParticleType[p]))
      bh_particles.push_back(p);

  std::sort(bh_particles.begin(), bh_particles.end(),
            BHAccretionParticleOrder(ParticleNumber));

  for (int bp = 0; bp < int(bh_particles.size()); bp++) {
    const int p = bh_particles[bp];

    if (BHAccretionIsNewlySeededThisPass(ParticleNumber[p]))
      continue;

    FLOAT bh_pos[MAX_DIMENSION];
    bh_pos[0] = ParticlePosition[0][p];
    bh_pos[1] = (GridRank > 1) ? ParticlePosition[1][p] : 0.0;
    bh_pos[2] = (GridRank > 2) ? ParticlePosition[2][p] : 0.0;

    double bh_mass_code = ParticleMass[p];

    const int child_covered =
      (SubgridPointer != NULL &&
       BHAccretionPositionCoveredByChild(SubgridPointer, bh_pos)) ? 1 : 0;
    if (child_covered) {
      if (BHAccretionVerbose >= 1) {
        double expected_full_mass = -1.0;
        double ratio = -1.0;
        BHAccretionNonAuthoritativeMassRatio(
          NumberOfParticleAttributes, ParticleAttribute, p, bh_mass_code,
          &expected_full_mass, &ratio);
        fprintf(logptr,
                "[BHACCR_SKIP_NON_AUTHORITATIVE] step=%d level=%d grid_id=%d "
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

    const int point_in_grid = this->PointInGrid(bh_pos) ? 1 : 0;
    if (!point_in_grid)
      continue;

    double expected_full_mass = -1.0;
    double mass_attribute_ratio = -1.0;
    if (BHAccretionNonAuthoritativeMassRatio(
          NumberOfParticleAttributes, ParticleAttribute, p, bh_mass_code,
          &expected_full_mass, &mass_attribute_ratio)) {
      if (BHAccretionVerbose >= 1) {
        fprintf(logptr,
                "[BHACCR_SKIP_NON_AUTHORITATIVE] step=%d level=%d grid_id=%d "
                "bh_id=%lld particle_index=%d ParticleMass=%.15e "
                "expected_full_mass=%.15e ratio=%.15e reason=mass_attribute_ratio\n",
                cycle_number, level, this->GetGridID(),
                (long long) ParticleNumber[p], p, bh_mass_code,
                expected_full_mass, mass_attribute_ratio);
      }
      continue;
    }

    const double reconcile_max_rel = 1.0e-4;
    const int has_mass_attrs =
      (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BH_FORMATION_MASS &&
       NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS &&
       ParticleAttribute[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS] != NULL &&
       ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS] != NULL) ? 1 : 0;
    const double particle_mass_before = double(ParticleMass[p]);
    const double formation_mass_code = has_mass_attrs ?
      double(ParticleAttribute[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS][p]) : -1.0;
    const double accreted_mass_code = has_mass_attrs ?
      double(ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS][p]) : -1.0;
    const double expected_mass = has_mass_attrs ?
      formation_mass_code + accreted_mass_code : -1.0;
    const double residual =
      (has_mass_attrs && isfinite(particle_mass_before) && isfinite(expected_mass)) ?
      particle_mass_before - expected_mass : -1.0;
    const double abs_residual = isfinite(residual) ? fabs(residual) : -1.0;
    const double mass_scale =
      (has_mass_attrs && isfinite(particle_mass_before) && isfinite(expected_mass)) ?
      max(fabs(particle_mass_before), fabs(expected_mass)) : -1.0;
    const double mass_sync_tolerance =
      (isfinite(mass_scale) && mass_scale > 0.0) ? 1.0e-8 * mass_scale : 1.0e-30;
    const double rel_residual =
      (isfinite(abs_residual) && isfinite(mass_scale) && mass_scale > 0.0) ?
      abs_residual / mass_scale : -1.0;
    const int mass_bookkeeping_valid =
      (has_mass_attrs &&
       isfinite(particle_mass_before) && particle_mass_before > 0.0 &&
       isfinite(formation_mass_code) && formation_mass_code > 0.0 &&
       isfinite(accreted_mass_code) && accreted_mass_code >= 0.0 &&
       isfinite(expected_mass) && expected_mass > 0.0 &&
       isfinite(residual) && isfinite(abs_residual) &&
       isfinite(rel_residual)) ? 1 : 0;

    if (!mass_bookkeeping_valid ||
        abs_residual > mass_sync_tolerance) {
      const int small_upward_mismatch =
        (mass_bookkeeping_valid &&
         particle_mass_before < expected_mass &&
         rel_residual <= reconcile_max_rel) ? 1 : 0;

      if (small_upward_mismatch) {
        ParticleMass[p] = float(expected_mass);
        bh_mass_code = expected_mass;

        fprintf(logptr,
                "[BHACCR_MASS_RECONCILE] rank=%d processor=%d step=%d level=%d "
                "grid_id=%d particle_index=%d bh_id=%lld "
                "ParticleMass_before=%.15e BHFormationMass=%.15e "
                "BHAccretedMass=%.15e expected_mass=%.15e residual=%.15e "
                "abs_residual=%.15e rel_residual=%.15e tolerance=%.15e "
                "reconcile_max_rel=%.15e ParticleMass_after=%.15e "
                "bh_mass_code_after=%.15e Time=%.15e "
                "action=particle_mass_set_to_expected\n",
                MyProcessorNumber, ProcessorNumber, cycle_number, level,
                this->GetGridID(), p, (long long) ParticleNumber[p],
                particle_mass_before, formation_mass_code, accreted_mass_code,
                expected_mass, residual, abs_residual, rel_residual,
                mass_sync_tolerance, reconcile_max_rel, double(ParticleMass[p]),
                bh_mass_code, double(Time));
        fflush(logptr);
      } else {
        const double precheck_dt_code = max(0.0, double(this->dtFixed));
        fprintf(stderr,
                "[BHACCR_MASS_MISMATCH_PRE] rank=%d processor=%d step=%d level=%d "
                "grid_id=%d particle_index=%d bh_id=%lld particle_type=%d "
                "NumberOfParticleAttributes=%d ParticleMass_before=%.15e "
                "BHFormationMass=%.15e BHAccretedMass=%.15e "
                "expected_mass=%.15e residual=%.15e abs_residual=%.15e "
                "rel_residual=%.15e tolerance=%.15e "
                "reconcile_max_rel=%.15e child_covered=%d point_in_grid=%d "
                "mass_attribute_ratio=%.15e authority_status=passed_authority_gates "
                "Time=%.15e dt_code=%.15e action=fail_before_gas_removal\n",
                MyProcessorNumber, ProcessorNumber, cycle_number, level,
                this->GetGridID(), p, (long long) ParticleNumber[p],
                ParticleType[p], NumberOfParticleAttributes,
                particle_mass_before, formation_mass_code, accreted_mass_code,
                expected_mass, residual, abs_residual, rel_residual,
                mass_sync_tolerance, reconcile_max_rel, child_covered,
                point_in_grid, mass_attribute_ratio, double(Time),
                precheck_dt_code);
        fflush(stderr);
        ENZO_FAIL("BHAccretion pre-update invariant failed: BH mass != BHFormationMass + BHAccretedMass.");
      }
    }

    const double t0_all = ReturnWallTime();

    const double bh_mass_msun = bh_mass_code * mass_to_msun;
    const double bh_vel_x = ParticleVelocity[0][p];
    const double bh_vel_y = (GridRank > 1) ? ParticleVelocity[1][p] : 0.0;
    const double bh_vel_z = (GridRank > 2) ? ParticleVelocity[2][p] : 0.0;

    int i0 = int(floor((ParticlePosition[0][p] - CellLeftEdge[0][0]) / cell_width));
    int j0 = int(floor((ParticlePosition[1][p] - CellLeftEdge[1][0]) / cell_width));
    int k0 = int(floor((ParticlePosition[2][p] - CellLeftEdge[2][0]) / cell_width));
    i0 = max(0, min(nx - 1, i0));
    j0 = max(0, min(ny - 1, j0));
    k0 = max(0, min(nz - 1, k0));

	    const int ilo = max(0, i0 - rcell);
	    const int ihi = min(nx - 1, i0 + rcell);
	    const int jlo = max(0, j0 - rcell);
	    const int jhi = min(ny - 1, j0 + rcell);
	    const int klo = max(0, k0 - rcell);
	    const int khi = min(nz - 1, k0 + rcell);

	    const double kernel_cells_per_radius = kernel_radius_code / cell_width;
	    const int kernel_exceeds_ghost =
	      (kernel_cells_per_radius > double(NumberOfGhostZones)) ? 1 : 0;
	    const int kernel_truncated_by_grid =
	      (i0 - rcell < 0 || i0 + rcell >= nx ||
	       j0 - rcell < 0 || j0 + rcell >= ny ||
	       k0 - rcell < 0 || k0 + rcell >= nz) ? 1 : 0;
	    const int kernel_truncated =
	      (kernel_exceeds_ghost || kernel_truncated_by_grid) ? 1 : 0;

	    long long n_kernel_cells_requested = 0;
	    if (rcell <= 4096) {
	      for (int kk = -rcell; kk <= rcell; kk++) {
	        const double dz = double(kk) * double(cell_width);
	        for (int jj = -rcell; jj <= rcell; jj++) {
	          const double dy = double(jj) * double(cell_width);
	          const double rem = r2 - dy*dy - dz*dz;
	          if (rem < 0.0)
	            continue;
	          const long long imax = (long long) floor(sqrt(rem) / cell_width);
	          n_kernel_cells_requested += 2*imax + 1;
	        }
	      }
	    } else {
	      n_kernel_cells_requested =
	        (long long) ((4.0*pi/3.0) * pow(double(rcell), 3.0));
	    }

	    BHAccretionChannelAccumulator hot = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	                                         0.0, 0.0, 0.0, 0};
	    BHAccretionChannelAccumulator cold = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	                                          0.0, 0.0, 0.0, 0};

	    int n_kernel_cells = 0;
	    int n_kernel_active_cells = 0;
	    int n_hot_cells = 0;
	    int n_cold_cells = 0;
	    int n_fallback_cells = 0;

    /* Deterministic k-j-i traversal, matching BH seeding kernels. */
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
          const double rho = BaryonField[DensNum][nindex];
          if (!isfinite(rho) || rho <= 0.0)
            continue;

          const double mcell = rho * cell_volume_code;
	          if (mcell <= 0.0)
	            continue;

	          n_kernel_cells++;
	          if (i >= isx && i <= iex && j >= isy && j <= iey &&
	              k >= isz && k <= iez)
	            n_kernel_active_cells++;

	          double vx = BaryonField[Vel1Num][nindex];
          double vy = BaryonField[Vel2Num][nindex];
          double vz = BaryonField[Vel3Num][nindex];
          if (HydroMethod == Zeus_Hydro) {
            if (i + 1 < nx)
              vx = 0.5*(BaryonField[Vel1Num][nindex] + BaryonField[Vel1Num][nindex+xo]);
            if (j + 1 < ny)
              vy = 0.5*(BaryonField[Vel2Num][nindex] + BaryonField[Vel2Num][nindex+yo]);
            if (k + 1 < nz)
              vz = 0.5*(BaryonField[Vel3Num][nindex] + BaryonField[Vel3Num][nindex+zo]);
          }

          const double T = temperature[nindex];
          double cs = 0.0;
          if (isfinite(T) && T > 0.0 && Mu > 0.0 && VelocityUnits > 0.0)
            cs = sqrt(Gamma * kboltz * T / (Mu * mh)) / VelocityUnits;

          const double t_dyn = (g_code > 0.0) ? 1.0 / sqrt(g_code * rho) : -1.0;
          const double t_cool = cooling_available ? cooling_time[nindex] : -1.0;

          int cold_cell = FALSE;
          if (cooling_available && isfinite(t_cool) && t_cool > 0.0 &&
              isfinite(t_dyn) && t_dyn > 0.0) {
            cold_cell = ((t_cool / t_dyn) < 1.0) ? TRUE : FALSE;
          } else {
            n_fallback_cells++;
            cold_cell = (isfinite(T) && T < BHAccretionTSplitFloor) ? TRUE : FALSE;
          }

          BHAccretionChannelAccumulator *acc = cold_cell ? &cold : &hot;
          acc->mass += mcell;
          acc->rho_sum += mcell * rho;
          acc->cs_sum += mcell * cs;
          acc->vx_sum += mcell * vx;
          acc->vy_sum += mcell * vy;
          acc->vz_sum += mcell * vz;
          acc->t_sum += mcell * T;
          acc->z_sum += mcell * metal_fraction[nindex];
          acc->ncells++;

          if (cold_cell) {
            const double xcell = CellLeftEdge[0][0] + (double(i) + 0.5)*cell_width;
            const double ycell = CellLeftEdge[1][0] + (double(j) + 0.5)*cell_width;
            const double zcell = CellLeftEdge[2][0] + (double(k) + 0.5)*cell_width;
            const double rx = xcell - ParticlePosition[0][p];
            const double ry = ycell - ParticlePosition[1][p];
            const double rz = zcell - ParticlePosition[2][p];
            const double dvx = vx - bh_vel_x;
            const double dvy = vy - bh_vel_y;
            const double dvz = vz - bh_vel_z;
            acc->lx += mcell * (ry*dvz - rz*dvy);
            acc->ly += mcell * (rz*dvx - rx*dvz);
            acc->lz += mcell * (rx*dvy - ry*dvx);
            n_cold_cells++;
          } else {
            n_hot_cells++;
          }
        }
      }
    }

    const double total_mass = hot.mass + cold.mass;
    const double f_hot = (total_mass > 0.0) ? hot.mass / total_mass : 0.0;
    const double f_cold = (total_mass > 0.0) ? cold.mass / total_mass : 0.0;

    double hot_rho_avg = 0.0, hot_cs_avg = 0.0, hot_vx_avg = 0.0, hot_vy_avg = 0.0;
    double hot_vz_avg = 0.0, hot_t_avg = 0.0;
    if (hot.mass > 0.0) {
      hot_rho_avg = hot.rho_sum / hot.mass;
      hot_cs_avg = hot.cs_sum / hot.mass;
      hot_vx_avg = hot.vx_sum / hot.mass;
      hot_vy_avg = hot.vy_sum / hot.mass;
      hot_vz_avg = hot.vz_sum / hot.mass;
      hot_t_avg = hot.t_sum / hot.mass;
    }

    double cold_rho_avg = 0.0, cold_cs_avg = 0.0, cold_vx_avg = 0.0, cold_vy_avg = 0.0;
    double cold_vz_avg = 0.0, cold_t_avg = 0.0, v_rot_cold = 0.0;
    if (cold.mass > 0.0) {
      cold_rho_avg = cold.rho_sum / cold.mass;
      cold_cs_avg = cold.cs_sum / cold.mass;
      cold_vx_avg = cold.vx_sum / cold.mass;
      cold_vy_avg = cold.vy_sum / cold.mass;
      cold_vz_avg = cold.vz_sum / cold.mass;
      cold_t_avg = cold.t_sum / cold.mass;
      const double lx_spec = cold.lx / cold.mass;
      const double ly_spec = cold.ly / cold.mass;
      const double lz_spec = cold.lz / cold.mass;
      if (kernel_radius_code > 0.0)
        v_rot_cold = sqrt(lx_spec*lx_spec + ly_spec*ly_spec + lz_spec*lz_spec) /
                     kernel_radius_code;
    }

    const double v_rel_hot = (hot.mass > 0.0) ?
      sqrt((bh_vel_x-hot_vx_avg)*(bh_vel_x-hot_vx_avg) +
           (bh_vel_y-hot_vy_avg)*(bh_vel_y-hot_vy_avg) +
           (bh_vel_z-hot_vz_avg)*(bh_vel_z-hot_vz_avg)) : 0.0;

    const double v_rel_cold = (cold.mass > 0.0) ?
      sqrt((bh_vel_x-cold_vx_avg)*(bh_vel_x-cold_vx_avg) +
           (bh_vel_y-cold_vy_avg)*(bh_vel_y-cold_vy_avg) +
           (bh_vel_z-cold_vz_avg)*(bh_vel_z-cold_vz_avg)) : 0.0;

    double alpha_boost = 0.0;
    double mdot_bondi_hot = 0.0;
    double mdot_hot_raw = 0.0;
    if (hot.mass > 0.0) {
      const double hot_denom = pow(hot_cs_avg*hot_cs_avg + v_rel_hot*v_rel_hot, 1.5);
      if (hot_denom > 0.0)
        mdot_bondi_hot = 4.0*pi*g_code*g_code*bh_mass_code*bh_mass_code*hot_rho_avg /
                         hot_denom;
      const double n_h_hot = hot_rho_avg * rho_to_nh;
      double alpha = 1.0;
      if (BHAccretionNHStar > 0.0f && n_h_hot > 0.0)
        alpha = pow(n_h_hot / BHAccretionNHStar, BHAccretionBeta);
      alpha = max(1.0, alpha);
      alpha = min(double(BHAccretionAlphaMax), alpha);
      alpha_boost = alpha;
      mdot_hot_raw = alpha_boost * mdot_bondi_hot;
      if (!isfinite(mdot_hot_raw) || mdot_hot_raw < 0.0)
        mdot_hot_raw = 0.0;
    }

    double f_am = 0.0;
    double mdot_bondi_cold = 0.0;
    double mdot_cold_raw = 0.0;
    if (cold.mass > 0.0) {
      const double cold_denom = pow(cold_cs_avg*cold_cs_avg + v_rel_cold*v_rel_cold, 1.5);
      if (cold_denom > 0.0)
        mdot_bondi_cold = 4.0*pi*g_code*g_code*bh_mass_code*bh_mass_code*cold_rho_avg /
                          cold_denom;
      if (v_rot_cold > 0.0)
        f_am = min(1.0, pow(cold_cs_avg / v_rot_cold, 3.0) / double(BHAccretionCVisc));
      else
        f_am = 1.0;
      if (f_am < 0.0)
        f_am = 0.0;
      mdot_cold_raw = f_am * mdot_bondi_cold;
      if (!isfinite(mdot_cold_raw) || mdot_cold_raw < 0.0)
        mdot_cold_raw = 0.0;
    }

    const double mdot_total_raw = mdot_hot_raw + mdot_cold_raw;
    const double bh_mass_cgs = bh_mass_code * mass_units;
    double mdot_edd_code = 0.0;
    if (BHAccretionRadiativeEfficiency > 0.0f && bh_mass_cgs > 0.0) {
      const double mdot_edd_cgs =
        4.0*pi*GravConst*bh_mass_cgs*mh /
        (BHAccretionRadiativeEfficiency * sigma_thompson * clight);
      mdot_edd_code = mdot_edd_cgs * TimeUnits / mass_units;
    }

    const double lambda_edd =
      (mdot_edd_code > 0.0) ? mdot_total_raw / mdot_edd_code : 0.0;
    const double f_edd = lambda_edd;

    double frac_cap = 1.0;
    if (mdot_edd_code > 0.0 && mdot_total_raw > 0.0) {
      if (BHFeedbackEddingtonFactor == 1.0f) {
        /* Exact Phase B path: no Eddington-factor multiply at default. */
        frac_cap = min(1.0, mdot_edd_code / mdot_total_raw);
      } else {
        const double mdot_edd_effective =
          double(BHFeedbackEddingtonFactor) * mdot_edd_code;
        frac_cap = min(1.0, mdot_edd_effective / mdot_total_raw);
      }
    }
    frac_cap = min(1.0, max(0.0, frac_cap));
    const int cap_active = (frac_cap < 1.0) ? 1 : 0;

    const double mdot_hot_capped = mdot_hot_raw * frac_cap;
    const double mdot_cold_capped = mdot_cold_raw * frac_cap;
    const double mdot_total_capped = mdot_hot_capped + mdot_cold_capped;

    /* Preserve Phase A field semantics (capped diagnostic split). */
    const double mdot_actual = mdot_total_capped;
    const double mdot_hot_actual = mdot_hot_capped;
    const double mdot_cold_actual = mdot_cold_capped;

    const double t_diag_done = ReturnWallTime();
    const double accretion_diag_wall_ms = 1000.0 * (t_diag_done - t0_all);

    const double dt_code = max(0.0, double(this->dtFixed));
    double dm_requested = mdot_total_capped * dt_code;
    if (!isfinite(dm_requested) || dm_requested < 0.0)
      dm_requested = 0.0;

    const int removal_mode = (BHAccretionRemovalMode == 1) ? 1 : 0;
    const int removal_radius_cells = max(0, BHAccretionRemovalRadius);
    const double removal_radius_code =
      (removal_mode == 1) ? 0.0 : double(removal_radius_cells) * double(cell_width);
    const int rcell_remove = (removal_mode == 1) ? 0 :
      max(0, int(ceil(removal_radius_code / cell_width)));
    const double remove_r2 = removal_radius_code * removal_radius_code;

	    std::vector<int> remove_index;
	    std::vector<double> remove_mass;
	    std::vector<double> remove_safe_mass;
	    std::vector<double> remove_velx;
	    std::vector<double> remove_vely;
	    std::vector<double> remove_velz;

    const int reserve_n = max(1, (2*rcell_remove+1)*(2*rcell_remove+1)*(2*rcell_remove+1));
	    remove_index.reserve(reserve_n);
	    remove_mass.reserve(reserve_n);
	    remove_safe_mass.reserve(reserve_n);
	    remove_velx.reserve(reserve_n);
	    remove_vely.reserve(reserve_n);
	    remove_velz.reserve(reserve_n);

	    double removal_kernel_gas = 0.0;
	    double removal_kernel_safe_gas = 0.0;
	    const double density_floor_code =
	      BHAccretionDensityFloorCode(tiny_number, SmallRho);
	    const double density_floor_mass = density_floor_code * cell_volume_code;
	    const int accretion_source_rejected = kernel_truncated;
	    int n_floor_rejected_cells = 0;

	    if (!accretion_source_rejected && removal_mode == 1) {
	      if (i0 >= isx && i0 <= iex && j0 >= isy && j0 <= iey && k0 >= isz && k0 <= iez) {
	        const int nindex = (k0*ny + j0)*nx + i0;
	        const double rho = BaryonField[DensNum][nindex];
	        if (rho > 0.0) {
	          const double mcell = rho * cell_volume_code;
	          if (mcell > 0.0) {
	            const double safe_mass =
	              BHAccretionSafeRemovableMass(mcell, density_floor_mass);
	            removal_kernel_gas += mcell;
	            removal_kernel_safe_gas += safe_mass;
	            if (safe_mass <= 0.0) {
	              n_floor_rejected_cells++;
	            } else {
	              double vx = BaryonField[Vel1Num][nindex];
	              double vy = BaryonField[Vel2Num][nindex];
	              double vz = BaryonField[Vel3Num][nindex];
	              if (HydroMethod == Zeus_Hydro) {
	                if (i0 + 1 < nx)
	                  vx = 0.5*(BaryonField[Vel1Num][nindex] + BaryonField[Vel1Num][nindex+xo]);
	                if (j0 + 1 < ny)
	                  vy = 0.5*(BaryonField[Vel2Num][nindex] + BaryonField[Vel2Num][nindex+yo]);
	                if (k0 + 1 < nz)
	                  vz = 0.5*(BaryonField[Vel3Num][nindex] + BaryonField[Vel3Num][nindex+zo]);
	              }
	              remove_index.push_back(nindex);
	              remove_mass.push_back(mcell);
	              remove_safe_mass.push_back(safe_mass);
	              remove_velx.push_back(vx);
	              remove_vely.push_back(vy);
	              remove_velz.push_back(vz);
	            }
	          }
	        }
	      }
	    } else if (!accretion_source_rejected) {
      const int ilo_r = max(0, i0 - rcell_remove);
      const int ihi_r = min(nx - 1, i0 + rcell_remove);
      const int jlo_r = max(0, j0 - rcell_remove);
      const int jhi_r = min(ny - 1, j0 + rcell_remove);
      const int klo_r = max(0, k0 - rcell_remove);
      const int khi_r = min(nz - 1, k0 + rcell_remove);

      for (int k = klo_r; k <= khi_r; k++) {
        const double dz = (double(k - k0)) * double(cell_width);
        for (int j = jlo_r; j <= jhi_r; j++) {
          const double dy = (double(j - j0)) * double(cell_width);
          for (int i = ilo_r; i <= ihi_r; i++) {
            const double dx = (double(i - i0)) * double(cell_width);
            const double dist2 = dx*dx + dy*dy + dz*dz;
            if (dist2 > remove_r2)
              continue;

            if (i < isx || i > iex || j < isy || j > iey || k < isz || k > iez)
              continue;

            const int nindex = (k*ny + j)*nx + i;
            const double rho = BaryonField[DensNum][nindex];
            if (rho <= 0.0)
              continue;

	            const double mcell = rho * cell_volume_code;
	            if (mcell <= 0.0)
	              continue;
	            const double safe_mass =
	              BHAccretionSafeRemovableMass(mcell, density_floor_mass);
	            removal_kernel_gas += mcell;
	            removal_kernel_safe_gas += safe_mass;
	            if (safe_mass <= 0.0) {
	              n_floor_rejected_cells++;
	              continue;
	            }

	            double vx = BaryonField[Vel1Num][nindex];
            double vy = BaryonField[Vel2Num][nindex];
            double vz = BaryonField[Vel3Num][nindex];
            if (HydroMethod == Zeus_Hydro) {
              if (i + 1 < nx)
                vx = 0.5*(BaryonField[Vel1Num][nindex] + BaryonField[Vel1Num][nindex+xo]);
              if (j + 1 < ny)
                vy = 0.5*(BaryonField[Vel2Num][nindex] + BaryonField[Vel2Num][nindex+yo]);
              if (k + 1 < nz)
                vz = 0.5*(BaryonField[Vel3Num][nindex] + BaryonField[Vel3Num][nindex+zo]);
            }

	            remove_index.push_back(nindex);
	            remove_mass.push_back(mcell);
	            remove_safe_mass.push_back(safe_mass);
	            remove_velx.push_back(vx);
	            remove_vely.push_back(vy);
	            remove_velz.push_back(vz);
	          }
	        }
	      }
    }

	    double dm_removed_total = 0.0;
	    int removal_gas_limited = 0;
	    double frac_gas = 1.0;
	    const double availability_tol =
	      1.0e-12 * max(1.0, max(removal_kernel_safe_gas, dm_requested));
	    if (dm_requested > 0.0) {
	      if (accretion_source_rejected) {
	        dm_removed_total = 0.0;
	        removal_gas_limited = 1;
	      } else if (removal_kernel_safe_gas + availability_tol >= dm_requested) {
	        dm_removed_total = dm_requested;
	      } else {
	        dm_removed_total = max(0.0, removal_kernel_safe_gas);
	        removal_gas_limited = 1;
	      }

      frac_gas = dm_removed_total / dm_requested;
      frac_gas = min(1.0, max(0.0, frac_gas));
    }

    std::vector<double> dm_removed_cell(remove_index.size(), 0.0);

    double removed_mass = 0.0;
    double delta_px = 0.0;
    double delta_py = 0.0;
    double delta_pz = 0.0;
    int removal_cells = 0;
    int n_sf_blocked_cells = 0;

	    double min_density_before = DBL_MAX;
	    double min_density_after = DBL_MAX;
	    double min_eint_before = DBL_MAX;
	    double min_eint_after = DBL_MAX;
	    double min_pressure_before = DBL_MAX;
	    double min_pressure_after = DBL_MAX;
	    double max_cell_removal_fraction = 0.0;
	    int n_cells_clamped = 0;
	    int n_energy_consistency_fail = 0;

	    if (dm_removed_total > 0.0 && !remove_index.empty() && removal_kernel_safe_gas > 0.0) {
	      double remaining_to_remove = dm_removed_total;
	      double remaining_available = removal_kernel_safe_gas;

	      for (int n = 0; n < int(remove_index.size()); n++) {
	        const double mcell = remove_safe_mass[n];

	        double dm = 0.0;
	        if (n == int(remove_index.size()) - 1)
          dm = remaining_to_remove;
        else if (remaining_available > 0.0)
          dm = remaining_to_remove * (mcell / remaining_available);

	        dm = max(0.0, dm);
	        dm = min(dm, mcell);

        dm_removed_cell[n] = dm;
        remaining_to_remove -= dm;
	        remaining_available -= mcell;
	      }

      if (remaining_to_remove != 0.0) {
        if (remaining_to_remove > 0.0) {
          for (int n = int(remove_index.size()) - 1;
               n >= 0 && remaining_to_remove > 0.0; n--) {
	            const double cap = remove_safe_mass[n] - dm_removed_cell[n];
            if (cap <= 0.0)
              continue;
            const double add = min(cap, remaining_to_remove);
            dm_removed_cell[n] += add;
            remaining_to_remove -= add;
          }
        } else {
          for (int n = int(remove_index.size()) - 1;
               n >= 0 && remaining_to_remove < 0.0; n--) {
            if (dm_removed_cell[n] <= 0.0)
              continue;
            const double sub = min(dm_removed_cell[n], -remaining_to_remove);
            dm_removed_cell[n] -= sub;
            remaining_to_remove += sub;
          }
        }
      }

	      for (int n = 0; n < int(remove_index.size()); n++) {
	        const int nindex = remove_index[n];
	        double dm = dm_removed_cell[n];
	        if (dm <= 0.0)
	          continue;

	        const double rho_old = BaryonField[DensNum][nindex];
	        const double mold = rho_old * cell_volume_code;
	        if (mold <= 0.0) {
	          n_floor_rejected_cells++;
	          continue;
	        }

	        double mnew = mold - dm;
	        const double mnew_floor =
	          max(density_floor_mass,
	              (1.0 - BHAccretionMaxCellRemovalFraction) * mold);
	        if (mnew < mnew_floor) {
	          mnew = mnew_floor;
	          dm = max(0.0, mold - mnew);
	          n_cells_clamped++;
	        }
	        if (dm <= 0.0 || mnew <= 0.0) {
	          n_floor_rejected_cells++;
	          continue;
	        }

	        const double rho_new = mnew / cell_volume_code;
	        const double vx = remove_velx[n];
	        const double vy = remove_vely[n];
	        const double vz = remove_velz[n];
	        const double v2 = vx*vx + vy*vy + vz*vz;
	        const int use_dual_energy =
	          (GENum >= 0 && DualEnergyFormalism) ? TRUE : FALSE;
	        const double ge_old =
	          (GENum >= 0) ? BaryonField[GENum][nindex] : -1.0;
	        const double te_old = BaryonField[TENum][nindex];
	        const double eint_old =
	          BHAccretionInternalEnergySpecific(te_old, ge_old, vx, vy, vz,
	                                            use_dual_energy != FALSE);
	        const double pressure_old =
	          BHAccretionPressureFromState(rho_old, te_old, ge_old, vx, vy, vz,
	                                       use_dual_energy != FALSE, Gamma);

	        if (isfinite(rho_old))
	          min_density_before = min(min_density_before, rho_old);
	        if (isfinite(eint_old))
	          min_eint_before = min(min_eint_before, eint_old);
	        if (isfinite(pressure_old))
	          min_pressure_before = min(min_pressure_before, pressure_old);

	        removal_cells++;

	        if (HydroMethod == PPM_DirectEuler) {
	          const double ge_floor =
	            BHAccretionGasEnergyFloor(rho_new, tiny_number, Gamma);
	          double eint_new = eint_old;
	          if (!isfinite(eint_new) || eint_new < ge_floor) {
	            eint_new = ge_floor;
	            n_cells_clamped++;
	          }
	          BaryonField[TENum][nindex] = float(eint_new + 0.5*v2);
	          if (GENum >= 0 && DualEnergyFormalism)
	            BaryonField[GENum][nindex] = float(eint_new);
	        } else if (HydroMethod == Zeus_Hydro) {
	          const double ge_floor =
	            BHAccretionGasEnergyFloor(rho_new, tiny_number, Gamma);
	          double eint_new = BaryonField[TENum][nindex];
	          if (!isfinite(eint_new) || eint_new < ge_floor) {
	            eint_new = ge_floor;
	            n_cells_clamped++;
	          }
	          BaryonField[TENum][nindex] = float(eint_new);
	        } else {
	          ENZO_FAIL("BHAccretion removal does not support RK Hydro or RK MHD.");
	        }

	        BaryonField[DensNum][nindex] = float(rho_new);

	        const double ge_after =
	          (GENum >= 0) ? BaryonField[GENum][nindex] : -1.0;
	        const double te_after = BaryonField[TENum][nindex];
	        const double eint_after =
	          BHAccretionInternalEnergySpecific(te_after, ge_after, vx, vy, vz,
	                                            use_dual_energy != FALSE);
	        const double pressure_after =
	          BHAccretionPressureFromState(rho_new, te_after, ge_after, vx, vy, vz,
	                                       use_dual_energy != FALSE, Gamma);

	        if (isfinite(rho_new))
	          min_density_after = min(min_density_after, rho_new);
	        if (isfinite(eint_after))
	          min_eint_after = min(min_eint_after, eint_after);
	        if (isfinite(pressure_after))
	          min_pressure_after = min(min_pressure_after, pressure_after);

	        const double consistency =
	          te_after - (eint_after + 0.5*v2);
	        const double consistency_tol =
	          1.0e-6 * max(1.0, fabs(te_after));
	        int cell_energy_update_consistent = TRUE;
	        if (HydroMethod == PPM_DirectEuler &&
	            (!isfinite(consistency) || fabs(consistency) > consistency_tol)) {
	          n_energy_consistency_fail++;
	          cell_energy_update_consistent = FALSE;
	        }

	        const double removal_fraction = dm / mold;
	        if (isfinite(removal_fraction))
	          max_cell_removal_fraction =
	            max(max_cell_removal_fraction, removal_fraction);

	        if (BHAccretionVerbose >= 2) {
	          fprintf(logptr,
	                  "[BHACCR_CELL] rank=%d step=%d level=%d grid_id=%d bh_id=%lld "
	                  "i=%d j=%d k=%d index=%d dm=%.15e removal_fraction=%.8g "
	                  "rho_before=%.15e rho_after=%.15e "
	                  "total_energy_before=%.15e total_energy_after=%.15e "
	                  "gas_energy_before=%.15e gas_energy_after=%.15e "
	                  "vx_before=%.15e vy_before=%.15e vz_before=%.15e "
	                  "vx_after=%.15e vy_after=%.15e vz_after=%.15e "
	                  "px_before=%.15e py_before=%.15e pz_before=%.15e "
	                  "px_after=%.15e py_after=%.15e pz_after=%.15e "
	                  "pressure_before=%.15e pressure_after=%.15e "
	                  "energy_update_consistent=%d\n",
	                  MyProcessorNumber, cycle_number, level, this->GetGridID(),
	                  (long long) ParticleNumber[p],
	                  nindex % nx, (nindex / nx) % ny, nindex / (nx*ny), nindex,
	                  dm, removal_fraction,
	                  rho_old, rho_new, te_old, te_after, ge_old, ge_after,
	                  vx, vy, vz, vx, vy, vz,
	                  rho_old*vx, rho_old*vy, rho_old*vz,
	                  rho_new*vx, rho_new*vy, rho_new*vz,
	                  pressure_old, pressure_after,
	                  cell_energy_update_consistent);
	        }

	        removed_mass += dm;
        delta_px += dm * remove_velx[n];
        delta_py += dm * remove_vely[n];
        delta_pz += dm * remove_velz[n];

        if (BaryonField[NumberOfBaryonFields] != NULL) {
          BaryonField[NumberOfBaryonFields][nindex] = 1.0f;
          n_sf_blocked_cells++;
        }
      }
    }

    dm_removed_total = removed_mass;

    if (dm_requested > 0.0) {
      frac_gas = dm_removed_total / dm_requested;
      frac_gas = min(1.0, max(0.0, frac_gas));
      const double post_flag_tol = 1.0e-8 * max(1.0, dm_requested);
      if (dm_removed_total + post_flag_tol < dm_requested)
        removal_gas_limited = 1;
	    } else {
	      frac_gas = 1.0;
	      removal_gas_limited = 0;
	    }

	    const double log_min_density_before =
	      (min_density_before == DBL_MAX) ? -1.0 : min_density_before;
	    const double log_min_density_after =
	      (min_density_after == DBL_MAX) ? -1.0 : min_density_after;
	    const double log_min_eint_before =
	      (min_eint_before == DBL_MAX) ? -1.0 : min_eint_before;
	    const double log_min_eint_after =
	      (min_eint_after == DBL_MAX) ? -1.0 : min_eint_after;
	    const double log_min_pressure_before =
	      (min_pressure_before == DBL_MAX) ? -1.0 : min_pressure_before;
	    const double log_min_pressure_after =
	      (min_pressure_after == DBL_MAX) ? -1.0 : min_pressure_after;

	    if (BHAccretionVerbose >= 1) {
	      fprintf(logptr,
	              "[BHACCR_SRC] rank=%d step=%d level=%d grid_id=%d bh_id=%lld "
	              "kernel_radius_over_dx=%.8g kernel_cells_requested=%lld "
	              "kernel_cells_valid=%d kernel_cells_active=%d kernel_truncated=%d "
	              "accretion_source_rejected=%d removal_cells=%d "
	              "removal_kernel_gas=%.15e removal_kernel_safe_gas=%.15e "
	              "dm_requested=%.15e dm_removed=%.15e "
	              "max_cell_removal_fraction=%.8g max_allowed_cell_removal_fraction=%.8g "
	              "min_density_before=%.15e min_density_after=%.15e "
	              "min_internal_energy_before=%.15e min_internal_energy_after=%.15e "
	              "min_pressure_before=%.15e min_pressure_after=%.15e "
	              "cells_clamped=%d cells_rejected_by_floor=%d "
	              "energy_update_consistent=%d\n",
	              MyProcessorNumber, cycle_number, level, this->GetGridID(),
	              (long long) ParticleNumber[p],
	              kernel_cells_per_radius, n_kernel_cells_requested,
	              n_kernel_cells, n_kernel_active_cells, kernel_truncated,
	              accretion_source_rejected, removal_cells,
	              removal_kernel_gas, removal_kernel_safe_gas,
	              dm_requested, dm_removed_total,
	              max_cell_removal_fraction,
	              BHAccretionMaxCellRemovalFraction,
	              log_min_density_before, log_min_density_after,
	              log_min_eint_before, log_min_eint_after,
	              log_min_pressure_before, log_min_pressure_after,
	              n_cells_clamped, n_floor_rejected_cells,
	              (n_energy_consistency_fail == 0) ? 1 : 0);
	    }

    const double mdot_realized =
      (dt_code > 0.0) ? (dm_removed_total / dt_code) : 0.0;
    const double realized_over_requested =
      (mdot_actual > 0.0) ? (mdot_realized / mdot_actual) : 0.0;

	    double mdot_hot_realized = mdot_hot_capped * frac_gas;
    double mdot_cold_realized = mdot_cold_capped * frac_gas;

    if (dt_code > 0.0) {
      const double realized_sum = mdot_hot_realized + mdot_cold_realized;
      const double rel_tol = 1.0e-10 * max(1.0, max(fabs(realized_sum), fabs(mdot_realized)));
      if (fabs(realized_sum - mdot_realized) > rel_tol)
        ENZO_FAIL("BHAccretion invariant failed: realized channel rates do not match dm_removed/dt.");
    }

    const double bh_mass_new_code = bh_mass_code + dm_removed_total;
    ParticleMass[p] = float(bh_mass_new_code);
    const double bh_mass_new_msun = bh_mass_new_code * mass_to_msun;

    const double dm_removed_msun = dm_removed_total * mass_to_msun;

    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_RESERVOIR_MASS) {
      float &v = ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_RESERVOIR_MASS][p];
      if (!isfinite(v) || v < 0.0f)
        v = 0.0f;
    }

    int bh_formation_mass_initialized = 0;
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BH_FORMATION_MASS) {
      float &v = ParticleAttribute[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS][p];
      if (!isfinite(v) || v <= 0.0f) {
        v = float(bh_mass_code);
        bh_formation_mass_initialized = 1;
      }
    }

    int bh_accreted_mass_initialized = 0;
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS) {
      float &v = ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS][p];
      if (!isfinite(v) || v < 0.0f) {
        v = 0.0f;
        bh_accreted_mass_initialized = 1;
      }
      v += float(dm_removed_total);
    }

    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_LAST_REDSHIFT) {
      float &v = ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_REDSHIFT][p];
      if (!isfinite(v))
        v = -1.0f;
      if (dm_removed_total > 0.0)
        v = float(zred);
    }

    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_LAST_EDD_RATIO)
      ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_EDD_RATIO][p] = float(f_edd);
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_ACTUAL)
      ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_ACTUAL][p] = float(mdot_actual);
    const int realized_attr_written =
      (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED &&
       ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED] != NULL) ? 1 : 0;
    if (realized_attr_written)
      ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED][p] =
        float(mdot_realized);

    const double mass_update_rhs = bh_mass_code + dm_removed_total;
    const double mass_update_scale = max(fabs(bh_mass_new_code), fabs(mass_update_rhs));
    const double mass_invariant_tol =
      (mass_update_scale > 0.0) ? 1.0e-10 * mass_update_scale : 1.0e-30;
    if (fabs(bh_mass_new_code - mass_update_rhs) > mass_invariant_tol)
      ENZO_FAIL("BHAccretion invariant failed: BH mass update mismatch.");

    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS &&
        NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BH_FORMATION_MASS) {
      const double formation_mass_code =
        ParticleAttribute[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS][p];
      const double accreted_mass_code =
        ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS][p];
      const double invariant_rhs = formation_mass_code + accreted_mass_code;
      const double mass_scale = max(fabs(bh_mass_new_code), fabs(invariant_rhs));
      const double invariant_tol = (mass_scale > 0.0) ? 1.0e-8 * mass_scale : 1.0e-30;
      const double invariant_residual = bh_mass_new_code - invariant_rhs;
      if (fabs(invariant_residual) > invariant_tol) {
        const double rel_residual =
          (mass_scale > 0.0) ? fabs(invariant_residual) / mass_scale : 0.0;
        fprintf(stderr,
                "[BHACCR_INVARIANT_FAIL] rank=%d processor=%d step=%d level=%d "
                "grid_id=%d particle_index=%d bh_id=%lld particle_type=%d "
                "NumberOfParticleAttributes=%d ParticleMass=%.15e "
                "bh_mass_code=%.15e bh_mass_new_code=%.15e "
                "BHFormationMass=%.15e BHAccretedMass=%.15e "
                "expected_mass=%.15e residual=%.15e abs_residual=%.15e "
                "rel_residual=%.15e tolerance=%.15e Time=%.15e dt=%.15e "
                "dm_requested=%.15e dm_removed=%.15e mdot_actual=%.15e "
                "mdot_realized=%.15e formation_mass_initialized=%d "
                "accreted_mass_initialized=%d pos_x=%.15e pos_y=%.15e pos_z=%.15e\n",
                MyProcessorNumber, ProcessorNumber, cycle_number, level,
                this->GetGridID(), p, (long long) ParticleNumber[p],
                ParticleType[p], NumberOfParticleAttributes,
                double(ParticleMass[p]), bh_mass_code, bh_mass_new_code,
                formation_mass_code, accreted_mass_code, invariant_rhs,
                invariant_residual, fabs(invariant_residual), rel_residual,
                invariant_tol, double(Time), dt_code, dm_requested,
                dm_removed_total, mdot_actual, mdot_realized,
                bh_formation_mass_initialized, bh_accreted_mass_initialized,
                double(ParticlePosition[0][p]),
                (GridRank > 1) ? double(ParticlePosition[1][p]) : 0.0,
                (GridRank > 2) ? double(ParticlePosition[2][p]) : 0.0);
        fflush(stderr);
        ENZO_FAIL("BHAccretion invariant failed: BH mass != BHFormationMass + BHAccretedMass.");
      }
    }

    if (BHAccretionVerbose >= 1) {
	      fprintf(logptr,
	              "[BHACCR_DEBUG] step=%d level=%d bh_id=%lld "
	              "Mdot_total_capped=%.15e dm_requested=%.15e "
	              "removal_kernel_gas=%.15e removal_kernel_safe_gas=%.15e "
	              "dm_removed=%.15e\n",
	              cycle_number, level, (long long) ParticleNumber[p],
	              mdot_total_capped, dm_requested, removal_kernel_gas,
	              removal_kernel_safe_gas, dm_removed_total);

      if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS &&
          NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BH_FORMATION_MASS) {
        const double formation_mass_code =
          ParticleAttribute[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS][p];
        const double accreted_mass_code =
          ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS][p];
        const double mass_sum = formation_mass_code + accreted_mass_code;
        const double mass_diff = bh_mass_new_code - mass_sum;
        fprintf(logptr,
                "[BHACCR_DEBUG] step=%d level=%d bh_id=%lld "
                "bh_mass_code=%.15e formation_mass_code=%.15e "
                "accreted_mass_code=%.15e sum=%.15e diff=%.15e\n",
                cycle_number, level, (long long) ParticleNumber[p],
                bh_mass_new_code, formation_mass_code, accreted_mass_code,
                mass_sum, mass_diff);
      }
    }

    const double delta_p_mag = sqrt(delta_px*delta_px + delta_py*delta_py + delta_pz*delta_pz);
    const double acc_ignored_dv =
      (bh_mass_new_code > 0.0) ? (delta_p_mag / bh_mass_new_code) : 0.0;
    const double acc_ignored_dv_kms = acc_ignored_dv * vel_to_cgs / 1.0e5;

    const double bh_speed = sqrt(bh_vel_x*bh_vel_x + bh_vel_y*bh_vel_y + bh_vel_z*bh_vel_z);
    const double p_bh_mag = bh_mass_code * bh_speed;
    const double acc_ignored_p_frac = delta_p_mag / (p_bh_mag + 1.0e-30);

    const int acc_momentum_warn =
      (acc_ignored_dv_kms > double(BHAccretionIgnoredDVWarn) ||
       acc_ignored_p_frac > double(BHAccretionIgnoredPFracWarn)) ? 1 : 0;

	    const double accretion_wall_ms = 1000.0 * (ReturnWallTime() - t0_all);

	    if (BHAccretionVerbose >= 1) {
	      if (kernel_cells_per_radius < 1.5)
	        fprintf(logptr,
                "[BHACCR_WARN] step=%d level=%d bh_id=%lld kernel_radius_over_dx=%.6g "
                "kernel under-resolved (<1.5 cells).\n",
                cycle_number, level, (long long) ParticleNumber[p],
                kernel_cells_per_radius);
      if (kernel_cells_per_radius > 3.0)
        fprintf(logptr,
                "[BHACCR_WARN] step=%d level=%d bh_id=%lld kernel_radius_over_dx=%.6g "
                "kernel may exceed ghost-zone support (>3 cells).\n",
	                cycle_number, level, (long long) ParticleNumber[p],
	                kernel_cells_per_radius);
	      if (accretion_source_rejected)
	        fprintf(logptr,
	                "[BHACCR_WARN] step=%d level=%d bh_id=%lld "
	                "accretion_source_rejected=1 kernel_radius_over_dx=%.6g "
	                "NumberOfGhostZones=%d kernel_truncated_by_grid=%d; "
	                "diagnostics retained but gas removal skipped.\n",
	                cycle_number, level, (long long) ParticleNumber[p],
	                kernel_cells_per_radius, NumberOfGhostZones,
	                kernel_truncated_by_grid);
	      if (n_kernel_cells > 0 && (2*n_fallback_cells > n_kernel_cells))
        fprintf(logptr,
                "[BHACCR_WARN] step=%d level=%d bh_id=%lld n_fallback_cells=%d "
                "n_kernel_cells=%d fallback fraction exceeds 50%%.\n",
                cycle_number, level, (long long) ParticleNumber[p],
                n_fallback_cells, n_kernel_cells);
      if (hot.mass > 0.0 && alpha_boost >= BHAccretionAlphaMax)
        fprintf(logptr,
                "[BHACCR_WARN] step=%d level=%d bh_id=%lld alpha_boost=%.6g "
                "at configured cap alpha_max=%.6g.\n",
                cycle_number, level, (long long) ParticleNumber[p],
                alpha_boost, double(BHAccretionAlphaMax));

      if (removal_gas_limited)
        fprintf(logptr,
                "[BHACCR_WARN] step=%d level=%d bh_id=%lld removal_gas_limited=1 "
                "dm_requested=%.8e dm_removed=%.8e frac_gas=%.8g\n",
                cycle_number, level, (long long) ParticleNumber[p],
                dm_requested, dm_removed_total, frac_gas);

      if (acc_momentum_warn)
        fprintf(logptr,
                "[BHACCR_WARN] step=%d level=%d bh_id=%lld acc_ignored_dv_kms=%.8g "
                "acc_ignored_p_frac=%.8g threshold_dv_kms=%.8g threshold_p_frac=%.8g\n",
                cycle_number, level, (long long) ParticleNumber[p],
                acc_ignored_dv_kms, acc_ignored_p_frac,
                double(BHAccretionIgnoredDVWarn), double(BHAccretionIgnoredPFracWarn));
    }

    if (BHAccretionVerbose >= 1) {
      fprintf(logptr,
              "[BHACCR] step=%d level=%d z=%.8g bh_id=%lld bh_mass=%.8g "
              "f_hot=%.8g f_cold=%.8g n_hot_cells=%d n_cold_cells=%d n_fallback_cells=%d "
              "rho_hot_avg=%.8e rho_cold_avg=%.8e T_hot_avg=%.8e T_cold_avg=%.8e "
              "cs_hot_avg=%.8e cs_cold_avg=%.8e V_rot_cold=%.8e "
              "v_rel_hot=%.8e v_rel_cold=%.8e "
              "Mdot_hot_raw=%.8e Mdot_cold_raw=%.8e Mdot_total_raw=%.8e "
              "Mdot_Edd=%.8e f_Edd=%.8e alpha_boost=%.8g f_AM=%.8g "
              "Mdot_actual=%.8e Mdot_hot_actual=%.8e Mdot_cold_actual=%.8e "
              "Mdot_realized=%.8e mdot_actual_code=%.15e mdot_realized_code=%.15e "
	              "cap_active=%d accretion_diag_wall_ms=%.4f "
	              "dm_requested=%.8e dm_removed=%.8e dm_removed_msun=%.8e "
	              "frac_cap=%.8g frac_gas=%.8g realized_over_requested=%.8e "
              "realized_attr_written=%d removal_gas_limited=%d "
	              "kernel_truncated=%d accretion_source_rejected=%d "
	              "removal_kernel_safe_gas=%.8e "
	              "Mdot_hot_realized=%.8e Mdot_cold_realized=%.8e bh_mass_new=%.8g "
              "removal_cells=%d n_sf_blocked_cells=%d "
              "acc_ignored_dv_kms=%.8g acc_ignored_p_frac=%.8g acc_momentum_warn=%d "
              "accretion_wall_ms=%.4f lambda_edd=%.8e edd_factor=%.8g\n",
              cycle_number, level, zred, (long long) ParticleNumber[p], bh_mass_msun,
              f_hot, f_cold, n_hot_cells, n_cold_cells, n_fallback_cells,
              hot_rho_avg * rho_to_cgs, cold_rho_avg * rho_to_cgs,
              hot_t_avg, cold_t_avg,
              hot_cs_avg * vel_to_cgs, cold_cs_avg * vel_to_cgs,
              v_rot_cold * vel_to_cgs,
              v_rel_hot * vel_to_cgs, v_rel_cold * vel_to_cgs,
              mdot_hot_raw * mdot_to_msunyr,
              mdot_cold_raw * mdot_to_msunyr,
              mdot_total_raw * mdot_to_msunyr,
              mdot_edd_code * mdot_to_msunyr, f_edd,
              alpha_boost, f_am,
              mdot_actual * mdot_to_msunyr,
              mdot_hot_actual * mdot_to_msunyr,
              mdot_cold_actual * mdot_to_msunyr,
              mdot_realized * mdot_to_msunyr,
              mdot_actual, mdot_realized,
              cap_active, accretion_diag_wall_ms,
	              dm_requested, dm_removed_total, dm_removed_msun,
	              frac_cap, frac_gas, realized_over_requested,
              realized_attr_written, removal_gas_limited,
	              kernel_truncated, accretion_source_rejected,
	              removal_kernel_safe_gas,
	              mdot_hot_realized * mdot_to_msunyr,
              mdot_cold_realized * mdot_to_msunyr,
              bh_mass_new_msun,
              removal_cells, n_sf_blocked_cells,
              acc_ignored_dv_kms, acc_ignored_p_frac, acc_momentum_warn,
              accretion_wall_ms, lambda_edd, double(BHFeedbackEddingtonFactor));
    }
  }

  return SUCCESS;
}
