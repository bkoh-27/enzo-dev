/***********************************************************************
/
/  GRID CLASS (BH ACCRETION DIAGNOSTICS, PHASE A)
/
/  PURPOSE:
/    Evaluate per-BH local gas diagnostics and log two-channel accretion
/    rates without removing gas or changing BH masses.
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include <float.h>
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

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL)
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  std::vector<float> temperature(size, 0.0f);
  if (this->ComputeTemperatureField(&temperature[0]) == FAIL)
    ENZO_FAIL("Error in grid->ComputeTemperatureField for BH accretion diagnostics.");

  const int cooling_available = (RadiativeCooling != 0);
  std::vector<float> cooling_time(size, 0.0f);
  if (cooling_available)
    if (this->ComputeCoolingTime(&cooling_time[0]) == FAIL)
      ENZO_FAIL("Error in grid->ComputeCoolingTime for BH accretion diagnostics.");

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

  float a_phys = 1.0f, h_param = 1.0f, box_kpch = 0.0f;
  const double kernel_radius_code =
    BHAccretionKernelRadiusCode(BHAccretionKernelRadius, Time, LengthUnits,
                                &a_phys, &h_param, &box_kpch);
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

  FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;

  for (int p = 0; p < NumberOfParticles; p++) {
    if (!BHAccretionIsBHType(ParticleType[p]))
      continue;

    if (!BHAccretionRunEveryTimestep && SubgridPointer != NULL) {
      FLOAT bh_pos[MAX_DIMENSION];
      bh_pos[0] = ParticlePosition[0][p];
      bh_pos[1] = (GridRank > 1) ? ParticlePosition[1][p] : 0.0;
      bh_pos[2] = (GridRank > 2) ? ParticlePosition[2][p] : 0.0;

      int covered_by_child = FALSE;
      for (HierarchyEntry *sub = SubgridPointer; sub != NULL;
           sub = sub->NextGridThisLevel) {
        if (sub->GridData != NULL && sub->GridData->PointInGrid(bh_pos)) {
          covered_by_child = TRUE;
          break;
        }
      }
      if (covered_by_child)
        continue;
    }

    const double t0 = ReturnWallTime();

    const double bh_mass_code = ParticleMass[p];
    if (!isfinite(bh_mass_code) || bh_mass_code <= 0.0)
      continue;

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

    BHAccretionChannelAccumulator hot = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                         0.0, 0.0, 0.0, 0};
    BHAccretionChannelAccumulator cold = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                          0.0, 0.0, 0.0, 0};

    int n_kernel_cells = 0;
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
        f_am = min(1.0, double(BHAccretionCVisc) * pow(cold_cs_avg / v_rot_cold, 3.0));
      else
        f_am = 1.0;
      if (f_am < 0.0)
        f_am = 0.0;
      mdot_cold_raw = f_am * mdot_bondi_cold;
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

    const double f_edd = (mdot_edd_code > 0.0) ? mdot_total_raw / mdot_edd_code : 0.0;
    double mdot_actual = 0.0;
    double mdot_hot_actual = 0.0;
    double mdot_cold_actual = 0.0;
    int cap_active = 0;
    if (mdot_edd_code > 0.0 && mdot_total_raw > mdot_edd_code) {
      const double frac = mdot_edd_code / mdot_total_raw;
      mdot_actual = mdot_edd_code;
      mdot_hot_actual = mdot_hot_raw * frac;
      mdot_cold_actual = mdot_cold_raw * frac;
      cap_active = 1;
    } else {
      mdot_actual = mdot_total_raw;
      mdot_hot_actual = mdot_hot_raw;
      mdot_cold_actual = mdot_cold_raw;
    }

    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS) {
      float &v = ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS][p];
      if (!isfinite(v) || v < 0.0f)
        v = 0.0f;
    }
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_RESERVOIR_MASS) {
      float &v = ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_RESERVOIR_MASS][p];
      if (!isfinite(v) || v < 0.0f)
        v = 0.0f;
    }
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_LAST_REDSHIFT) {
      float &v = ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_REDSHIFT][p];
      if (!isfinite(v))
        v = -1.0f;
    }
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_LAST_EDD_RATIO)
      ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_EDD_RATIO][p] = float(f_edd);
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BH_FORMATION_MASS) {
      float &v = ParticleAttribute[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS][p];
      if (!isfinite(v) || v <= 0.0f)
        v = float(bh_mass_msun);
    }

    const double accretion_diag_wall_ms = 1000.0 * (ReturnWallTime() - t0);

    if (BHAccretionVerbose >= 1) {
      const double kernel_cells_per_radius = kernel_radius_code / cell_width;
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
    }

    if (BHAccretionVerbose >= 1) {
      fprintf(logptr,
              "[BHACCR] step=%d level=%d z=%.8g bh_id=%lld bh_mass=%.8g "
              "f_hot=%.8g f_cold=%.8g n_hot_cells=%d n_cold_cells=%d n_fallback_cells=%d "
              "rho_hot_avg=%.8e rho_cold_avg=%.8e T_hot_avg=%.8e T_cold_avg=%.8e "
              "cs_hot_avg=%.8e V_rot_cold=%.8e v_rel_hot=%.8e v_rel_cold=%.8e "
              "Mdot_hot_raw=%.8e Mdot_cold_raw=%.8e Mdot_total_raw=%.8e "
              "Mdot_Edd=%.8e f_Edd=%.8e alpha_boost=%.8g f_AM=%.8g "
              "Mdot_actual=%.8e Mdot_hot_actual=%.8e Mdot_cold_actual=%.8e "
              "cap_active=%d accretion_diag_wall_ms=%.4f\n",
              cycle_number, level, zred, (long long) ParticleNumber[p], bh_mass_msun,
              f_hot, f_cold, n_hot_cells, n_cold_cells, n_fallback_cells,
              hot_rho_avg * rho_to_cgs, cold_rho_avg * rho_to_cgs,
              hot_t_avg, cold_t_avg,
              hot_cs_avg * vel_to_cgs, v_rot_cold * vel_to_cgs,
              v_rel_hot * vel_to_cgs, v_rel_cold * vel_to_cgs,
              mdot_hot_raw * mdot_to_msunyr,
              mdot_cold_raw * mdot_to_msunyr,
              mdot_total_raw * mdot_to_msunyr,
              mdot_edd_code * mdot_to_msunyr, f_edd,
              alpha_boost, f_am,
              mdot_actual * mdot_to_msunyr,
              mdot_hot_actual * mdot_to_msunyr,
              mdot_cold_actual * mdot_to_msunyr,
              cap_active, accretion_diag_wall_ms);
    }
  }

  return SUCCESS;
}
