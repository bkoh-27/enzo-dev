/***********************************************************************
/
/  GRID CLASS (BH FEEDBACK DIAGNOSTICS, PHASE A)
/
/  PURPOSE:
/    Compute diagnostics for thermal and kinetic BH feedback budgets and
/    kernel statistics without modifying grid fields or feedback attributes.
/
/  PHASE A GUARANTEES:
/    - No energy/momentum deposition.
/    - No baryon-field writes.
/    - No mutation of BH feedback reservoir attributes.
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <map>

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

static int BHFeedbackIsBHType(int type)
{
  return (type == PARTICLE_TYPE_MBH || type == PARTICLE_TYPE_BLACK_HOLE);
}

static int BHFeedbackIsNewlySeededThisPass(PINT particle_id)
{
  return (particle_id == INT_UNDEFINED);
}

struct BHFeedbackParticleOrder {
  PINT *ParticleIDs;
  BHFeedbackParticleOrder(PINT *ids) : ParticleIDs(ids) { }
  bool operator()(int a, int b) const {
    if (ParticleIDs[a] < ParticleIDs[b])
      return true;
    if (ParticleIDs[a] > ParticleIDs[b])
      return false;
    return (a < b);
  }
};

static double BHFeedbackKernelRadiusCode(float KernelRadiusPhysKpc,
                                         FLOAT time,
                                         float LengthUnits)
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
    box_kpch = LengthUnits / kpc_cm;
    if (box_kpch > 0.0f)
      kernel_radius_code = KernelRadiusPhysKpc / box_kpch;
  }

  return kernel_radius_code;
}

int grid::BHFeedbackHandler(HierarchyEntry* SubgridPointer,
                            int level, int cycle_number,
                            float dtLevelAbove)
{
  (void) dtLevelAbove;

  if (!BHFeedbackMethod)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0 || NumberOfParticles <= 0)
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL)
    ENZO_FAIL("Error in IdentifyPhysicalQuantities for BH feedback.");

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  std::vector<float> temperature(size, 0.0f);
  if (this->ComputeTemperatureField(&temperature[0]) == FAIL)
    ENZO_FAIL("Error in grid->ComputeTemperatureField for BH feedback.");

  float DensityUnits = 1.0f, LengthUnits = 1.0f, TemperatureUnits = 1.0f;
  float TimeUnits = 1.0f, VelocityUnits = 1.0f;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL)
    ENZO_FAIL("Error in GetUnits for BH feedback.");

  const float cell_width = float(CellWidth[0][0]);
  if (cell_width <= 0.0f)
    return SUCCESS;

  const double kernel_radius_code =
    BHFeedbackKernelRadiusCode(BHFeedbackKernelRadius, Time, LengthUnits);
  if (kernel_radius_code <= 0.0)
    return SUCCESS;

  FLOAT a = 1.0, dadt = 0.0;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  const double zred = ComovingCoordinates ?
    double((1.0 + InitialRedshift) / a - 1.0) : 0.0;

  const double dt_code = max(0.0, double(this->dtFixed));
  const double dt_cgs = dt_code * TimeUnits;
  const double mass_units = DensityUnits * pow(LengthUnits, 3.0);
  const double mass_rate_to_cgs = (TimeUnits > 0.0f) ? mass_units / TimeUnits : 0.0;
  const double mass_to_msun = mass_units / SolarMass;
  const double cell_volume_code = pow(cell_width, 3.0);

  const int nx = GridDimension[0];
  const int ny = GridDimension[1];
  const int nz = GridDimension[2];
  const int xo = 1;
  const int yo = nx;
  const int zo = nx*ny;

  const int isx = this->GetGridStartIndex(0);
  const int isy = this->GetGridStartIndex(1);
  const int isz = this->GetGridStartIndex(2);
  const int iex = this->GetGridEndIndex(0);
  const int iey = this->GetGridEndIndex(1);
  const int iez = this->GetGridEndIndex(2);

  std::vector<int> bh_particles;
  bh_particles.reserve(NumberOfParticles);
  for (int p = 0; p < NumberOfParticles; p++)
    if (BHFeedbackIsBHType(ParticleType[p]))
      bh_particles.push_back(p);

  std::sort(bh_particles.begin(), bh_particles.end(),
            BHFeedbackParticleOrder(ParticleNumber));

  FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;
  static int warned_method_two = FALSE;
  static std::map<PINT, double> feedback_reservoir_diag_cache;

  for (size_t ip = 0; ip < bh_particles.size(); ip++) {
    const double t0_all = ReturnWallTime();
    const int p = bh_particles[ip];

    FLOAT bh_pos[MAX_DIMENSION] = {ParticlePosition[0][p],
                                   ParticlePosition[1][p],
                                   ParticlePosition[2][p]};
    if (!this->PointInGrid(bh_pos))
      continue;

    if (!BHAccretionRunEveryTimestep && SubgridPointer != NULL)
      continue;

    const double bh_mass_code = ParticleMass[p];
    const double bh_mass_msun = bh_mass_code * mass_to_msun;

    if (BHFeedbackIsNewlySeededThisPass(ParticleNumber[p])) {
      if (BHFeedbackVerbose >= 1) {
        const double feedback_wall_ms = 1000.0 * (ReturnWallTime() - t0_all);
        fprintf(logptr,
                "[BHFDBK] step=%d level=%d z=%.8g bh_id=%lld bh_mass=%.8g "
                "feedback_mode=OFF f_Edd=-1 L_feedback=-1 E_requested=-1 "
                "reservoir_before=-1 reservoir_after=-1 burst_diag=0 "
                "E_deposited=0 p_requested=-1 p_deposited=0 "
                "feedback_kernel_cells=0 feedback_kernel_active_cells=0 "
                "T_before_mean=-1 T_after_mean=-1 n_sf_blocked_feedback=0 "
                "newly_seeded_skip=1 feedback_wall_ms=%.4f\n",
                cycle_number, level, zred, (long long) ParticleNumber[p], bh_mass_msun,
                feedback_wall_ms);
      }
      continue;
    }

    double f_edd = 0.0;
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_LAST_EDD_RATIO) {
      f_edd = ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_EDD_RATIO][p];
      if (!isfinite(f_edd))
        f_edd = 0.0;
    }

    double mdot_actual_code = 0.0;
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_ACTUAL) {
      mdot_actual_code = ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_ACTUAL][p];
      if (!isfinite(mdot_actual_code) || mdot_actual_code < 0.0)
        mdot_actual_code = 0.0;
    }

    const char *feedback_mode = (f_edd > BHFeedbackModeThreshold) ? "THERMAL" : "KINETIC";
    if (BHFeedbackMethod == 2 && BHFeedbackVerbose >= 1 && !warned_method_two) {
      fprintf(logptr,
              "[BHFDBK_WARN] step=%d level=%d BHFeedbackMethod=2 requested; "
              "active two-mode deposition requires Phase C. Running diagnostics-only framework.\n",
              cycle_number, level);
      warned_method_two = TRUE;
    }

    const double mdot_cgs = mdot_actual_code * mass_rate_to_cgs;
    const double L_feedback = BHAccretionRadiativeEfficiency * mdot_cgs * clight * clight;
    const double E_requested = BHFeedbackThermalEfficiency * L_feedback * dt_cgs;
    const double p_requested =
      BHFeedbackKineticEfficiency * mdot_cgs *
      (double(BHFeedbackWindVelocity) * 1.0e5) * dt_cgs;

    double reservoir_before = 0.0;
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHFDBK_ENERGY_RESERVOIR) {
      reservoir_before = ParticleAttribute[PARTICLE_ATTRIBUTE_BHFDBK_ENERGY_RESERVOIR][p];
      if (!isfinite(reservoir_before) || reservoir_before < 0.0)
        reservoir_before = 0.0;
    }
    const PINT bh_id = ParticleNumber[p];
    const double reservoir_diag_before = feedback_reservoir_diag_cache[bh_id];
    const double reservoir_after = reservoir_diag_before + E_requested;
    feedback_reservoir_diag_cache[bh_id] = reservoir_after;
    const int burst_diag = (reservoir_after >= BHFeedbackMinEnergyBurst) ? 1 : 0;

    const int i0 = int((bh_pos[0] - GridLeftEdge[0]) / cell_width);
    const int j0 = int((bh_pos[1] - GridLeftEdge[1]) / cell_width);
    const int k0 = int((bh_pos[2] - GridLeftEdge[2]) / cell_width);
    const int rcell = max(0, int(ceil(kernel_radius_code / cell_width)));
    const double r2max = kernel_radius_code * kernel_radius_code;

    int n_kernel_cells = 0;
    int n_kernel_active_cells = 0;
    double kernel_gas_mass = 0.0;
    double temp_mass_sum = 0.0;
    double temp_weight_sum = 0.0;

    for (int k = max(0, k0-rcell); k <= min(nz-1, k0+rcell); k++)
      for (int j = max(0, j0-rcell); j <= min(ny-1, j0+rcell); j++)
        for (int i = max(0, i0-rcell); i <= min(nx-1, i0+rcell); i++) {
          const double dx = (double(i)+0.5)*cell_width + GridLeftEdge[0] - bh_pos[0];
          const double dy = (double(j)+0.5)*cell_width + GridLeftEdge[1] - bh_pos[1];
          const double dz = (double(k)+0.5)*cell_width + GridLeftEdge[2] - bh_pos[2];
          const double r2 = dx*dx + dy*dy + dz*dz;
          if (r2 > r2max)
            continue;

          n_kernel_cells++;
          if (i >= isx && i <= iex && j >= isy && j <= iey && k >= isz && k <= iez)
            n_kernel_active_cells++;

          const int n = k*zo + j*yo + i*xo;
          const double rho = BaryonField[DensNum][n];
          if (rho <= 0.0)
            continue;
          const double mcell = rho * cell_volume_code;
          if (mcell <= 0.0)
            continue;
          kernel_gas_mass += mcell;
          temp_mass_sum += mcell * temperature[n];
          temp_weight_sum += mcell;
        }

    const double T_before_mean = (temp_weight_sum > 0.0) ? (temp_mass_sum / temp_weight_sum) : 0.0;
    const double T_after_mean = T_before_mean;
    const double E_deposited = 0.0;
    const double p_deposited = 0.0;
    const int n_sf_blocked_feedback = 0;

    const double kernel_cells_per_radius = kernel_radius_code / cell_width;
    if (BHFeedbackVerbose >= 1) {
      if (kernel_cells_per_radius < 1.5)
        fprintf(logptr,
                "[BHFDBK_WARN] step=%d level=%d bh_id=%lld kernel_radius_over_dx=%.6g "
                "kernel under-resolved (<1.5 cells).\n",
                cycle_number, level, (long long) ParticleNumber[p], kernel_cells_per_radius);
      if (kernel_cells_per_radius > 3.0)
        fprintf(logptr,
                "[BHFDBK_WARN] step=%d level=%d bh_id=%lld kernel_radius_over_dx=%.6g "
                "kernel may exceed ghost-zone support (>3 cells).\n",
                cycle_number, level, (long long) ParticleNumber[p], kernel_cells_per_radius);

      const double feedback_wall_ms = 1000.0 * (ReturnWallTime() - t0_all);
      fprintf(logptr,
              "[BHFDBK] step=%d level=%d z=%.8g bh_id=%lld bh_mass=%.8g "
              "feedback_mode=%s f_Edd=%.8e L_feedback=%.8e E_requested=%.8e "
              "reservoir_before=%.8e reservoir_after=%.8e burst_diag=%d "
              "E_deposited=%.8e p_requested=%.8e p_deposited=%.8e "
              "feedback_kernel_cells=%d feedback_kernel_active_cells=%d "
              "feedback_kernel_gas_msun=%.8e "
              "T_before_mean=%.8e T_after_mean=%.8e n_sf_blocked_feedback=%d "
              "newly_seeded_skip=0 feedback_wall_ms=%.4f\n",
              cycle_number, level, zred, (long long) ParticleNumber[p], bh_mass_msun,
              feedback_mode, f_edd, L_feedback, E_requested,
              reservoir_before, reservoir_after, burst_diag,
              E_deposited, p_requested, p_deposited,
              n_kernel_cells, n_kernel_active_cells,
              kernel_gas_mass * mass_to_msun,
              T_before_mean, T_after_mean, n_sf_blocked_feedback,
              feedback_wall_ms);
    }
  }

  return SUCCESS;
}
