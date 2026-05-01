/***********************************************************************
/
/  GRID CLASS (BH FEEDBACK, PHASE B THERMAL COUPLING)
/
/  PURPOSE:
/    Compute BH feedback diagnostics, accumulate the thermal feedback
/    reservoir, and deposit burst thermal energy into active-zone gas cells.
/
/  PHASE B LIMITS:
/    - Thermal feedback only.
/    - Kinetic-mode BHs are logged as KINETIC_INACTIVE and do not couple.
/
************************************************************************/

#include <stdio.h>
#include <math.h>
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

struct BHFeedbackKernelCell {
  int index;
  double mass_code;
  BHFeedbackKernelCell(int n, double m) : index(n), mass_code(m) { }
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

/* The checkpointed particle attribute stores the feedback reservoir in code
   energy units. All calculations and BHFDBK log fields use CGS ergs. Keeping
   the storage conversion here prevents 4-byte ParticleAttribute builds from
   overflowing on normal reservoir values around 1e41 erg. */
static double BHFeedbackReservoirStoredToCGS(double reservoir_code,
                                             double energy_units)
{
  if (!isfinite(reservoir_code) || reservoir_code < 0.0 ||
      energy_units <= 0.0)
    return 0.0;
  return reservoir_code * energy_units;
}

static float BHFeedbackReservoirCGSToStored(double reservoir_cgs,
                                            double energy_units)
{
  if (!isfinite(reservoir_cgs) || reservoir_cgs < 0.0 ||
      energy_units <= 0.0)
    return 0.0f;

  const double reservoir_code = reservoir_cgs / energy_units;
  if (!isfinite(reservoir_code) || reservoir_code < 0.0)
    return 0.0f;

  return float(reservoir_code);
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

  /* Ensure the transient SF mask exists for feedback-heated cell blocking.
     Accretion may have already allocated and marked it earlier this pass. */
  if (BaryonField[NumberOfBaryonFields] == NULL) {
    this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
    for (HierarchyEntry *sg = SubgridPointer; sg != NULL; sg = sg->NextGridThisLevel)
      this->ZeroSolutionUnderSubgrid(sg->GridData, ZERO_UNDER_SUBGRID_FIELD);
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL)
    ENZO_FAIL("Error in IdentifyPhysicalQuantities for BH feedback.");

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  std::vector<float> temperature(size, 0.0f);

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
  const double mdot_to_msunyr =
    (TimeUnits > 0.0f) ? mass_units * yr_s / (SolarMass * TimeUnits) : 0.0;
  const double energy_units = mass_units * VelocityUnits * VelocityUnits;
  if (energy_units <= 0.0)
    ENZO_FAIL("BH feedback requires positive code energy units.");
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
  static int warned_realized_missing = FALSE;
  static int warned_realized_anomalous = FALSE;

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
                "reservoir_before=-1 burst_diag=0 "
                "E_deposited=0 p_requested=-1 p_deposited=0 "
                "feedback_kernel_cells=0 feedback_kernel_active_cells=0 "
                "T_before_mean=-1 T_after_mean=-1 n_sf_blocked_feedback=0 "
                "reservoir_after_accum=-1 reservoir_final=-1 burst_occurred=0 "
                "kernel_gas_low=0 kernel_gas_zero=0 dT_mean=0 "
                "sum_deposited=0 deposit_rel_err=0 mass_weight_rel_err=0 "
                "newly_seeded_skip=1 feedback_wall_ms=%.4f "
                "cumul_reservoir_in_cgs=%.8e cumul_reservoir_out_cgs=%.8e "
                "conservation_residual_cgs=%.8e\n",
                cycle_number, level, zred, (long long) ParticleNumber[p], bh_mass_msun,
                feedback_wall_ms, 0.0, 0.0, 0.0);
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

    double mdot_realized_code = 0.0;
    int realized_attr_present = 0;
    int realized_attr_anomalous = 0;
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED &&
        ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED] != NULL) {
      realized_attr_present = 1;
      mdot_realized_code =
        ParticleAttribute[PARTICLE_ATTRIBUTE_BHACCR_LAST_MDOT_REALIZED][p];
      if (!isfinite(mdot_realized_code) || mdot_realized_code < 0.0) {
        realized_attr_anomalous = 1;
        mdot_realized_code = 0.0;
      }
    }

    const int thermal_mode = (f_edd > BHFeedbackModeThreshold) ? 1 : 0;
    const char *feedback_mode = thermal_mode ? "THERMAL" : "KINETIC_INACTIVE";
    const char *feedback_mdot_basis = "requested_actual";
    if (BHFeedbackMethod == 2 && BHFeedbackVerbose >= 1 && !warned_method_two) {
      fprintf(logptr,
              "[BHFDBK_WARN] step=%d level=%d BHFeedbackMethod=2 requested; "
              "kinetic feedback requires Phase C. Running Phase B thermal-only coupling; "
              "low-f_Edd BHs are KINETIC_INACTIVE.\n",
              cycle_number, level);
      warned_method_two = TRUE;
    }

    const double mdot_cgs = mdot_actual_code * mass_rate_to_cgs;
    const double L_feedback = BHAccretionRadiativeEfficiency * mdot_cgs * clight * clight;
    const double mdot_realized_cgs = mdot_realized_code * mass_rate_to_cgs;
    const double L_feedback_realized_basis =
      BHAccretionRadiativeEfficiency * mdot_realized_cgs * clight * clight;

    double reservoir_before = 0.0;
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHFDBK_ENERGY_RESERVOIR) {
      const float reservoir_before_code =
        ParticleAttribute[PARTICLE_ATTRIBUTE_BHFDBK_ENERGY_RESERVOIR][p];
      reservoir_before =
        BHFeedbackReservoirStoredToCGS(reservoir_before_code, energy_units);
    }

    double E_requested = 0.0;
    double p_requested = 0.0;
    double reservoir_after_accum = reservoir_before;
    double reservoir_final = reservoir_before;
    double E_burst = 0.0;
    int burst_diag = 0;
    int burst_occurred = 0;

    if (thermal_mode && mdot_actual_code > 0.0) {
      E_requested = BHFeedbackThermalEfficiency * L_feedback * dt_cgs;
      if (!isfinite(E_requested) || E_requested < 0.0)
        E_requested = 0.0;
      const float E_requested_code =
        BHFeedbackReservoirCGSToStored(E_requested, energy_units);
#ifndef WINDS
      if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BH_CUMUL_RESERVOIR_IN)
        ParticleAttribute[PARTICLE_ATTRIBUTE_BH_CUMUL_RESERVOIR_IN][p] +=
          E_requested_code;
#endif
      reservoir_after_accum = reservoir_before + E_requested;
      reservoir_final = reservoir_after_accum;
      if (reservoir_after_accum >= BHFeedbackMinEnergyBurst) {
        E_burst = reservoir_after_accum;
        reservoir_final = 0.0;
        burst_diag = 1;
      }
    }
    double E_requested_realized_basis = 0.0;
    if (thermal_mode && mdot_realized_code > 0.0) {
      E_requested_realized_basis =
        BHFeedbackThermalEfficiency * L_feedback_realized_basis * dt_cgs;
      if (!isfinite(E_requested_realized_basis) ||
          E_requested_realized_basis < 0.0)
        E_requested_realized_basis = 0.0;
    }

    if (this->ComputeTemperatureField(&temperature[0]) == FAIL)
      ENZO_FAIL("Error in grid->ComputeTemperatureField for BH feedback.");

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
    std::vector<BHFeedbackKernelCell> feedback_cells;
    feedback_cells.reserve(max(1, (2*rcell+1)*(2*rcell+1)*(2*rcell+1)));

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
          if (i < isx || i > iex || j < isy || j > iey || k < isz || k > iez)
            continue;

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
          feedback_cells.push_back(BHFeedbackKernelCell(n, mcell));
        }

    const double T_before_mean = (temp_weight_sum > 0.0) ? (temp_mass_sum / temp_weight_sum) : 0.0;
    double T_after_mean = T_before_mean;
    double E_deposited = 0.0;
    const double p_deposited = 0.0;
    int n_sf_blocked_feedback = 0;
    int kernel_gas_low = 0;
    int kernel_gas_zero = 0;
    double sum_deposited = 0.0;
    double deposit_rel_err = 0.0;
    double mass_weight_rel_err = 0.0;

    if (E_burst > 0.0) {
      if (feedback_cells.empty() || kernel_gas_mass <= 0.0) {
        kernel_gas_zero = 1;
        reservoir_final += E_burst;
        E_burst = 0.0;
        burst_diag = 1;
      } else {
        const double kernel_gas_msun = kernel_gas_mass * mass_to_msun;
        if (kernel_gas_msun < BHFeedbackKernelMassWarnThreshold)
          kernel_gas_low = 1;

        double remaining_energy = E_burst;
        double remaining_mass = kernel_gas_mass;
        for (int c = 0; c < int(feedback_cells.size()); c++) {
          const int n = feedback_cells[c].index;
          const double mcell = feedback_cells[c].mass_code;
          if (mcell <= 0.0)
            continue;

          double dE = 0.0;
          if (c == int(feedback_cells.size()) - 1 || remaining_mass <= 0.0)
            dE = remaining_energy;
          else
            dE = remaining_energy * (mcell / remaining_mass);

          dE = max(0.0, dE);
          remaining_energy -= dE;
          remaining_mass -= mcell;
          sum_deposited += dE;

          if (E_burst > 0.0 && kernel_gas_mass > 0.0) {
            const double expected_fraction = mcell / kernel_gas_mass;
            const double actual_fraction = dE / E_burst;
            const double frac_err = fabs(actual_fraction - expected_fraction);
            if (frac_err > mass_weight_rel_err)
              mass_weight_rel_err = frac_err;
          }

          if (dE <= 0.0)
            continue;

          const double de_specific = dE / (mcell * energy_units);
          if (!isfinite(de_specific) || de_specific < 0.0)
            ENZO_FAIL("BH feedback produced invalid thermal specific-energy increment.");

          if (GENum >= 0 && DualEnergyFormalism) {
            const double newGE = BaryonField[GENum][n] + de_specific;
            BaryonField[GENum][n] = float(newGE);
            BaryonField[TENum][n] = float(newGE);
            if (HydroMethod != Zeus_Hydro) {
              for (int dim = 0; dim < GridRank; dim++)
                BaryonField[TENum][n] +=
                  float(0.5 * BaryonField[Vel1Num+dim][n] *
                        BaryonField[Vel1Num+dim][n]);
            }
          } else {
            BaryonField[TENum][n] += float(de_specific);
          }

          if (BaryonField[NumberOfBaryonFields] != NULL) {
            BaryonField[NumberOfBaryonFields][n] = 1.0f;
            n_sf_blocked_feedback++;
          }
        }

        deposit_rel_err = fabs(sum_deposited - E_burst) / E_burst;
        if (deposit_rel_err > 1.0e-10)
          ENZO_FAIL("BH feedback burst energy conservation violation.");

        E_deposited = E_burst;
        burst_occurred = (E_deposited > 0.0) ? 1 : 0;
        reservoir_final = 0.0;

        if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHFDBK_LAST_REDSHIFT)
          ParticleAttribute[PARTICLE_ATTRIBUTE_BHFDBK_LAST_REDSHIFT][p] = float(zred);

        if (this->ComputeTemperatureField(&temperature[0]) == FAIL)
          ENZO_FAIL("Error in grid->ComputeTemperatureField after BH feedback deposition.");

        double temp_after_mass_sum = 0.0;
        double temp_after_weight_sum = 0.0;
        for (int c = 0; c < int(feedback_cells.size()); c++) {
          const int n = feedback_cells[c].index;
          const double mcell = feedback_cells[c].mass_code;
          temp_after_mass_sum += mcell * temperature[n];
          temp_after_weight_sum += mcell;
        }
        T_after_mean = (temp_after_weight_sum > 0.0) ?
          (temp_after_mass_sum / temp_after_weight_sum) : 0.0;
      }
    }

    if (!isfinite(reservoir_final) || reservoir_final < 0.0)
      reservoir_final = 0.0;

    const float reservoir_after_accum_code =
      BHFeedbackReservoirCGSToStored(reservoir_after_accum, energy_units);
    const float reservoir_final_code =
      BHFeedbackReservoirCGSToStored(reservoir_final, energy_units);

#ifndef WINDS
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BH_CUMUL_RESERVOIR_OUT) {
      const double out_increment =
        double(reservoir_after_accum_code) - double(reservoir_final_code);
      if (out_increment > 0.0)
        ParticleAttribute[PARTICLE_ATTRIBUTE_BH_CUMUL_RESERVOIR_OUT][p] +=
          float(out_increment);
    }
#endif

    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHFDBK_ENERGY_RESERVOIR) {
      ParticleAttribute[PARTICLE_ATTRIBUTE_BHFDBK_ENERGY_RESERVOIR][p] =
        reservoir_final_code;
    }

    double cumul_reservoir_in_cgs = 0.0;
    double cumul_reservoir_out_cgs = 0.0;
    double conservation_residual_cgs = 0.0;
#ifndef WINDS
    if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BH_CUMUL_RESERVOIR_IN &&
        NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BH_CUMUL_RESERVOIR_OUT) {
      const float cumul_in =
        ParticleAttribute[PARTICLE_ATTRIBUTE_BH_CUMUL_RESERVOIR_IN][p];
      const float cumul_out =
        ParticleAttribute[PARTICLE_ATTRIBUTE_BH_CUMUL_RESERVOIR_OUT][p];
      const float reservoir_current =
        (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHFDBK_ENERGY_RESERVOIR) ?
        ParticleAttribute[PARTICLE_ATTRIBUTE_BHFDBK_ENERGY_RESERVOIR][p] :
        reservoir_final_code;
      const double residual_code =
        double(cumul_in) - double(cumul_out) - double(reservoir_current);

      cumul_reservoir_in_cgs =
        BHFeedbackReservoirStoredToCGS(cumul_in, energy_units);
      cumul_reservoir_out_cgs =
        BHFeedbackReservoirStoredToCGS(cumul_out, energy_units);
      conservation_residual_cgs =
        (isfinite(residual_code) && energy_units > 0.0) ?
        residual_code * energy_units : 0.0;

      if (cumul_reservoir_in_cgs > 0.0 &&
          fabs(conservation_residual_cgs) > 1.0e-5 * cumul_reservoir_in_cgs)
        fprintf(stderr, "WARNING: BH feedback conservation residual %.6e "
                "exceeds 1e-5 of cumulative input %.6e\n",
                conservation_residual_cgs, cumul_reservoir_in_cgs);
    }
#endif

    const double kernel_cells_per_radius = kernel_radius_code / cell_width;
    const double dT_mean = T_after_mean - T_before_mean;
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
      if (kernel_gas_low)
        fprintf(logptr,
                "[BHFDBK_WARN] step=%d level=%d bh_id=%lld kernel_gas_msun=%.8e "
                "below BHFeedbackKernelMassWarnThreshold=%.8e Msun; depositing full burst.\n",
                cycle_number, level, (long long) ParticleNumber[p],
                kernel_gas_mass * mass_to_msun,
                double(BHFeedbackKernelMassWarnThreshold));
      if (kernel_gas_zero)
        fprintf(logptr,
                "[BHFDBK_WARN] step=%d level=%d bh_id=%lld kernel_gas_zero=1 "
                "active-zone gas unavailable; returning burst energy to reservoir.\n",
                cycle_number, level, (long long) ParticleNumber[p]);
      if (dT_mean > 1.0e10)
        fprintf(logptr,
                "[BHFDBK_WARN] step=%d level=%d bh_id=%lld dT_mean=%.8e K "
                "extreme feedback temperature jump.\n",
                cycle_number, level, (long long) ParticleNumber[p], dT_mean);
      if (!realized_attr_present && BHAccretionMethod && !warned_realized_missing) {
        fprintf(logptr,
                "[BHFDBK_WARN] step=%d level=%d bh_id=%lld "
                "bhaccr_last_mdot_realized_missing=1; keeping D2-1 feedback "
                "basis requested_actual.\n",
                cycle_number, level, (long long) ParticleNumber[p]);
        warned_realized_missing = TRUE;
      }
      if (realized_attr_anomalous && !warned_realized_anomalous) {
        fprintf(logptr,
                "[BHFDBK_WARN] step=%d level=%d bh_id=%lld "
                "bhaccr_last_mdot_realized_anomalous=1; value reset to zero "
                "for D2-1 diagnostics.\n",
                cycle_number, level, (long long) ParticleNumber[p]);
        warned_realized_anomalous = TRUE;
      }

      const double feedback_wall_ms = 1000.0 * (ReturnWallTime() - t0_all);
      fprintf(logptr,
              "[BHFDBK] step=%d level=%d z=%.8g bh_id=%lld bh_mass=%.8g "
              "feedback_mode=%s f_Edd=%.8e L_feedback=%.8e E_requested=%.8e "
              "feedback_mdot_basis=%s mdot_actual_code=%.15e mdot_realized_code=%.15e "
              "mdot_actual_msunyr=%.8e mdot_realized_msunyr=%.8e "
              "realized_attr_present=%d realized_attr_anomalous=%d "
              "L_feedback_realized_basis=%.8e E_requested_realized_basis=%.8e "
              "reservoir_before=%.8e burst_diag=%d "
              "E_deposited=%.8e p_requested=%.8e p_deposited=%.8e "
              "feedback_kernel_cells=%d feedback_kernel_active_cells=%d "
              "feedback_kernel_gas_msun=%.8e "
              "T_before_mean=%.8e T_after_mean=%.8e n_sf_blocked_feedback=%d "
              "reservoir_after_accum=%.8e reservoir_final=%.8e burst_occurred=%d "
              "kernel_gas_low=%d kernel_gas_zero=%d dT_mean=%.8e "
              "sum_deposited=%.8e deposit_rel_err=%.8e mass_weight_rel_err=%.8e "
              "newly_seeded_skip=0 feedback_wall_ms=%.4f "
              "cumul_reservoir_in_cgs=%.8e cumul_reservoir_out_cgs=%.8e "
              "conservation_residual_cgs=%.8e\n",
              cycle_number, level, zred, (long long) ParticleNumber[p], bh_mass_msun,
              feedback_mode, f_edd, L_feedback, E_requested,
              feedback_mdot_basis, mdot_actual_code, mdot_realized_code,
              mdot_actual_code * mdot_to_msunyr,
              mdot_realized_code * mdot_to_msunyr,
              realized_attr_present, realized_attr_anomalous,
              L_feedback_realized_basis, E_requested_realized_basis,
              reservoir_before, burst_diag,
              E_deposited, p_requested, p_deposited,
              n_kernel_cells, n_kernel_active_cells,
              kernel_gas_mass * mass_to_msun,
              T_before_mean, T_after_mean, n_sf_blocked_feedback,
              reservoir_after_accum, reservoir_final, burst_occurred,
              kernel_gas_low, kernel_gas_zero, dT_mean,
              sum_deposited, deposit_rel_err, mass_weight_rel_err,
              feedback_wall_ms, cumul_reservoir_in_cgs,
              cumul_reservoir_out_cgs, conservation_residual_cgs);
    }
  }

  return SUCCESS;
}
