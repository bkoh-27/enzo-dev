/***********************************************************************
/
/  GRID CLASS (BH DYNAMICAL-FRICTION DIAGNOSTICS, PHASE DF-0)
/
/  PURPOSE:
/    Compute and log a Tremmel-style DM-only collisionless dynamical-
/    friction estimator for BH particles. This phase is diagnostics-only:
/    no force, no particle-state updates, no persistent attributes, and no
/    timestep limiter.
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

static int BHDFIsBHType(int type)
{
  return (type == PARTICLE_TYPE_MBH || type == PARTICLE_TYPE_BLACK_HOLE);
}

static int BHDFIsDarkMatter(int type)
{
  return (type == PARTICLE_TYPE_DARK_MATTER);
}

static int BHDFIsStar(int type)
{
  return (type == PARTICLE_TYPE_STAR);
}

static int BHDFIsNewlySeededThisPass(PINT particle_id)
{
  return (particle_id == INT_UNDEFINED);
}

struct BHDFParticleOrder {
  PINT *ParticleIDs;
  BHDFParticleOrder(PINT *ids) : ParticleIDs(ids) { }
  bool operator()(int a, int b) const {
    if (ParticleIDs[a] < ParticleIDs[b])
      return true;
    if (ParticleIDs[a] > ParticleIDs[b])
      return false;
    return (a < b);
  }
};

static int BHDFPositionCoveredByChild(HierarchyEntry *SubgridPointer,
                                      FLOAT *pos)
{
  for (HierarchyEntry *sub = SubgridPointer; sub != NULL;
       sub = sub->NextGridThisLevel)
    if (sub->GridData != NULL && sub->GridData->PointInGrid(pos))
      return TRUE;

  return FALSE;
}

static int BHDFMassFromAttributes(int num_particle_attributes,
                                  float *particle_attribute[],
                                  int p,
                                  double *full_mass)
{
  if (full_mass != NULL)
    *full_mass = -1.0;

  if (num_particle_attributes <= PARTICLE_ATTRIBUTE_BH_FORMATION_MASS ||
      num_particle_attributes <= PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS ||
      particle_attribute[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS] == NULL ||
      particle_attribute[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS] == NULL)
    return FALSE;

  const double formation_mass =
    particle_attribute[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS][p];
  const double accreted_mass =
    particle_attribute[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS][p];

  if (!isfinite(formation_mass) || formation_mass <= 0.0 ||
      !isfinite(accreted_mass) || accreted_mass < 0.0)
    return FALSE;

  const double mass = formation_mass + accreted_mass;
  if (!isfinite(mass) || mass <= 0.0)
    return FALSE;

  if (full_mass != NULL)
    *full_mass = mass;
  return TRUE;
}

static int BHDFNonAuthoritativeMassRatio(
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

  if (!isfinite(particle_mass) || particle_mass <= 0.0)
    return FALSE;

  double full_mass = -1.0;
  if (!BHDFMassFromAttributes(num_particle_attributes, particle_attribute, p,
                              &full_mass))
    return FALSE;

  const double mass_ratio = full_mass / particle_mass;
  if (expected_full_mass != NULL)
    *expected_full_mass = full_mass;
  if (ratio != NULL)
    *ratio = mass_ratio;

  return (isfinite(mass_ratio) && mass_ratio > 1.5) ? TRUE : FALSE;
}

static double BHDFKernelRadiusCode(float KernelRadiusPhysKpc,
                                   FLOAT time,
                                   float LengthUnits)
{
  double kernel_radius_code = 0.0;

  if (ComovingCoordinates) {
    FLOAT a = 1.0, dadt = 0.0;
    CosmologyComputeExpansionFactor(time, &a, &dadt);
    float a_phys = float(a / (1.0 + InitialRedshift));

    float h_param = HubbleConstantNow;
    if (h_param > 10.0f)
      h_param *= 0.01f;
    if (h_param <= 0.0f)
      h_param = 1.0f;
    if (a_phys <= 0.0f)
      a_phys = 1.0f;

    const float box_kpch = ComovingBoxSize * 1000.0f;
    const float kernel_comoving_kpch = KernelRadiusPhysKpc * h_param / a_phys;
    if (box_kpch > 0.0f)
      kernel_radius_code = kernel_comoving_kpch / box_kpch;
  } else {
    const float box_kpc = LengthUnits / kpc_cm;
    if (box_kpc > 0.0f)
      kernel_radius_code = KernelRadiusPhysKpc / box_kpc;
  }

  return kernel_radius_code;
}

struct BHDFRow {
  int processor;
  int invocation_seq;
  int grid_id;
  int step;
  int level;
  int grid_rank;
  double z;
  long long bh_id;
  int particle_index;
  double bh_mass_msun;
  int active_reposition_enabled;
  int authoritative;
  double bh_pos[MAX_DIMENSION];
  double bh_vel_kms[MAX_DIMENSION];
  double kernel_radius_code;
  double kernel_radius_phys_kpc;
  double kernel_radius_over_dx;
  int kernel_complete;
  int kernel_truncated_by_grid;
  int neighbor_particles_total;
  int n_dm_in_kernel;
  int n_star_excluded;
  int n_bh_excluded;
  int n_slow;
  double mass_slow_code;
  double rho_slow_cgs;
  double v_bh_rel_kms;
  double v_com_kms[MAX_DIMENSION];
  double sigma_1d_kms;
  double b_min_phys_kpc;
  double b_max_phys_kpc;
  double ln_lambda;
  int ln_lambda_zero;
  int ln_lambda_large;
  double a_df_cgs[MAX_DIMENSION];
  double a_df_mag_cgs;
  double dt_df_myr;
  double velocity_kick_cap_fraction;
  int cap_triggered;
  int under_resolved;
  int rejection_reason;
  double df_wall_ms;
};

static void BHDFInitRow(BHDFRow *row)
{
  row->processor = 0;
  row->invocation_seq = 0;
  row->grid_id = 0;
  row->step = 0;
  row->level = 0;
  row->grid_rank = 0;
  row->z = 0.0;
  row->bh_id = -1;
  row->particle_index = -1;
  row->bh_mass_msun = -1.0;
  row->active_reposition_enabled = 0;
  row->authoritative = 1;
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    row->bh_pos[dim] = 0.0;
    row->bh_vel_kms[dim] = 0.0;
    row->v_com_kms[dim] = 0.0;
    row->a_df_cgs[dim] = 0.0;
  }
  row->kernel_radius_code = 0.0;
  row->kernel_radius_phys_kpc = 0.0;
  row->kernel_radius_over_dx = 0.0;
  row->kernel_complete = 1;
  row->kernel_truncated_by_grid = 0;
  row->neighbor_particles_total = 0;
  row->n_dm_in_kernel = 0;
  row->n_star_excluded = 0;
  row->n_bh_excluded = 0;
  row->n_slow = 0;
  row->mass_slow_code = 0.0;
  row->rho_slow_cgs = 0.0;
  row->v_bh_rel_kms = 0.0;
  row->sigma_1d_kms = 0.0;
  row->b_min_phys_kpc = -1.0;
  row->b_max_phys_kpc = 0.0;
  row->ln_lambda = 0.0;
  row->ln_lambda_zero = 1;
  row->ln_lambda_large = 0;
  row->a_df_mag_cgs = 0.0;
  row->dt_df_myr = -1.0;
  row->velocity_kick_cap_fraction = 0.0;
  row->cap_triggered = 0;
  row->under_resolved = 0;
  row->rejection_reason = 0;
  row->df_wall_ms = 0.0;
}

static void BHDFEmitDetail(FILE *logptr,
                           const BHDFRow &row,
                           int detail_rows_emitted,
                           int detail_rows_complete)
{
  fprintf(logptr,
          "[BHDF] schema_version=1 "
          "processor=%d invocation_seq=%d grid_id=%d step=%d level=%d "
          "grid_rank=%d z=%.8g "
          "bh_id=%lld particle_index=%d bh_mass_msun=%.8e "
          "active_reposition_enabled=%d mode=diagnostic applied=0 "
          "authoritative=%d "
          "bh_pos_x=%.15g bh_pos_y=%.15g bh_pos_z=%.15g "
          "bh_vel_x_kms=%.8g bh_vel_y_kms=%.8g bh_vel_z_kms=%.8g "
          "kernel_radius_code=%.8g kernel_radius_phys_kpc=%.8g "
          "kernel_radius_over_dx=%.8g "
          "kernel_complete=%d kernel_truncated_by_grid=%d "
          "neighbor_particles_total=%d n_dm_in_kernel=%d "
          "n_star_excluded=%d n_bh_excluded=%d "
          "n_slow=%d mass_slow_code=%.8e rho_slow_cgs=%.8e "
          "v_bh_rel_kms=%.8g "
          "v_CoM_x_kms=%.8g v_CoM_y_kms=%.8g v_CoM_z_kms=%.8g "
          "sigma_1d_kms=%.8g "
          "b_min_phys_kpc=%.8g b_max_phys_kpc=%.8g "
          "ln_lambda=%.8g ln_lambda_zero=%d ln_lambda_large=%d "
          "a_DF_x_cgs=%.8e a_DF_y_cgs=%.8e a_DF_z_cgs=%.8e "
          "a_DF_mag_cgs=%.8e "
          "a_grav_available=0 a_grav_cgs=-1.0 "
          "a_DF_over_a_grav_available=0 a_DF_over_a_grav=-1.0 "
          "dt_DF_Myr=%.8g velocity_kick_cap_fraction=%.8g "
          "cap_triggered=%d under_resolved=%d rejection_reason=%d "
          "df_wall_ms=%.4f detail_rows_emitted=%d "
          "detail_rows_complete=%d\n",
          row.processor, row.invocation_seq, row.grid_id, row.step, row.level,
          row.grid_rank, row.z,
          row.bh_id, row.particle_index, row.bh_mass_msun,
          row.active_reposition_enabled, row.authoritative,
          row.bh_pos[0], row.bh_pos[1], row.bh_pos[2],
          row.bh_vel_kms[0], row.bh_vel_kms[1], row.bh_vel_kms[2],
          row.kernel_radius_code, row.kernel_radius_phys_kpc,
          row.kernel_radius_over_dx,
          row.kernel_complete, row.kernel_truncated_by_grid,
          row.neighbor_particles_total, row.n_dm_in_kernel,
          row.n_star_excluded, row.n_bh_excluded,
          row.n_slow, row.mass_slow_code, row.rho_slow_cgs,
          row.v_bh_rel_kms,
          row.v_com_kms[0], row.v_com_kms[1], row.v_com_kms[2],
          row.sigma_1d_kms,
          row.b_min_phys_kpc, row.b_max_phys_kpc,
          row.ln_lambda, row.ln_lambda_zero, row.ln_lambda_large,
          row.a_df_cgs[0], row.a_df_cgs[1], row.a_df_cgs[2],
          row.a_df_mag_cgs,
          row.dt_df_myr, row.velocity_kick_cap_fraction,
          row.cap_triggered, row.under_resolved, row.rejection_reason,
          row.df_wall_ms, detail_rows_emitted, detail_rows_complete);
}

static void BHDFEmitGuardWarning(FILE *logptr, const BHDFRow &row)
{
  fprintf(logptr,
          "[BHDF_WARN] schema_version=1 processor=%d invocation_seq=%d "
          "grid_id=%d step=%d level=%d bh_id=%lld warning=guard_skip "
          "rejection_reason=%d authoritative=%d\n",
          row.processor, row.invocation_seq, row.grid_id, row.step, row.level,
          row.bh_id, row.rejection_reason, row.authoritative);
}

struct BHDFCandidateResult {
  int n_dm_in_kernel;
  int n_slow;
  double rho_slow_cgs;
  double v_rel_kms;
  double sigma_1d_kms;
  double b_min_phys_kpc;
  double b_max_phys_kpc;
  double ln_lambda;
  double candidate_a_DF_cgs[MAX_DIMENSION];
  double candidate_a_DF_mag_cgs;
  double dt_DF_Myr;
  double dt_level_over_dt_DF;
  double requested_delta_v_code;
  double cap_fraction_requested;
  int kernel_complete;
  int kernel_truncated_by_grid;
  int computed;
};

static void BHDFInitCandidateResult(BHDFCandidateResult *candidate,
                                    double kernel_radius_phys_kpc)
{
  candidate->n_dm_in_kernel = 0;
  candidate->n_slow = 0;
  candidate->rho_slow_cgs = 0.0;
  candidate->v_rel_kms = 0.0;
  candidate->sigma_1d_kms = 0.0;
  candidate->b_min_phys_kpc = 0.0;
  candidate->b_max_phys_kpc =
    (isfinite(kernel_radius_phys_kpc) && kernel_radius_phys_kpc > 0.0) ?
    kernel_radius_phys_kpc : 0.0;
  candidate->ln_lambda = 0.0;
  for (int dim = 0; dim < MAX_DIMENSION; dim++)
    candidate->candidate_a_DF_cgs[dim] = 0.0;
  candidate->candidate_a_DF_mag_cgs = 0.0;
  candidate->dt_DF_Myr = 0.0;
  candidate->dt_level_over_dt_DF = 0.0;
  candidate->requested_delta_v_code = 0.0;
  candidate->cap_fraction_requested = 0.0;
  candidate->kernel_complete = 1;
  candidate->kernel_truncated_by_grid = 0;
  candidate->computed = 0;
}

static BHDFCandidateResult BHDFDefaultCandidateWithKernel(
  double kernel_radius_phys_kpc, int kernel_complete)
{
  BHDFCandidateResult candidate;
  BHDFInitCandidateResult(&candidate, kernel_radius_phys_kpc);
  candidate.kernel_complete = kernel_complete ? 1 : 0;
  candidate.kernel_truncated_by_grid = candidate.kernel_complete ? 0 : 1;
  return candidate;
}

static BHDFCandidateResult BHDFComputeCandidate(
  int number_of_particles, int *particle_type, float *particle_mass,
  FLOAT *particle_position[], float *particle_velocity[],
  int bh_index, const double bh_pos[], double bh_mass_code,
  double kernel_radius_code, double kernel_radius_phys_kpc,
  double dt_level_myr, double mass_units, double length_units,
  double velocity_units, double cell_volume_code,
  const FLOAT *grid_left_edge, const FLOAT *grid_right_edge)
{
  BHDFCandidateResult candidate;
  BHDFInitCandidateResult(&candidate, kernel_radius_phys_kpc);

  if (number_of_particles <= 0 || particle_type == NULL ||
      particle_mass == NULL || particle_position == NULL ||
      particle_velocity == NULL || bh_index < 0 ||
      bh_index >= number_of_particles || !isfinite(bh_mass_code) ||
      bh_mass_code <= 0.0 || !isfinite(kernel_radius_code) ||
      kernel_radius_code <= 0.0 || !isfinite(kernel_radius_phys_kpc) ||
      kernel_radius_phys_kpc <= 0.0 || !isfinite(mass_units) ||
      mass_units <= 0.0 || !isfinite(length_units) ||
      length_units <= 0.0 || !isfinite(velocity_units) ||
      velocity_units <= 0.0 || !isfinite(cell_volume_code) ||
      cell_volume_code <= 0.0)
    return candidate;

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    if (particle_position[dim] == NULL || particle_velocity[dim] == NULL)
      return candidate;
  }

  const double kernel_r2 = kernel_radius_code * kernel_radius_code;
  const double radius_cgs = kernel_radius_code * length_units;
  const double kernel_volume_cgs =
    (4.0 * M_PI / 3.0) * radius_cgs * radius_cgs * radius_cgs;
  const double b_max_phys_cm = candidate.b_max_phys_kpc * kpc_cm;

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    if (grid_left_edge != NULL && grid_right_edge != NULL &&
        (bh_pos[dim] - kernel_radius_code < double(grid_left_edge[dim]) ||
         bh_pos[dim] + kernel_radius_code > double(grid_right_edge[dim])))
      candidate.kernel_complete = 0;
  }
  candidate.kernel_truncated_by_grid =
    candidate.kernel_complete ? 0 : 1;
  const int kernel_complete = candidate.kernel_complete;

  double total_dm_mass = 0.0;
  double mass_slow_code = 0.0;
  double com_vel_code[MAX_DIMENSION] = {0.0, 0.0, 0.0};
  std::vector<int> dm_neighbors;
  dm_neighbors.reserve(number_of_particles);

  for (int q = 0; q < number_of_particles; q++) {
    double dr2 = 0.0;
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      const double dr = double(particle_position[dim][q]) - bh_pos[dim];
      dr2 += dr * dr;
    }
    if (!isfinite(dr2))
      continue;
    if (dr2 > kernel_r2)
      continue;

    if (!BHDFIsDarkMatter(particle_type[q]))
      continue;

    const double dm_mass = double(particle_mass[q]);
    if (!isfinite(dm_mass) || dm_mass <= 0.0)
      continue;

    candidate.n_dm_in_kernel++;
    total_dm_mass += dm_mass;
    dm_neighbors.push_back(q);
    for (int dim = 0; dim < MAX_DIMENSION; dim++)
      com_vel_code[dim] += dm_mass * double(particle_velocity[dim][q]);
  }

  if (total_dm_mass <= 0.0)
    return BHDFDefaultCandidateWithKernel(
      kernel_radius_phys_kpc, kernel_complete);

  for (int dim = 0; dim < MAX_DIMENSION; dim++)
    com_vel_code[dim] /= total_dm_mass;

  double v_rel_code_vec[MAX_DIMENSION] = {0.0, 0.0, 0.0};
  double v_rel_code2 = 0.0;
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    v_rel_code_vec[dim] =
      double(particle_velocity[dim][bh_index]) - com_vel_code[dim];
    v_rel_code2 += v_rel_code_vec[dim] * v_rel_code_vec[dim];
  }

  if (!isfinite(v_rel_code2) || v_rel_code2 <= 0.0)
    return BHDFDefaultCandidateWithKernel(
      kernel_radius_phys_kpc, kernel_complete);

  const double v_rel_code = sqrt(v_rel_code2);
  const double v_rel_cgs = v_rel_code * velocity_units;
  candidate.v_rel_kms = v_rel_cgs / 1.0e5;

  double sigma2_num = 0.0;
  for (int ni = 0; ni < int(dm_neighbors.size()); ni++) {
    const int q = dm_neighbors[ni];
    const double dm_mass = double(particle_mass[q]);
    double dv2 = 0.0;
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      const double dv = double(particle_velocity[dim][q]) - com_vel_code[dim];
      dv2 += dv * dv;
    }
    if (!isfinite(dv2))
      return BHDFDefaultCandidateWithKernel(
        kernel_radius_phys_kpc, kernel_complete);
    sigma2_num += dm_mass * dv2;
    if (sqrt(dv2) < v_rel_code) {
      candidate.n_slow++;
      mass_slow_code += dm_mass;
    }
  }

  if (!isfinite(sigma2_num) || sigma2_num < 0.0)
    return BHDFDefaultCandidateWithKernel(
      kernel_radius_phys_kpc, kernel_complete);

  const double sigma_1d_code = sqrt(sigma2_num / (3.0 * total_dm_mass));
  if (!isfinite(sigma_1d_code))
    return BHDFDefaultCandidateWithKernel(
      kernel_radius_phys_kpc, kernel_complete);
  candidate.sigma_1d_kms = sigma_1d_code * velocity_units / 1.0e5;

  // MASSAUDIT0/MASSFIX0: native Enzo DM ParticleMass is density-like.
  // Convert the stored slow-particle sum to true code mass only at this
  // endpoint so mass-weighted velocity, v_rel, and sigma paths remain
  // unchanged. BH masses are absolute and are not scaled here.
  const double mass_slow_code_true = mass_slow_code * cell_volume_code;
  if (kernel_volume_cgs > 0.0 && isfinite(mass_slow_code_true))
    candidate.rho_slow_cgs =
      mass_slow_code_true * mass_units / kernel_volume_cgs;

  const double bh_mass_cgs = bh_mass_code * mass_units;
  const double sigma_1d_cgs = sigma_1d_code * velocity_units;
  const double bmin_denom =
    v_rel_cgs * v_rel_cgs + sigma_1d_cgs * sigma_1d_cgs;
  if (bmin_denom > 0.0 && isfinite(bh_mass_cgs) && bh_mass_cgs > 0.0) {
    const double b_min_phys_cm = GravConst * bh_mass_cgs / bmin_denom;
    candidate.b_min_phys_kpc = b_min_phys_cm / kpc_cm;
    if (b_min_phys_cm > 0.0 && b_min_phys_cm < b_max_phys_cm)
      candidate.ln_lambda = log(b_max_phys_cm / b_min_phys_cm);
  }

  if (!isfinite(candidate.ln_lambda) || candidate.ln_lambda <= 0.0 ||
      !isfinite(candidate.rho_slow_cgs) || candidate.rho_slow_cgs <= 0.0 ||
      v_rel_cgs <= 0.0)
    return BHDFDefaultCandidateWithKernel(
      kernel_radius_phys_kpc, kernel_complete);

  const double prefactor =
    -4.0 * M_PI * GravConst * GravConst * bh_mass_cgs *
    candidate.rho_slow_cgs * candidate.ln_lambda /
    (v_rel_cgs * v_rel_cgs * v_rel_cgs);
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    const double v_rel_dim_cgs = v_rel_code_vec[dim] * velocity_units;
    candidate.candidate_a_DF_cgs[dim] = prefactor * v_rel_dim_cgs;
    candidate.candidate_a_DF_mag_cgs +=
      candidate.candidate_a_DF_cgs[dim] *
      candidate.candidate_a_DF_cgs[dim];
  }
  candidate.candidate_a_DF_mag_cgs =
    sqrt(candidate.candidate_a_DF_mag_cgs);

  if (!isfinite(candidate.candidate_a_DF_mag_cgs) ||
      candidate.candidate_a_DF_mag_cgs <= 0.0)
    return BHDFDefaultCandidateWithKernel(
      kernel_radius_phys_kpc, kernel_complete);

  candidate.dt_DF_Myr =
    (v_rel_cgs / candidate.candidate_a_DF_mag_cgs) / Myr_s;
  if (isfinite(candidate.dt_DF_Myr) && candidate.dt_DF_Myr > 0.0)
    candidate.dt_level_over_dt_DF = dt_level_myr / candidate.dt_DF_Myr;

  candidate.requested_delta_v_code =
    (candidate.candidate_a_DF_mag_cgs / velocity_units) *
    (dt_level_myr * Myr_s);
  if (v_rel_code > 0.0)
    candidate.cap_fraction_requested =
      candidate.requested_delta_v_code / v_rel_code;

  if (!isfinite(candidate.dt_DF_Myr) ||
      !isfinite(candidate.dt_level_over_dt_DF) ||
      !isfinite(candidate.requested_delta_v_code) ||
      candidate.requested_delta_v_code < 0.0 ||
      !isfinite(candidate.cap_fraction_requested) ||
      candidate.cap_fraction_requested < 0.0) {
    return BHDFDefaultCandidateWithKernel(
      kernel_radius_phys_kpc, kernel_complete);
  }

  candidate.computed = 1;
  return candidate;
}

static void BHDFEmitActiveDryRunDetail(
  FILE *logptr, int processor, int invocation_seq, int grid_id, int step,
  int level, long long bh_id, int particle_index, double bh_mass_msun,
  const double bh_pos[], const double bh_vel_kms[],
  double kernel_radius_code, double kernel_radius_phys_kpc,
  double kernel_radius_over_dx, double dt_level_myr,
  const BHDFCandidateResult &candidate)
{
  fprintf(logptr,
          "[BHDF_ACTIVE] schema_version=%d mode=active_dry_run "
          "active_method=2 active_policy=dry_run "
          "processor=%d invocation_seq=%d grid_id=%d step=%d level=%d "
          "bh_id=%lld particle_index=%d bh_mass_msun=%.8e "
          "bh_pos_x=%.15g bh_pos_y=%.15g bh_pos_z=%.15g "
          "bh_vel_x_kms=%.8g bh_vel_y_kms=%.8g bh_vel_z_kms=%.8g "
          "kernel_radius_code=%.8g kernel_radius_phys_kpc=%.8g "
          "kernel_radius_over_dx=%.8g kernel_complete=%d "
          "kernel_truncated_by_grid=%d n_dm_in_kernel=%d n_slow=%d "
          "rho_slow_cgs=%.8e v_rel_kms=%.8g sigma_1d_kms=%.8g "
          "b_min_phys_kpc=%.8g b_max_phys_kpc=%.8g ln_lambda=%.8g "
          "candidate_a_DF_x_cgs=%.8e candidate_a_DF_y_cgs=%.8e "
          "candidate_a_DF_z_cgs=%.8e candidate_a_DF_mag_cgs=%.8e "
          "limited_a_DF_x_cgs=0 limited_a_DF_y_cgs=0 "
          "limited_a_DF_z_cgs=0 limited_a_DF_mag_cgs=0 "
          "dt_level_Myr=%.8g dt_DF_Myr=%.8g dt_level_over_dt_DF=%.8g "
          "requested_delta_v_code=%.8g limited_delta_v_code=0 "
          "applied_delta_v_code=0 cap_fraction_requested=%.8g "
          "cap_fraction_applied=0 per_substep_cap=%.8g "
          "per_level_cap=%.8g subcycle_count=0 max_subcycles=%d "
          "suppression_reason_code=1 suppression_reason_name=diagnostics_only "
          "active_applied=0 kick_applied=0 local_gravity_available=1 "
          "local_gravity_magnitude_cgs=0 a_DF_over_a_grav=0 "
          "local_gravity_guard_passed=0 momentum_policy=bh_only "
          "conservation_bookkeeping_available=0 "
          "delta_p_bh_x_code=0 delta_p_bh_y_code=0 delta_p_bh_z_code=0 "
          "delta_p_wake_x_code=0 delta_p_wake_y_code=0 "
          "delta_p_wake_z_code=0 active_validation_flags=dry_run_no_physics\n",
          BHDynamicalFrictionActiveSchemaVersion,
          processor, invocation_seq, grid_id, step, level,
          bh_id, particle_index, bh_mass_msun,
          bh_pos[0], bh_pos[1], bh_pos[2],
          bh_vel_kms[0], bh_vel_kms[1], bh_vel_kms[2],
          kernel_radius_code, kernel_radius_phys_kpc,
          kernel_radius_over_dx, candidate.kernel_complete,
          candidate.kernel_truncated_by_grid, candidate.n_dm_in_kernel,
          candidate.n_slow, candidate.rho_slow_cgs, candidate.v_rel_kms,
          candidate.sigma_1d_kms, candidate.b_min_phys_kpc,
          candidate.b_max_phys_kpc, candidate.ln_lambda,
          candidate.candidate_a_DF_cgs[0],
          candidate.candidate_a_DF_cgs[1],
          candidate.candidate_a_DF_cgs[2],
          candidate.candidate_a_DF_mag_cgs,
          dt_level_myr, candidate.dt_DF_Myr,
          candidate.dt_level_over_dt_DF,
          candidate.requested_delta_v_code,
          candidate.cap_fraction_requested,
          BHDynamicalFrictionActiveMaxKickFraction,
          BHDynamicalFrictionActiveMaxLevelKickFraction,
          BHDynamicalFrictionActiveMaxSubcycles);
}

static void BHDFEmitActiveDryRunSummary(
  FILE *logptr, int processor, int invocation_seq, int grid_id, int step,
  int level, int eligible_bh, int detail_rows_emitted,
  int detail_rows_complete)
{
  fprintf(logptr,
          "[BHDF_ACTIVE_SUMMARY] schema_version=%d "
          "processor=%d invocation_seq=%d grid_id=%d step=%d level=%d "
          "eligible_bh=%d detail_rows_emitted=%d active_rows=%d "
          "applied_kicks=0 suppressed_rows=%d cap_limited_rows=0 "
          "subcycled_rows=0 local_gravity_failed_rows=0 "
          "under_resolved_rows=0 no_slow_rows=0 kernel_incomplete_rows=0 "
          "low_vrel_suppressed_rows=0 lnLambda_suppressed_rows=0 "
          "dtDF_suppressed_rows=0 max_subcycle_suppressed_rows=0 "
          "reposition_conflict_rows=0 newly_seeded_suppressed_rows=0 "
          "non_authoritative_mass_rows=0 invalid_units_rows=0 "
          "momentum_policy_unsupported_rows=0 "
          "invalid_candidate_acceleration_rows=0 diagnostics_only_rows=%d "
          "cap_fraction_extreme_rows=0 local_gravity_ratio_exceeded_rows=0 "
          "summary_accounting_passed=1 detail_rows_complete=%d "
          "active_rows_consistent=1 applied_kicks_consistent=1\n",
          BHDynamicalFrictionActiveSchemaVersion,
          processor, invocation_seq, grid_id, step, level,
          eligible_bh, detail_rows_emitted, detail_rows_emitted,
          detail_rows_emitted, detail_rows_emitted, detail_rows_complete);
}

static int BHDFActiveDryRunSkeleton(
  FILE *logptr, int processor, int invocation_seq, int grid_id, int step,
  int level, int number_of_particles, int *particle_type,
  PINT *particle_number, float *particle_mass,
  FLOAT *particle_position[], float *particle_velocity[],
  int number_of_particle_attributes, float *particle_attribute[],
  double mass_units, double length_units, double velocity_units,
  double cell_volume_code, const FLOAT *grid_left_edge,
  const FLOAT *grid_right_edge, double kernel_radius_code,
  double kernel_radius_phys_kpc, double kernel_radius_over_dx,
  double dt_level_myr)
{
  if (logptr == NULL)
    logptr = stdout;
  if (BHDynamicalFrictionActiveVerbose <= 0)
    return SUCCESS;

  std::vector<int> bh_particles;
  if (number_of_particles > 0 && particle_type != NULL) {
    bh_particles.reserve(number_of_particles);
    for (int p = 0; p < number_of_particles; p++)
      if (BHDFIsBHType(particle_type[p]))
        bh_particles.push_back(p);
  }

  if (particle_number != NULL)
    std::sort(bh_particles.begin(), bh_particles.end(),
              BHDFParticleOrder(particle_number));

  const int eligible_bh = int(bh_particles.size());
  int detail_rows_emitted = 0;
  const double mass_to_msun =
    (isfinite(mass_units) && mass_units > 0.0) ? mass_units / SolarMass : 0.0;
  const double velocity_to_kms =
    (isfinite(velocity_units) && velocity_units > 0.0) ?
    velocity_units / 1.0e5 : 0.0;

  if (BHDynamicalFrictionActiveVerbose >= 2) {
    for (int bp = 0; bp < eligible_bh; bp++) {
      const int p = bh_particles[bp];
      double bh_pos[MAX_DIMENSION] = {0.0, 0.0, 0.0};
      double bh_vel_kms[MAX_DIMENSION] = {0.0, 0.0, 0.0};
      BHDFCandidateResult candidate;
      BHDFInitCandidateResult(&candidate, kernel_radius_phys_kpc);

      for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        if (particle_position[dim] != NULL)
          bh_pos[dim] = double(particle_position[dim][p]);
        if (particle_velocity[dim] != NULL)
          bh_vel_kms[dim] =
            double(particle_velocity[dim][p]) * velocity_to_kms;
      }

      double bh_mass_code = 0.0;
      if (particle_mass != NULL &&
          isfinite(double(particle_mass[p])) && double(particle_mass[p]) > 0.0)
        bh_mass_code = double(particle_mass[p]);

      double attribute_mass = -1.0;
      if (BHDFMassFromAttributes(number_of_particle_attributes,
                                 particle_attribute, p, &attribute_mass))
        bh_mass_code = attribute_mass;

      double bh_mass_msun = 0.0;
      if (isfinite(bh_mass_code) && bh_mass_code > 0.0 &&
          isfinite(mass_to_msun) && mass_to_msun > 0.0)
        bh_mass_msun = bh_mass_code * mass_to_msun;

      const long long bh_id = (particle_number != NULL) ?
        (long long) particle_number[p] : -1;
      const int newly_seeded =
        (particle_number != NULL &&
         BHDFIsNewlySeededThisPass(particle_number[p])) ? TRUE : FALSE;

      int point_in_grid = TRUE;
      if (grid_left_edge != NULL && grid_right_edge != NULL) {
        for (int dim = 0; dim < MAX_DIMENSION; dim++) {
          if (bh_pos[dim] < double(grid_left_edge[dim]) ||
              bh_pos[dim] > double(grid_right_edge[dim]))
            point_in_grid = FALSE;
        }
      }

      const int invalid_mass =
        (!isfinite(bh_mass_code) || bh_mass_code <= 0.0 ||
         particle_mass == NULL || !isfinite(double(particle_mass[p])) ||
         double(particle_mass[p]) <= 0.0) ? TRUE : FALSE;

      if (!newly_seeded && point_in_grid && !invalid_mass) {
        candidate = BHDFComputeCandidate(
          number_of_particles, particle_type, particle_mass,
          particle_position, particle_velocity, p, bh_pos, bh_mass_code,
          kernel_radius_code, kernel_radius_phys_kpc, dt_level_myr,
          mass_units, length_units, velocity_units, cell_volume_code,
          grid_left_edge, grid_right_edge);
      }

      BHDFEmitActiveDryRunDetail(
        logptr, processor, invocation_seq, grid_id, step, level,
        bh_id, p, bh_mass_msun, bh_pos, bh_vel_kms,
        kernel_radius_code, kernel_radius_phys_kpc,
        kernel_radius_over_dx, dt_level_myr, candidate);
      detail_rows_emitted++;
    }
  }

  const int detail_rows_complete =
    (BHDynamicalFrictionActiveVerbose >= 2 || eligible_bh == 0) ? 1 : 0;
  BHDFEmitActiveDryRunSummary(
    logptr, processor, invocation_seq, grid_id, step, level,
    eligible_bh, detail_rows_emitted, detail_rows_complete);

  return SUCCESS;
}

static double BHDFMethod3CleanReal(double value)
{
  if (!isfinite(value))
    return 0.0;
  return (value == 0.0) ? 0.0 : value;
}

static const char *BHDFMethod3ReasonName(int code)
{
  switch (code) {
  case 0: return "valid";
  case 1: return "shadow_no_apply";
  case 2: return "under_resolved";
  case 3: return "no_slow_particles";
  case 4: return "kernel_truncated";
  case 5: return "ln_lambda_zero";
  case 6: return "guard_skip";
  case 7: return "invalid_units";
  case 8: return "invalid_bh_mass";
  case 9: return "non_authoritative_mass";
  case 10: return "newly_seeded";
  case 11: return "outside_grid";
  case 12: return "child_covered";
  case 13: return "reposition_conflict";
  case 14: return "invalid_candidate_acceleration";
  case 15: return "kernel_incomplete";
  case 16: return "ln_lambda_invalid";
  case 17: return "low_vrel";
  case 18: return "low_sigma";
  case 19: return "dt_df_extreme";
  case 20: return "local_gravity_unavailable";
  case 21: return "local_gravity_ratio_exceeded";
  case 22: return "momentum_policy_unsupported";
  case 23: return "cap_fraction_extreme";
  case 24: return "max_subcycles_exceeded";
  case 25: return "apply_not_authorized";
  default: return "invalid_candidate_acceleration";
  }
}

static const char *BHDFMethod3PrimaryReasonPrecedence()
{
  return "valid|shadow_no_apply|under_resolved|no_slow_particles|"
    "kernel_truncated|ln_lambda_zero|guard_skip|invalid_units|"
    "invalid_bh_mass|non_authoritative_mass|newly_seeded|outside_grid|"
    "child_covered|reposition_conflict|invalid_candidate_acceleration|"
    "kernel_incomplete|ln_lambda_invalid|low_vrel|low_sigma|"
    "dt_df_extreme|local_gravity_unavailable|local_gravity_ratio_exceeded|"
    "momentum_policy_unsupported|cap_fraction_extreme|"
    "max_subcycles_exceeded|apply_not_authorized";
}

struct BHDFMethod3Summary {
  int eligible_bh;
  int detail_rows_emitted;
  int cap_limited_rows;
  int reason_counts[26];
};

static void BHDFMethod3InitSummary(BHDFMethod3Summary *summary)
{
  summary->eligible_bh = 0;
  summary->detail_rows_emitted = 0;
  summary->cap_limited_rows = 0;
  for (int i = 0; i < 26; i++)
    summary->reason_counts[i] = 0;
}

static void BHDFMethod3NeighborPartition(
  int number_of_particles, int *particle_type, FLOAT *particle_position[],
  const double bh_pos[], double kernel_radius_code,
  int *n_star_excluded, int *n_bh_excluded)
{
  *n_star_excluded = 0;
  *n_bh_excluded = 0;
  if (number_of_particles <= 0 || particle_type == NULL ||
      particle_position == NULL || !isfinite(kernel_radius_code) ||
      kernel_radius_code <= 0.0)
    return;

  for (int dim = 0; dim < MAX_DIMENSION; dim++)
    if (particle_position[dim] == NULL)
      return;

  const double kernel_r2 = kernel_radius_code * kernel_radius_code;
  for (int q = 0; q < number_of_particles; q++) {
    double dr2 = 0.0;
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      const double dr = double(particle_position[dim][q]) - bh_pos[dim];
      dr2 += dr * dr;
    }
    if (!isfinite(dr2) || dr2 > kernel_r2)
      continue;
    if (BHDFIsStar(particle_type[q]))
      (*n_star_excluded)++;
    else if (BHDFIsBHType(particle_type[q]))
      (*n_bh_excluded)++;
  }
}

static void BHDFMethod3VRelVectorKms(
  int number_of_particles, int *particle_type, float *particle_mass,
  FLOAT *particle_position[], float *particle_velocity[], int bh_index,
  const double bh_pos[], double kernel_radius_code, double velocity_units,
  const BHDFCandidateResult &candidate, double v_rel_kms_vec[])
{
  for (int dim = 0; dim < MAX_DIMENSION; dim++)
    v_rel_kms_vec[dim] = 0.0;
  if (!candidate.computed || number_of_particles <= 0 ||
      particle_type == NULL || particle_mass == NULL ||
      particle_position == NULL || particle_velocity == NULL ||
      bh_index < 0 || bh_index >= number_of_particles ||
      !isfinite(kernel_radius_code) || kernel_radius_code <= 0.0 ||
      !isfinite(velocity_units) || velocity_units <= 0.0)
    return;
  for (int dim = 0; dim < MAX_DIMENSION; dim++)
    if (particle_position[dim] == NULL || particle_velocity[dim] == NULL)
      return;

  const double kernel_r2 = kernel_radius_code * kernel_radius_code;
  double total_dm_mass = 0.0;
  double com_vel_code[MAX_DIMENSION] = {0.0, 0.0, 0.0};
  for (int q = 0; q < number_of_particles; q++) {
    double dr2 = 0.0;
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      const double dr = double(particle_position[dim][q]) - bh_pos[dim];
      dr2 += dr * dr;
    }
    if (!isfinite(dr2) || dr2 > kernel_r2)
      continue;
    if (!BHDFIsDarkMatter(particle_type[q]))
      continue;
    const double dm_mass = double(particle_mass[q]);
    if (!isfinite(dm_mass) || dm_mass <= 0.0)
      continue;
    total_dm_mass += dm_mass;
    for (int dim = 0; dim < MAX_DIMENSION; dim++)
      com_vel_code[dim] += dm_mass * double(particle_velocity[dim][q]);
  }
  if (total_dm_mass <= 0.0)
    return;
  for (int dim = 0; dim < MAX_DIMENSION; dim++)
    com_vel_code[dim] /= total_dm_mass;

  double mag_kms2 = 0.0;
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    v_rel_kms_vec[dim] =
      (double(particle_velocity[dim][bh_index]) - com_vel_code[dim]) *
      velocity_units / 1.0e5;
    mag_kms2 += v_rel_kms_vec[dim] * v_rel_kms_vec[dim];
  }
  if (!isfinite(mag_kms2) || mag_kms2 <= 0.0)
    return;
  const double mag_kms = sqrt(mag_kms2);
  if (isfinite(candidate.v_rel_kms) && candidate.v_rel_kms > 0.0) {
    const double scale = candidate.v_rel_kms / mag_kms;
    for (int dim = 0; dim < MAX_DIMENSION; dim++)
      v_rel_kms_vec[dim] *= scale;
  }
}

static int BHDFMethod3SelectReason(
  int invalid_units, int invalid_bh_mass, int non_authoritative_mass,
  int newly_seeded, int child_covered, int under_resolved, int no_slow,
  int kernel_incomplete, int invalid_candidate_acceleration,
  int ln_lambda_invalid, int ln_lambda_zero, int low_vrel, int low_sigma,
  int dt_df_extreme, int cap_fraction_extreme)
{
  if (invalid_units)
    return 7;
  if (invalid_bh_mass)
    return 8;
  if (non_authoritative_mass)
    return 9;
  if (newly_seeded)
    return 10;
  if (child_covered)
    return 12;
  if (under_resolved)
    return 2;
  if (no_slow)
    return 3;
  if (kernel_incomplete)
    return 15;
  if (invalid_candidate_acceleration)
    return 14;
  if (ln_lambda_invalid)
    return 16;
  if (ln_lambda_zero)
    return 5;
  if (low_vrel)
    return 17;
  if (low_sigma)
    return 18;
  if (dt_df_extreme)
    return 19;
  if (cap_fraction_extreme)
    return 23;
  return 1;
}

static void BHDFMethod3EmitDetail(
  FILE *logptr, int processor, int invocation_seq, int grid_id, int step,
  int level, double zred, long long bh_id, int particle_index,
  int detail_row_index, double dx_code, double mass_units,
  double length_units, double time_units, double velocity_units,
  double kernel_volume_cgs, double bh_mass_code, double bh_mass_msun,
  double mass_slow_code, double kernel_radius_code,
  double kernel_radius_phys_kpc, double kernel_radius_over_dx,
  const BHDFCandidateResult &candidate, int under_resolved, int no_slow,
  int reason_code, int neighbor_particles_total, int n_star_excluded,
  int n_bh_excluded, const double v_rel_kms_vec[],
  double dt_level_myr, double active_dt_code, double limited_scale,
  double limited_delta_v_code, int cap_triggered,
  double cap_fraction_applied)
{
  const char *reason_name = BHDFMethod3ReasonName(reason_code);
  fprintf(logptr,
          "[BHDF_ACTIVE3] schema_version=3 active_method=3 mode=shadow "
          "active_policy=formula0_shadow_design0 active_apply_param=0 "
          "shadow_mode=1 processor=%d invocation_seq=%d grid_id=%d "
          "step=%d level=%d z=%.17g bh_id=%lld particle_index=%d "
          "detail_row_index=%d dx_code=%.17g mass_units=%.17g "
          "length_units=%.17g time_units=%.17g velocity_units=%.17g "
          "kernel_volume_cgs=%.17g bh_mass_code=%.17g "
          "bh_mass_msun=%.17g mass_slow_code=%.17g "
          "kernel_radius_code=%.17g kernel_radius_phys_kpc=%.17g "
          "kernel_radius_over_dx=%.17g kernel_complete=%d "
          "kernel_truncated_by_grid=%d n_dm_in_kernel=%d n_slow=%d "
          "under_resolved=%d no_slow=%d rejection_reason_code=%d "
          "rejection_reason_name=%s neighbor_particles_total=%d "
          "n_star_excluded=%d n_bh_excluded=%d rho_slow_cgs=%.17g "
          "v_rel_x_kms=%.17g v_rel_y_kms=%.17g v_rel_z_kms=%.17g "
          "v_rel_kms=%.17g sigma_1d_kms=%.17g b_min_phys_kpc=%.17g "
          "b_max_phys_kpc=%.17g ln_lambda=%.17g "
          "candidate_a_DF_x_cgs=%.17g candidate_a_DF_y_cgs=%.17g "
          "candidate_a_DF_z_cgs=%.17g candidate_a_DF_mag_cgs=%.17g "
          "dt_DF_Myr=%.17g dt_level_Myr=%.17g "
          "dt_level_over_dt_DF=%.17g requested_delta_v_code=%.17g "
          "limited_a_DF_x_cgs=%.17g limited_a_DF_y_cgs=%.17g "
          "limited_a_DF_z_cgs=%.17g limited_a_DF_mag_cgs=%.17g "
          "limited_delta_v_code=%.17g cap_fraction_requested=%.17g "
          "cap_fraction_applied=%.17g cap_triggered=%d "
          "per_substep_cap=%.17g per_level_cap=%.17g subcycle_count=1 "
          "subcycle_dt_code=%.17g active_dt_code=%.17g max_subcycles=%d "
          "local_gravity_required=0 local_gravity_available=0 "
          "local_gravity_magnitude_cgs=-1 a_DF_over_a_grav=-1 "
          "local_gravity_guard_passed=0 applied_delta_v_code=0 "
          "applied_delta_v_x_code=0 applied_delta_v_y_code=0 "
          "applied_delta_v_z_code=0 active_applied=0 kick_applied=0 "
          "momentum_policy=none delta_p_bh_x_code=0 "
          "delta_p_bh_y_code=0 delta_p_bh_z_code=0 "
          "delta_p_wake_x_code=0 delta_p_wake_y_code=0 "
          "delta_p_wake_z_code=0 suppression_reason_code=%d "
          "suppression_reason_name=%s active_validation_flags=none "
          "primary_reason_precedence=%s\n",
          processor, invocation_seq, grid_id, step, level, zred, bh_id,
          particle_index, detail_row_index, dx_code, mass_units, length_units,
          time_units, velocity_units, kernel_volume_cgs, bh_mass_code,
          bh_mass_msun, mass_slow_code, kernel_radius_code,
          kernel_radius_phys_kpc, kernel_radius_over_dx,
          candidate.kernel_complete, candidate.kernel_truncated_by_grid,
          candidate.n_dm_in_kernel, candidate.n_slow, under_resolved,
          no_slow, reason_code, reason_name, neighbor_particles_total,
          n_star_excluded, n_bh_excluded,
          BHDFMethod3CleanReal(candidate.rho_slow_cgs),
          BHDFMethod3CleanReal(v_rel_kms_vec[0]),
          BHDFMethod3CleanReal(v_rel_kms_vec[1]),
          BHDFMethod3CleanReal(v_rel_kms_vec[2]),
          BHDFMethod3CleanReal(candidate.v_rel_kms),
          BHDFMethod3CleanReal(candidate.sigma_1d_kms),
          BHDFMethod3CleanReal(candidate.b_min_phys_kpc),
          BHDFMethod3CleanReal(candidate.b_max_phys_kpc),
          BHDFMethod3CleanReal(candidate.ln_lambda),
          BHDFMethod3CleanReal(candidate.candidate_a_DF_cgs[0]),
          BHDFMethod3CleanReal(candidate.candidate_a_DF_cgs[1]),
          BHDFMethod3CleanReal(candidate.candidate_a_DF_cgs[2]),
          BHDFMethod3CleanReal(candidate.candidate_a_DF_mag_cgs),
          BHDFMethod3CleanReal(candidate.dt_DF_Myr),
          BHDFMethod3CleanReal(dt_level_myr),
          BHDFMethod3CleanReal(candidate.dt_level_over_dt_DF),
          BHDFMethod3CleanReal(candidate.requested_delta_v_code),
          BHDFMethod3CleanReal(candidate.candidate_a_DF_cgs[0] *
                               limited_scale),
          BHDFMethod3CleanReal(candidate.candidate_a_DF_cgs[1] *
                               limited_scale),
          BHDFMethod3CleanReal(candidate.candidate_a_DF_cgs[2] *
                               limited_scale),
          BHDFMethod3CleanReal(candidate.candidate_a_DF_mag_cgs *
                               limited_scale),
          BHDFMethod3CleanReal(limited_delta_v_code),
          BHDFMethod3CleanReal(candidate.cap_fraction_requested),
          BHDFMethod3CleanReal(cap_fraction_applied),
          cap_triggered,
          double(BHDynamicalFrictionActiveMaxKickFraction),
          double(BHDynamicalFrictionActiveMaxLevelKickFraction),
          BHDFMethod3CleanReal(active_dt_code),
          BHDFMethod3CleanReal(active_dt_code),
          BHDynamicalFrictionActiveMaxSubcycles,
          reason_code, reason_name, BHDFMethod3PrimaryReasonPrecedence());
}

static void BHDFMethod3EmitSummary(
  FILE *logptr, int processor, int invocation_seq, int grid_id, int step,
  int level, const BHDFMethod3Summary &summary)
{
  const int shadow_no_apply_rows = summary.reason_counts[1];
  const int applied_kicks = 0;
  const int shadow_rows = shadow_no_apply_rows;
  const int suppressed_rows = summary.detail_rows_emitted - shadow_rows;
  const int low_vrel_rows = summary.reason_counts[17];
  const int low_sigma_rows = summary.reason_counts[18];
  const int local_gravity_rows =
    summary.reason_counts[20] + summary.reason_counts[21];
  const int summary_accounting_passed =
    (summary.detail_rows_emitted >= 0 &&
     suppressed_rows >= 0 &&
     shadow_rows + suppressed_rows == summary.detail_rows_emitted) ? 1 : 0;
  const int detail_rows_complete =
    (summary.detail_rows_emitted == summary.eligible_bh) ? 1 : 0;

  fprintf(logptr,
          "[BHDF_ACTIVE3_SUMMARY] schema_version=3 active_method=3 "
          "mode=shadow processor=%d invocation_seq=%d grid_id=%d step=%d "
          "level=%d eligible_bh=%d detail_rows_emitted=%d active_rows=%d "
          "applied_kicks=%d suppressed_rows=%d shadow_rows=%d "
          "shadow_no_apply_rows=%d cap_limited_rows=%d "
          "local_gravity_suppressed_rows=%d low_sigma_suppressed_rows=%d "
          "low_vrel_suppressed_rows=%d momentum_delta_bh_x_sum=0 "
          "momentum_delta_bh_y_sum=0 momentum_delta_bh_z_sum=0 "
          "momentum_delta_wake_x_sum=0 momentum_delta_wake_y_sum=0 "
          "momentum_delta_wake_z_sum=0 summary_accounting_passed=%d "
          "detail_rows_complete=%d under_resolved_rows=%d "
          "no_slow_particles_rows=%d kernel_truncated_rows=%d "
          "ln_lambda_zero_rows=%d guard_skip_rows=%d invalid_units_rows=%d "
          "invalid_bh_mass_rows=%d non_authoritative_mass_rows=%d "
          "newly_seeded_rows=%d outside_grid_rows=%d child_covered_rows=%d "
          "reposition_conflict_rows=%d invalid_candidate_acceleration_rows=%d "
          "kernel_incomplete_rows=%d ln_lambda_invalid_rows=%d "
          "low_vrel_rows=%d low_sigma_rows=%d dt_df_extreme_rows=%d "
          "local_gravity_unavailable_rows=%d "
          "local_gravity_ratio_exceeded_rows=%d "
          "momentum_policy_unsupported_rows=%d cap_fraction_extreme_rows=%d "
          "max_subcycles_exceeded_rows=%d apply_not_authorized_rows=%d\n",
          processor, invocation_seq, grid_id, step, level,
          summary.eligible_bh, summary.detail_rows_emitted,
          summary.detail_rows_emitted, applied_kicks, suppressed_rows,
          shadow_rows, shadow_no_apply_rows, summary.cap_limited_rows,
          local_gravity_rows, low_sigma_rows, low_vrel_rows,
          summary_accounting_passed, detail_rows_complete,
          summary.reason_counts[2], summary.reason_counts[3],
          summary.reason_counts[4], summary.reason_counts[5],
          summary.reason_counts[6], summary.reason_counts[7],
          summary.reason_counts[8], summary.reason_counts[9],
          summary.reason_counts[10], summary.reason_counts[11],
          summary.reason_counts[12], summary.reason_counts[13],
          summary.reason_counts[14], summary.reason_counts[15],
          summary.reason_counts[16], summary.reason_counts[17],
          summary.reason_counts[18], summary.reason_counts[19],
          summary.reason_counts[20], summary.reason_counts[21],
          summary.reason_counts[22], summary.reason_counts[23],
          summary.reason_counts[24], summary.reason_counts[25]);
}

static int BHDFMethod3ShadowSkeleton(
  FILE *logptr, HierarchyEntry *SubgridPointer, int processor,
  int invocation_seq, int grid_id, int step, int level, double zred,
  int number_of_particles, int *particle_type, PINT *particle_number,
  float *particle_mass, FLOAT *particle_position[],
  float *particle_velocity[], int number_of_particle_attributes,
  float *particle_attribute[], double mass_units, double length_units,
  double time_units, double velocity_units, double cell_volume_code,
  const FLOAT *grid_left_edge, const FLOAT *grid_right_edge,
  double dx_code, double kernel_radius_code,
  double kernel_radius_phys_kpc, double kernel_radius_over_dx,
  double kernel_volume_cgs, double dt_level_myr, double active_dt_code,
  int invalid_units)
{
  if (logptr == NULL)
    logptr = stdout;

  std::vector<int> bh_particles;
  if (number_of_particles > 0 && particle_type != NULL) {
    bh_particles.reserve(number_of_particles);
    for (int p = 0; p < number_of_particles; p++)
      if (BHDFIsBHType(particle_type[p]))
        bh_particles.push_back(p);
  }
  if (particle_number != NULL)
    std::sort(bh_particles.begin(), bh_particles.end(),
              BHDFParticleOrder(particle_number));

  BHDFMethod3Summary summary;
  BHDFMethod3InitSummary(&summary);
  summary.eligible_bh = int(bh_particles.size());

  const double mass_to_msun =
    (isfinite(mass_units) && mass_units > 0.0) ? mass_units / SolarMass : 0.0;

  for (int bp = 0; bp < int(bh_particles.size()); bp++) {
    const int p = bh_particles[bp];
    double bh_pos[MAX_DIMENSION] = {0.0, 0.0, 0.0};
    for (int dim = 0; dim < MAX_DIMENSION; dim++)
      if (particle_position != NULL && particle_position[dim] != NULL)
        bh_pos[dim] = double(particle_position[dim][p]);

    FLOAT bh_pos_float[MAX_DIMENSION] = {
      FLOAT(bh_pos[0]), FLOAT(bh_pos[1]), FLOAT(bh_pos[2])
    };

    double bh_mass_code = 0.0;
    if (particle_mass != NULL && isfinite(double(particle_mass[p])) &&
        double(particle_mass[p]) > 0.0)
      bh_mass_code = double(particle_mass[p]);

    double attribute_mass = -1.0;
    if (BHDFMassFromAttributes(number_of_particle_attributes,
                               particle_attribute, p, &attribute_mass))
      bh_mass_code = attribute_mass;

    const int invalid_bh_mass =
      (!isfinite(bh_mass_code) || bh_mass_code <= 0.0) ? 1 : 0;
    double expected_full_mass = -1.0;
    double mass_attribute_ratio = -1.0;
    const int non_authoritative_mass =
      (particle_mass != NULL &&
       BHDFNonAuthoritativeMassRatio(
         number_of_particle_attributes, particle_attribute, p,
         double(particle_mass[p]), &expected_full_mass,
         &mass_attribute_ratio)) ? 1 : 0;
    const int newly_seeded =
      (BHDynamicalFrictionActiveSuppressNewSeeds != 0 &&
       particle_number != NULL &&
       BHDFIsNewlySeededThisPass(particle_number[p])) ? 1 : 0;
    const int point_in_grid =
      (grid_left_edge != NULL && grid_right_edge != NULL) ?
      ((bh_pos[0] >= double(grid_left_edge[0]) &&
        bh_pos[0] <= double(grid_right_edge[0]) &&
        bh_pos[1] >= double(grid_left_edge[1]) &&
        bh_pos[1] <= double(grid_right_edge[1]) &&
        bh_pos[2] >= double(grid_left_edge[2]) &&
        bh_pos[2] <= double(grid_right_edge[2])) ? 1 : 0) : 1;
    const int child_covered =
      (SubgridPointer != NULL &&
       BHDFPositionCoveredByChild(SubgridPointer, bh_pos_float)) ? 1 : 0;
    const int mass_or_grid_not_authoritative =
      (non_authoritative_mass || !point_in_grid) ? 1 : 0;

    BHDFCandidateResult candidate;
    BHDFInitCandidateResult(&candidate, kernel_radius_phys_kpc);
    if (!invalid_units && !invalid_bh_mass &&
        !mass_or_grid_not_authoritative && !newly_seeded && !child_covered) {
      candidate = BHDFComputeCandidate(
        number_of_particles, particle_type, particle_mass, particle_position,
        particle_velocity, p, bh_pos, bh_mass_code, kernel_radius_code,
        kernel_radius_phys_kpc, dt_level_myr, mass_units, length_units,
        velocity_units, cell_volume_code, grid_left_edge, grid_right_edge);
    }

    int n_star_excluded = 0;
    int n_bh_excluded = 0;
    BHDFMethod3NeighborPartition(
      number_of_particles, particle_type, particle_position, bh_pos,
      kernel_radius_code, &n_star_excluded, &n_bh_excluded);
    const int neighbor_particles_total =
      candidate.n_dm_in_kernel + n_star_excluded + n_bh_excluded;

    double v_rel_kms_vec[MAX_DIMENSION] = {0.0, 0.0, 0.0};
    BHDFMethod3VRelVectorKms(
      number_of_particles, particle_type, particle_mass, particle_position,
      particle_velocity, p, bh_pos, kernel_radius_code, velocity_units,
      candidate, v_rel_kms_vec);

    const int under_resolved =
      (candidate.n_dm_in_kernel <
       BHDynamicalFrictionActiveMinParticles) ? 1 : 0;
    const int no_slow =
      (candidate.n_slow <
       BHDynamicalFrictionActiveMinSlowParticles) ? 1 : 0;
    const int kernel_incomplete =
      (BHDynamicalFrictionActiveRequireKernelComplete != 0 &&
       candidate.kernel_complete == 0) ? 1 : 0;
    const int invalid_candidate_acceleration =
      (!isfinite(candidate.candidate_a_DF_mag_cgs) ||
       candidate.candidate_a_DF_mag_cgs < 0.0 ||
       !isfinite(candidate.candidate_a_DF_cgs[0]) ||
       !isfinite(candidate.candidate_a_DF_cgs[1]) ||
       !isfinite(candidate.candidate_a_DF_cgs[2])) ? 1 : 0;
    const int ln_lambda_invalid =
      (!isfinite(candidate.ln_lambda) || candidate.ln_lambda < 0.0) ? 1 : 0;
    const int ln_lambda_zero =
      (!ln_lambda_invalid && candidate.ln_lambda <= 0.0) ? 1 : 0;
    const int low_vrel =
      (candidate.computed && BHDynamicalFrictionActiveVrelFloorKmS > 0.0f &&
       candidate.v_rel_kms < double(BHDynamicalFrictionActiveVrelFloorKmS)) ?
      1 : 0;
    const int low_sigma =
      (candidate.computed && BHDynamicalFrictionActiveSigmaFloorKmS > 0.0f &&
       candidate.sigma_1d_kms <
       double(BHDynamicalFrictionActiveSigmaFloorKmS)) ? 1 : 0;
    const int dt_df_extreme =
      (candidate.computed &&
       candidate.dt_level_over_dt_DF >
       double(BHDynamicalFrictionActiveMaxDtOverDtDF)) ? 1 : 0;

    double limited_delta_v_code = candidate.requested_delta_v_code;
    int cap_triggered = 0;
    double limited_scale = 1.0;
    if (!isfinite(limited_delta_v_code) || limited_delta_v_code < 0.0)
      limited_delta_v_code = 0.0;
    const double velocity_kms_to_code =
      (isfinite(velocity_units) && velocity_units > 0.0) ?
      1.0e5 / velocity_units : 0.0;
    const double v_rel_code = candidate.v_rel_kms * velocity_kms_to_code;
    if (isfinite(v_rel_code) && v_rel_code > 0.0) {
      const double substep_limit =
        double(BHDynamicalFrictionActiveMaxKickFraction) * v_rel_code;
      const double level_limit =
        double(BHDynamicalFrictionActiveMaxLevelKickFraction) * v_rel_code;
      const double cap_limit = min(substep_limit, level_limit);
      if (isfinite(cap_limit) && cap_limit >= 0.0 &&
          limited_delta_v_code > cap_limit) {
        limited_delta_v_code = cap_limit;
        cap_triggered = 1;
      }
    }
    if (candidate.requested_delta_v_code > 0.0 &&
        isfinite(candidate.requested_delta_v_code))
      limited_scale = limited_delta_v_code / candidate.requested_delta_v_code;
    else
      limited_scale = 0.0;
    const int cap_fraction_extreme =
      (!isfinite(candidate.cap_fraction_requested) ||
       candidate.cap_fraction_requested < 0.0 ||
       !isfinite(limited_delta_v_code) || limited_delta_v_code < 0.0) ? 1 : 0;

    const int reason_code = BHDFMethod3SelectReason(
      invalid_units, invalid_bh_mass, mass_or_grid_not_authoritative,
      newly_seeded, child_covered, under_resolved, no_slow,
      kernel_incomplete, invalid_candidate_acceleration, ln_lambda_invalid,
      ln_lambda_zero, low_vrel, low_sigma, dt_df_extreme,
      cap_fraction_extreme);

    double mass_slow_code = 0.0;
    if (isfinite(candidate.rho_slow_cgs) && candidate.rho_slow_cgs > 0.0 &&
        isfinite(kernel_volume_cgs) && kernel_volume_cgs > 0.0 &&
        isfinite(mass_units) && mass_units > 0.0)
      mass_slow_code =
        candidate.rho_slow_cgs * kernel_volume_cgs / mass_units;

    BHDFMethod3EmitDetail(
      logptr, processor, invocation_seq, grid_id, step, level, zred,
      (particle_number != NULL) ? (long long) particle_number[p] : -1, p,
      summary.detail_rows_emitted, dx_code, mass_units, length_units,
      time_units, velocity_units, kernel_volume_cgs, bh_mass_code,
      bh_mass_code * mass_to_msun, mass_slow_code, kernel_radius_code,
      kernel_radius_phys_kpc, kernel_radius_over_dx, candidate,
      under_resolved, no_slow, reason_code, neighbor_particles_total,
      n_star_excluded, n_bh_excluded, v_rel_kms_vec, dt_level_myr,
      active_dt_code, limited_scale, limited_delta_v_code, cap_triggered,
      0.0);

    summary.detail_rows_emitted++;
    if (reason_code >= 0 && reason_code < 26)
      summary.reason_counts[reason_code]++;
    if (cap_triggered)
      summary.cap_limited_rows++;
  }

  BHDFMethod3EmitSummary(
    logptr, processor, invocation_seq, grid_id, step, level, summary);
  return SUCCESS;
}

int grid::BHDynamicalFrictionHandler(HierarchyEntry* SubgridPointer,
                                     int level,
                                     int cycle_number,
                                     float dtLevelAbove)
{
  (void) dtLevelAbove;
  const double dt_code = max(0.0, double(this->dtFixed));

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (BHDynamicalFrictionMethod == 3) {
    if (GridRank != 3)
      return SUCCESS;

    FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;
    static int bhdf_active3_invocation_seq_counter = 0;
    const int bhdf_active3_invocation_seq =
      ++bhdf_active3_invocation_seq_counter;

    FLOAT a = 1.0, dadt = 0.0;
    if (ComovingCoordinates)
      CosmologyComputeExpansionFactor(Time, &a, &dadt);
    const double zred = ComovingCoordinates ?
      double((1.0 + InitialRedshift) / a - 1.0) : 0.0;

    float DensityUnits = 1.0f, LengthUnits = 1.0f;
    float TemperatureUnits = 1.0f, TimeUnits = 1.0f, VelocityUnits = 1.0f;
    const int units_valid =
      (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                &TimeUnits, &VelocityUnits, Time) != FAIL);

    const int cell_widths_available =
      (CellWidth[0] != NULL && CellWidth[1] != NULL &&
       CellWidth[2] != NULL) ? 1 : 0;
    const double dx_code =
      cell_widths_available ? double(CellWidth[0][0]) : 0.0;
    const double cell_volume_code =
      cell_widths_available ?
      double(CellWidth[0][0]) * double(CellWidth[1][0]) *
      double(CellWidth[2][0]) : 0.0;
    const double mass_units =
      double(DensityUnits) * double(LengthUnits) *
      double(LengthUnits) * double(LengthUnits);
    const double kernel_radius_code =
      BHDFKernelRadiusCode(BHDynamicalFrictionKernelRadius, Time,
                           LengthUnits);
    const double kernel_radius_phys_kpc =
      isfinite(double(BHDynamicalFrictionKernelRadius)) ?
      double(BHDynamicalFrictionKernelRadius) : 0.0;
    const double kernel_radius_over_dx =
      (dx_code > 0.0) ? kernel_radius_code / dx_code : 0.0;
    const double radius_cgs = kernel_radius_code * double(LengthUnits);
    const double kernel_volume_cgs =
      (4.0 * M_PI / 3.0) * radius_cgs * radius_cgs * radius_cgs;
    const double dt_level_myr =
      (TimeUnits > 0.0f) ? dt_code * double(TimeUnits) / Myr_s : 0.0;
    const int invalid_units =
      (!units_valid || !isfinite(mass_units) || mass_units <= 0.0 ||
       !isfinite(double(LengthUnits)) || double(LengthUnits) <= 0.0 ||
       !isfinite(double(TimeUnits)) || double(TimeUnits) <= 0.0 ||
       !isfinite(double(VelocityUnits)) || double(VelocityUnits) <= 0.0 ||
       !isfinite(cell_volume_code) || cell_volume_code <= 0.0 ||
       !isfinite(dx_code) || dx_code <= 0.0 ||
       !isfinite(kernel_radius_code) || kernel_radius_code <= 0.0 ||
       !isfinite(kernel_volume_cgs) || kernel_volume_cgs <= 0.0) ? 1 : 0;

    return BHDFMethod3ShadowSkeleton(
      logptr, SubgridPointer, MyProcessorNumber, bhdf_active3_invocation_seq,
      this->GetGridID(), cycle_number, level, zred,
      NumberOfParticles, ParticleType, ParticleNumber, ParticleMass,
      ParticlePosition, ParticleVelocity,
      NumberOfParticleAttributes, ParticleAttribute,
      invalid_units ? 1.0 : mass_units,
      invalid_units ? 1.0 : double(LengthUnits),
      invalid_units ? 1.0 : double(TimeUnits),
      invalid_units ? 1.0 : double(VelocityUnits),
      invalid_units ? 1.0 : cell_volume_code,
      GridLeftEdge, GridRightEdge,
      invalid_units ? 1.0 : dx_code,
      invalid_units ? 1.0 : kernel_radius_code,
      invalid_units ? 1.0 : kernel_radius_phys_kpc,
      invalid_units ? 1.0 : kernel_radius_over_dx,
      invalid_units ? 1.0 : kernel_volume_cgs,
      invalid_units ? 0.0 : dt_level_myr,
      dt_code, invalid_units);
  }

  if (BHDynamicalFrictionMethod == 2) {
    if (BHDynamicalFrictionActiveVerbose <= 0)
      return SUCCESS;

    if (GridRank != 3)
      return SUCCESS;

    const double cell_volume_code =
      double(CellWidth[0][0]) * double(CellWidth[1][0]) *
      double(CellWidth[2][0]);

    FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;
    static int bhdf_active_invocation_seq_counter = 0;
    const int bhdf_active_invocation_seq =
      ++bhdf_active_invocation_seq_counter;

    float DensityUnits = 1.0f, LengthUnits = 1.0f, TemperatureUnits = 1.0f;
    float TimeUnits = 1.0f, VelocityUnits = 1.0f;
    const int units_valid =
      (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                &TimeUnits, &VelocityUnits, Time) != FAIL);

    double mass_units = 0.0;
    double kernel_radius_code = 0.0;
    double kernel_radius_over_dx = 0.0;
    double dt_level_myr = 0.0;
    if (units_valid) {
      mass_units = double(DensityUnits) *
        double(LengthUnits) * double(LengthUnits) * double(LengthUnits);
      kernel_radius_code =
        BHDFKernelRadiusCode(BHDynamicalFrictionKernelRadius, Time,
                             LengthUnits);
      const double cell_width = double(CellWidth[0][0]);
      if (cell_width > 0.0)
        kernel_radius_over_dx = kernel_radius_code / cell_width;
      if (TimeUnits > 0.0f)
        dt_level_myr = dt_code * double(TimeUnits) / Myr_s;
    }

    const double kernel_radius_phys_kpc =
      isfinite(double(BHDynamicalFrictionKernelRadius)) ?
      double(BHDynamicalFrictionKernelRadius) : 0.0;

    return BHDFActiveDryRunSkeleton(
      logptr, MyProcessorNumber, bhdf_active_invocation_seq,
      this->GetGridID(), cycle_number, level,
      NumberOfParticles, ParticleType, ParticleNumber, ParticleMass,
      ParticlePosition, ParticleVelocity,
      NumberOfParticleAttributes, ParticleAttribute,
      mass_units, double(LengthUnits), double(VelocityUnits),
      cell_volume_code,
      GridLeftEdge, GridRightEdge, kernel_radius_code,
      kernel_radius_phys_kpc, kernel_radius_over_dx, dt_level_myr);
  }

  if (NumberOfParticles <= 0)
    return SUCCESS;

  FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;

  static int bhdf_invocation_seq_counter = 0;
  const int bhdf_invocation_seq = ++bhdf_invocation_seq_counter;
  const int bhdf_processor = MyProcessorNumber;
  const int bhdf_grid_id = this->GetGridID();
  const int active_reposition_enabled = (BHRepositionMethod > 0) ? 1 : 0;

  FLOAT a = 1.0, dadt = 0.0;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  const double zred = ComovingCoordinates ?
    double((1.0 + InitialRedshift) / a - 1.0) : 0.0;

  if (GridRank != 3) {
    if (BHDynamicalFrictionVerbose >= 1) {
      fprintf(logptr,
              "[BHDF_WARN] schema_version=1 processor=%d invocation_seq=%d "
              "grid_id=%d step=%d level=%d warning=unsupported_grid_rank "
              "unsupported_grid_rank=1 grid_rank=%d\n",
              bhdf_processor, bhdf_invocation_seq, bhdf_grid_id,
              cycle_number, level, GridRank);
      fprintf(logptr,
              "[BHDF_SUMMARY] schema_version=1 processor=%d invocation_seq=%d "
              "grid_id=%d step=%d level=%d grid_rank=%d "
              "unsupported_grid_rank=1 active_reposition_enabled=%d "
              "rows=0 eligible_bh=0 valid_or_conditioned_rows=0 "
              "diagnostics_only_rows=0 active_rows=0 applied_kicks=0 "
              "skipped_child_covered=0 skipped_newly_seeded=0 "
              "skipped_non_authoritative_mass=0 no_local_neighbors_count=0 "
              "under_resolved_count=0 kernel_truncated_by_grid_count=0 "
              "no_slow_particles_count=0 ln_lambda_zero_count=0 "
              "ln_lambda_large_count=0 max_a_DF_cgs=0 min_dt_DF_Myr=-1.0 "
              "total_df_wall_ms=0 max_bh_df_wall_ms=0 "
              "detail_rows_emitted=0 detail_rows_complete=1\n",
              bhdf_processor, bhdf_invocation_seq, bhdf_grid_id,
              cycle_number, level, GridRank, active_reposition_enabled);
    }
    return SUCCESS;
  }

  float DensityUnits = 1.0f, LengthUnits = 1.0f, TemperatureUnits = 1.0f;
  float TimeUnits = 1.0f, VelocityUnits = 1.0f;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL)
    ENZO_FAIL("Error in GetUnits.");

  const float cell_width = float(CellWidth[0][0]);
  if (cell_width <= 0.0f)
    return SUCCESS;
  const double cell_volume_code =
    double(CellWidth[0][0]) * double(CellWidth[1][0]) *
    double(CellWidth[2][0]);
  if (!isfinite(cell_volume_code) || cell_volume_code <= 0.0)
    return SUCCESS;

  const double kernel_radius_code =
    BHDFKernelRadiusCode(BHDynamicalFrictionKernelRadius, Time, LengthUnits);
  if (kernel_radius_code <= 0.0)
    return SUCCESS;

  const double kernel_radius_phys_kpc = BHDynamicalFrictionKernelRadius;
  const double kernel_radius_over_dx = kernel_radius_code / double(cell_width);
  const double kernel_r2 = kernel_radius_code * kernel_radius_code;
  const double radius_cgs = kernel_radius_code * double(LengthUnits);
  const double kernel_volume_cgs =
    (4.0 * M_PI / 3.0) * radius_cgs * radius_cgs * radius_cgs;
  const double mass_units = double(DensityUnits) *
    double(LengthUnits) * double(LengthUnits) * double(LengthUnits);
  const double mass_to_msun = mass_units / SolarMass;
  const double b_max_phys_kpc = BHDynamicalFrictionKernelRadius;
  const double b_max_phys_cm = b_max_phys_kpc * kpc_cm;

  std::vector<int> bh_particles;
  bh_particles.reserve(NumberOfParticles);
  for (int p = 0; p < NumberOfParticles; p++)
    if (BHDFIsBHType(ParticleType[p]))
      bh_particles.push_back(p);

  if (bh_particles.empty())
    return SUCCESS;

  std::sort(bh_particles.begin(), bh_particles.end(),
            BHDFParticleOrder(ParticleNumber));

  if (BHDynamicalFrictionVerbose >= 1 && active_reposition_enabled) {
    fprintf(logptr,
            "[BHDF_WARN] schema_version=1 processor=%d invocation_seq=%d "
            "grid_id=%d step=%d level=%d warning=active_reposition_enabled "
            "active_reposition_enabled=1 BHRepositionMethod=%d\n",
            bhdf_processor, bhdf_invocation_seq, bhdf_grid_id,
            cycle_number, level, BHRepositionMethod);
  }

  int rows = 0;
  const int eligible_bh = int(bh_particles.size());
  int valid_or_conditioned_rows = 0;
  int diagnostics_only_rows = 0;
  int skipped_child_covered = 0;
  int skipped_newly_seeded = 0;
  int skipped_non_authoritative_mass = 0;
  int no_local_neighbors_count = 0;
  int under_resolved_count = 0;
  int kernel_truncated_by_grid_count = 0;
  int no_slow_particles_count = 0;
  int ln_lambda_zero_count = 0;
  int ln_lambda_large_count = 0;
  int detail_rows_emitted = 0;
  double max_a_df_cgs = 0.0;
  double min_dt_df_myr = -1.0;
  double total_df_wall_ms = 0.0;
  double max_bh_df_wall_ms = 0.0;
  const int detail_rows_complete_flag =
    (BHDynamicalFrictionVerbose >= 2) ? 1 : 0;

  for (int bp = 0; bp < int(bh_particles.size()); bp++) {
    const int p = bh_particles[bp];
    const double t0 = ReturnWallTime();

    BHDFRow row;
    BHDFInitRow(&row);
    row.processor = bhdf_processor;
    row.invocation_seq = bhdf_invocation_seq;
    row.grid_id = bhdf_grid_id;
    row.step = cycle_number;
    row.level = level;
    row.grid_rank = GridRank;
    row.z = zred;
    row.bh_id = (long long) ParticleNumber[p];
    row.particle_index = p;
    row.active_reposition_enabled = active_reposition_enabled;
    row.kernel_radius_code = kernel_radius_code;
    row.kernel_radius_phys_kpc = kernel_radius_phys_kpc;
    row.kernel_radius_over_dx = kernel_radius_over_dx;
    row.b_max_phys_kpc = b_max_phys_kpc;

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      row.bh_pos[dim] = double(ParticlePosition[dim][p]);
      row.bh_vel_kms[dim] =
        double(ParticleVelocity[dim][p]) * double(VelocityUnits) / 1.0e5;
    }
    FLOAT bh_pos_float[MAX_DIMENSION] = {
      FLOAT(row.bh_pos[0]), FLOAT(row.bh_pos[1]), FLOAT(row.bh_pos[2])
    };

    const double particle_mass_code = double(ParticleMass[p]);
    double expected_full_mass = -1.0;
    double mass_attribute_ratio = -1.0;
    double bh_mass_code = particle_mass_code;
    if (BHDFMassFromAttributes(NumberOfParticleAttributes, ParticleAttribute,
                               p, &expected_full_mass))
      bh_mass_code = expected_full_mass;
    row.bh_mass_msun = bh_mass_code * mass_to_msun;

    rows++;

    if (!isfinite(particle_mass_code) || particle_mass_code <= 0.0 ||
        !isfinite(bh_mass_code) || bh_mass_code <= 0.0 ||
        !this->PointInGrid(bh_pos_float)) {
      row.authoritative = 0;
      row.rejection_reason = 7;
      skipped_non_authoritative_mass++;
      row.df_wall_ms = 1000.0 * (ReturnWallTime() - t0);
      total_df_wall_ms += row.df_wall_ms;
      max_bh_df_wall_ms = max(max_bh_df_wall_ms, row.df_wall_ms);
      if (BHDynamicalFrictionVerbose >= 1)
        BHDFEmitGuardWarning(logptr, row);
      if (BHDynamicalFrictionVerbose >= 2) {
        BHDFEmitDetail(logptr, row, detail_rows_emitted + 1,
                       detail_rows_complete_flag);
        detail_rows_emitted++;
      }
      continue;
    }

    if (SubgridPointer != NULL &&
        BHDFPositionCoveredByChild(SubgridPointer, bh_pos_float)) {
      row.authoritative = 0;
      row.rejection_reason = 8;
      skipped_child_covered++;
      row.df_wall_ms = 1000.0 * (ReturnWallTime() - t0);
      total_df_wall_ms += row.df_wall_ms;
      max_bh_df_wall_ms = max(max_bh_df_wall_ms, row.df_wall_ms);
      if (BHDynamicalFrictionVerbose >= 1)
        BHDFEmitGuardWarning(logptr, row);
      if (BHDynamicalFrictionVerbose >= 2) {
        BHDFEmitDetail(logptr, row, detail_rows_emitted + 1,
                       detail_rows_complete_flag);
        detail_rows_emitted++;
      }
      continue;
    }

    if (BHDFIsNewlySeededThisPass(ParticleNumber[p])) {
      row.rejection_reason = 6;
      skipped_newly_seeded++;
      row.df_wall_ms = 1000.0 * (ReturnWallTime() - t0);
      total_df_wall_ms += row.df_wall_ms;
      max_bh_df_wall_ms = max(max_bh_df_wall_ms, row.df_wall_ms);
      if (BHDynamicalFrictionVerbose >= 1)
        BHDFEmitGuardWarning(logptr, row);
      if (BHDynamicalFrictionVerbose >= 2) {
        BHDFEmitDetail(logptr, row, detail_rows_emitted + 1,
                       detail_rows_complete_flag);
        detail_rows_emitted++;
      }
      continue;
    }

    if (BHDFNonAuthoritativeMassRatio(
          NumberOfParticleAttributes, ParticleAttribute, p,
          particle_mass_code, &expected_full_mass, &mass_attribute_ratio)) {
      row.authoritative = 0;
      row.rejection_reason = 7;
      skipped_non_authoritative_mass++;
      row.df_wall_ms = 1000.0 * (ReturnWallTime() - t0);
      total_df_wall_ms += row.df_wall_ms;
      max_bh_df_wall_ms = max(max_bh_df_wall_ms, row.df_wall_ms);
      if (BHDynamicalFrictionVerbose >= 1)
        BHDFEmitGuardWarning(logptr, row);
      if (BHDynamicalFrictionVerbose >= 2) {
        BHDFEmitDetail(logptr, row, detail_rows_emitted + 1,
                       detail_rows_complete_flag);
        detail_rows_emitted++;
      }
      continue;
    }

    valid_or_conditioned_rows++;
    diagnostics_only_rows++;

    row.kernel_complete = 1;
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      if (row.bh_pos[dim] - kernel_radius_code < double(GridLeftEdge[dim]) ||
          row.bh_pos[dim] + kernel_radius_code > double(GridRightEdge[dim]))
        row.kernel_complete = 0;
    }
    row.kernel_truncated_by_grid = row.kernel_complete ? 0 : 1;

    double total_dm_mass = 0.0;
    double com_vel_code[MAX_DIMENSION] = {0.0, 0.0, 0.0};
    std::vector<int> dm_neighbors;
    dm_neighbors.reserve(NumberOfParticles);

    for (int q = 0; q < NumberOfParticles; q++) {
      double dr2 = 0.0;
      for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        const double dr =
          double(ParticlePosition[dim][q]) - row.bh_pos[dim];
        dr2 += dr * dr;
      }
      if (dr2 > kernel_r2)
        continue;

      if (BHDFIsDarkMatter(ParticleType[q]) || BHDFIsStar(ParticleType[q]) ||
          BHDFIsBHType(ParticleType[q]))
        row.neighbor_particles_total++;

      if (BHDFIsStar(ParticleType[q])) {
        row.n_star_excluded++;
        continue;
      }
      if (BHDFIsBHType(ParticleType[q])) {
        row.n_bh_excluded++;
        continue;
      }
      if (!BHDFIsDarkMatter(ParticleType[q]))
        continue;

      const double dm_mass = double(ParticleMass[q]);
      if (!isfinite(dm_mass) || dm_mass <= 0.0)
        continue;

      row.n_dm_in_kernel++;
      total_dm_mass += dm_mass;
      dm_neighbors.push_back(q);
      for (int dim = 0; dim < MAX_DIMENSION; dim++)
        com_vel_code[dim] += dm_mass * double(ParticleVelocity[dim][q]);
    }

    if (total_dm_mass > 0.0) {
      for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        com_vel_code[dim] /= total_dm_mass;
        row.v_com_kms[dim] = com_vel_code[dim] *
          double(VelocityUnits) / 1.0e5;
      }
    } else {
      no_local_neighbors_count++;
      for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        com_vel_code[dim] = double(ParticleVelocity[dim][p]);
        row.v_com_kms[dim] = row.bh_vel_kms[dim];
      }
    }

    double v_rel_code_vec[MAX_DIMENSION] = {0.0, 0.0, 0.0};
    double v_rel_code2 = 0.0;
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      v_rel_code_vec[dim] =
        double(ParticleVelocity[dim][p]) - com_vel_code[dim];
      v_rel_code2 += v_rel_code_vec[dim] * v_rel_code_vec[dim];
    }
    const double v_rel_code = sqrt(v_rel_code2);
    const double v_rel_cgs = v_rel_code * double(VelocityUnits);
    row.v_bh_rel_kms = v_rel_cgs / 1.0e5;

    double sigma2_num = 0.0;
    for (int ni = 0; ni < int(dm_neighbors.size()); ni++) {
      const int q = dm_neighbors[ni];
      const double dm_mass = double(ParticleMass[q]);
      double dv2 = 0.0;
      for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        const double dv =
          double(ParticleVelocity[dim][q]) - com_vel_code[dim];
        dv2 += dv * dv;
      }
      sigma2_num += dm_mass * dv2;
      if (sqrt(dv2) < v_rel_code) {
        row.n_slow++;
        row.mass_slow_code += dm_mass;
      }
    }

    double sigma_1d_code = 0.0;
    if (total_dm_mass > 0.0)
      sigma_1d_code = sqrt(sigma2_num / (3.0 * total_dm_mass));
    row.sigma_1d_kms = sigma_1d_code * double(VelocityUnits) / 1.0e5;

    // MASSAUDIT0/MASSFIX0: native Enzo DM ParticleMass is density-like.
    // Convert the accumulated stored slow-particle density to true code mass
    // at the endpoint before density/acceleration/cap calculations. This
    // preserves the mass-weighted velocity chain; BH masses are absolute.
    row.mass_slow_code *= cell_volume_code;
    if (kernel_volume_cgs > 0.0)
      row.rho_slow_cgs =
        row.mass_slow_code * mass_units / kernel_volume_cgs;

    const double bh_mass_cgs = bh_mass_code * mass_units;
    const double sigma_1d_cgs = sigma_1d_code * double(VelocityUnits);
    const double bmin_denom =
      v_rel_cgs * v_rel_cgs + sigma_1d_cgs * sigma_1d_cgs;
    if (bmin_denom > 0.0 && isfinite(bh_mass_cgs) && bh_mass_cgs > 0.0) {
      const double b_min_phys_cm = GravConst * bh_mass_cgs / bmin_denom;
      row.b_min_phys_kpc = b_min_phys_cm / kpc_cm;
      if (b_min_phys_cm > 0.0 && b_min_phys_cm < b_max_phys_cm)
        row.ln_lambda = log(b_max_phys_cm / b_min_phys_cm);
    }
    if (!isfinite(row.ln_lambda) || row.ln_lambda < 0.0)
      row.ln_lambda = 0.0;
    row.ln_lambda_zero = (row.ln_lambda <= 0.0) ? 1 : 0;
    row.ln_lambda_large = (row.ln_lambda > 20.0) ? 1 : 0;

    if (v_rel_cgs > 0.0 && row.rho_slow_cgs > 0.0 &&
        row.ln_lambda > 0.0) {
      const double prefactor =
        -4.0 * M_PI * GravConst * GravConst * bh_mass_cgs *
        row.rho_slow_cgs * row.ln_lambda /
        (v_rel_cgs * v_rel_cgs * v_rel_cgs);
      for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        const double v_rel_dim_cgs =
          v_rel_code_vec[dim] * double(VelocityUnits);
        row.a_df_cgs[dim] = prefactor * v_rel_dim_cgs;
        row.a_df_mag_cgs += row.a_df_cgs[dim] * row.a_df_cgs[dim];
      }
      row.a_df_mag_cgs = sqrt(row.a_df_mag_cgs);
    }

    if (row.a_df_mag_cgs > 0.0 && v_rel_cgs > 0.0)
      row.dt_df_myr = (v_rel_cgs / row.a_df_mag_cgs) / Myr_s;

    const double hypothetical_delta_v_code =
      (row.a_df_mag_cgs / double(VelocityUnits)) *
      double(TimeUnits) * dt_code;
    if (v_rel_code > 0.0)
      row.velocity_kick_cap_fraction =
        hypothetical_delta_v_code / v_rel_code;
    row.cap_triggered = (row.velocity_kick_cap_fraction > 0.5) ? 1 : 0;

    row.under_resolved =
      (row.n_dm_in_kernel < BHDynamicalFrictionMinParticles) ? 1 : 0;
    const int no_slow =
      (row.n_slow < BHDynamicalFrictionMinSlowParticles) ? 1 : 0;

    if (row.under_resolved)
      row.rejection_reason = 5;
    else if (no_slow)
      row.rejection_reason = 2;
    else if (row.ln_lambda_zero)
      row.rejection_reason = 3;
    else
      row.rejection_reason = 0;

    under_resolved_count += row.under_resolved ? 1 : 0;
    no_slow_particles_count += no_slow ? 1 : 0;
    kernel_truncated_by_grid_count += row.kernel_truncated_by_grid ? 1 : 0;
    ln_lambda_zero_count += row.ln_lambda_zero ? 1 : 0;
    ln_lambda_large_count += row.ln_lambda_large ? 1 : 0;

    max_a_df_cgs = max(max_a_df_cgs, row.a_df_mag_cgs);
    if (row.dt_df_myr > 0.0)
      min_dt_df_myr =
        (min_dt_df_myr < 0.0) ? row.dt_df_myr :
        min(min_dt_df_myr, row.dt_df_myr);

    row.df_wall_ms = 1000.0 * (ReturnWallTime() - t0);
    total_df_wall_ms += row.df_wall_ms;
    max_bh_df_wall_ms = max(max_bh_df_wall_ms, row.df_wall_ms);

    if (BHDynamicalFrictionVerbose >= 1) {
      if (row.under_resolved)
        fprintf(logptr,
                "[BHDF_WARN] schema_version=1 processor=%d invocation_seq=%d "
                "grid_id=%d step=%d level=%d bh_id=%lld "
                "warning=under_resolved n_dm_in_kernel=%d min_required=%d\n",
                bhdf_processor, bhdf_invocation_seq, bhdf_grid_id,
                cycle_number, level, row.bh_id, row.n_dm_in_kernel,
                BHDynamicalFrictionMinParticles);
      if (no_slow)
        fprintf(logptr,
                "[BHDF_WARN] schema_version=1 processor=%d invocation_seq=%d "
                "grid_id=%d step=%d level=%d bh_id=%lld "
                "warning=no_slow_particles n_slow=%d min_required=%d "
                "v_bh_rel_kms=%.8g\n",
                bhdf_processor, bhdf_invocation_seq, bhdf_grid_id,
                cycle_number, level, row.bh_id, row.n_slow,
                BHDynamicalFrictionMinSlowParticles, row.v_bh_rel_kms);
      if (row.kernel_truncated_by_grid)
        fprintf(logptr,
                "[BHDF_WARN] schema_version=1 processor=%d invocation_seq=%d "
                "grid_id=%d step=%d level=%d bh_id=%lld "
                "warning=kernel_truncated kernel_complete=0 "
                "kernel_truncated_by_grid=1\n",
                bhdf_processor, bhdf_invocation_seq, bhdf_grid_id,
                cycle_number, level, row.bh_id);
      if (row.ln_lambda_zero)
        fprintf(logptr,
                "[BHDF_WARN] schema_version=1 processor=%d invocation_seq=%d "
                "grid_id=%d step=%d level=%d bh_id=%lld "
                "warning=ln_lambda_zero ln_lambda_zero=1 "
                "b_min_phys_kpc=%.8g b_max_phys_kpc=%.8g "
                "rejection_reason=3\n",
                bhdf_processor, bhdf_invocation_seq, bhdf_grid_id,
                cycle_number, level, row.bh_id,
                row.b_min_phys_kpc, row.b_max_phys_kpc);
      if (row.ln_lambda_large)
        fprintf(logptr,
                "[BHDF_WARN] schema_version=1 processor=%d invocation_seq=%d "
                "grid_id=%d step=%d level=%d bh_id=%lld "
                "warning=ln_lambda_large ln_lambda_large=1 ln_lambda=%.8g\n",
                bhdf_processor, bhdf_invocation_seq, bhdf_grid_id,
                cycle_number, level, row.bh_id, row.ln_lambda);
      if (row.cap_triggered)
        fprintf(logptr,
                "[BHDF_WARN] schema_version=1 processor=%d invocation_seq=%d "
                "grid_id=%d step=%d level=%d bh_id=%lld "
                "warning=cap_triggered velocity_kick_cap_fraction=%.8g "
                "cap_triggered=1 dt_code=%.8g\n",
                bhdf_processor, bhdf_invocation_seq, bhdf_grid_id,
                cycle_number, level, row.bh_id,
                row.velocity_kick_cap_fraction, dt_code);
    }

    if (BHDynamicalFrictionVerbose >= 2) {
      BHDFEmitDetail(logptr, row, detail_rows_emitted + 1,
                     detail_rows_complete_flag);
      detail_rows_emitted++;
    }
  }

  const int detail_rows_complete =
    (BHDynamicalFrictionVerbose >= 2 && detail_rows_emitted == rows) ? 1 : 0;

  if (BHDynamicalFrictionVerbose >= 1) {
    fprintf(logptr,
            "[BHDF_SUMMARY] schema_version=1 processor=%d invocation_seq=%d "
            "grid_id=%d step=%d level=%d grid_rank=%d "
            "unsupported_grid_rank=0 active_reposition_enabled=%d "
            "rows=%d eligible_bh=%d valid_or_conditioned_rows=%d "
            "diagnostics_only_rows=%d active_rows=0 applied_kicks=0 "
            "skipped_child_covered=%d skipped_newly_seeded=%d "
            "skipped_non_authoritative_mass=%d no_local_neighbors_count=%d "
            "under_resolved_count=%d kernel_truncated_by_grid_count=%d "
            "no_slow_particles_count=%d ln_lambda_zero_count=%d "
            "ln_lambda_large_count=%d max_a_DF_cgs=%.8e "
            "min_dt_DF_Myr=%.8g total_df_wall_ms=%.4f "
            "max_bh_df_wall_ms=%.4f detail_rows_emitted=%d "
            "detail_rows_complete=%d\n",
            bhdf_processor, bhdf_invocation_seq, bhdf_grid_id,
            cycle_number, level, GridRank, active_reposition_enabled,
            rows, eligible_bh, valid_or_conditioned_rows,
            diagnostics_only_rows, skipped_child_covered,
            skipped_newly_seeded, skipped_non_authoritative_mass,
            no_local_neighbors_count, under_resolved_count,
            kernel_truncated_by_grid_count, no_slow_particles_count,
            ln_lambda_zero_count, ln_lambda_large_count, max_a_df_cgs,
            min_dt_df_myr, total_df_wall_ms, max_bh_df_wall_ms,
            detail_rows_emitted, detail_rows_complete);
  }

  return SUCCESS;
}
