/***********************************************************************
/
/  MBH MAKER 2 (standalone BH seeding candidate caller)
/
/  PURPOSE:
/    Run BH candidate kernel, apply linked-cell distance exclusion,
/    evaluate Phase 2 spherical kernel quantities, enforce deterministic
/    winner selection across MPI ranks, and create one MBH.
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <vector>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "Star.h"
#include "StarParticleData.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

int star_maker_bh_seed(int *nx, int *ny, int *nz, int *ibuff, int *imethod,
                       float *d, float *dm, float *temp, float *u, float *v, float *w,
                       float *cooltime, float *r, float *metal,
                       float *dx, float *d1, float *t1,
                       float *odthresh, float *metalthresh, float *tempthresh,
                       int *veldivcrit, int *thermalcrit, int *selfboundcrit,
                       int *requirefinestlevel, int *requirelocalpeak,
                       int *ncand, int *cand_index, float *cand_density,
                       int *diag);

int BHSeedIsActive();
void BHSeedAccumulateKernelDiagnostics(int ncand, const int diag[8]);
void BHSeedAddDistBlocked(long long nblocked);
void BHSeedAddMassGate(long long ngated);
void BHSeedAddLegacyCellMassWouldFail(long long nshadow);
void BHSeedAddKernelEvaluationDiagnostics(long long nevaluated,
                                          long long ntruncated,
                                          long long ngenclosedmass);
void BHSeedNoteDMFieldAvailability(int available);
int BHSeedCandidateBlocked(FLOAT xpos, FLOAT ypos, FLOAT zpos);
void BHSeedAccumulateLocalBest(int has, float density,
                               FLOAT xpos, FLOAT ypos, FLOAT zpos,
                               int flat_index, grid *grid_ptr,
                               int DensNum, int Vel1Num, int Vel2Num, int Vel3Num,
                               float bh_mass_code, float cell_width_code,
                               float patch_mass_msun,
                               float patch_metallicity,
                               float patch_density_peak,
                               int kernel_complete,
                               float host_dm_density);

struct BHSeedKernelEval {
  float patch_mass_msun;
  float patch_metallicity;
  float patch_density_peak;
  int kernel_complete;
};

static void BHSeedEvaluateKernel(grid *ThisGrid,
                                 int index,
                                 int nx, int ny, int nz,
                                 const float *density,
                                 const float *metal_fraction,
                                 float cell_width,
                                 double cell_mass_msun_unit,
                                 float patch_radius,
                                 BHSeedKernelEval *kernel)
{
  const int i0 = index % nx;
  const int j0 = (index / nx) % ny;
  const int k0 = index / (nx*ny);

  const int isx = ThisGrid->GetGridStartIndex(0);
  const int isy = ThisGrid->GetGridStartIndex(1);
  const int isz = ThisGrid->GetGridStartIndex(2);
  const int iex = ThisGrid->GetGridEndIndex(0);
  const int iey = ThisGrid->GetGridEndIndex(1);
  const int iez = ThisGrid->GetGridEndIndex(2);

  const FLOAT xleft = ThisGrid->GetGridLeftEdge(0) - isx*cell_width;
  const FLOAT yleft = ThisGrid->GetGridLeftEdge(1) - isy*cell_width;
  const FLOAT zleft = ThisGrid->GetGridLeftEdge(2) - isz*cell_width;
  const FLOAT xright = ThisGrid->GetGridRightEdge(0) + (nx-1-iex)*cell_width;
  const FLOAT yright = ThisGrid->GetGridRightEdge(1) + (ny-1-iey)*cell_width;
  const FLOAT zright = ThisGrid->GetGridRightEdge(2) + (nz-1-iez)*cell_width;

  const FLOAT x0 = xleft + (FLOAT(i0) + 0.5f)*cell_width;
  const FLOAT y0 = yleft + (FLOAT(j0) + 0.5f)*cell_width;
  const FLOAT z0 = zleft + (FLOAT(k0) + 0.5f)*cell_width;

  double min_dist = x0 - xleft;
  min_dist = min(min_dist, double(xright - x0));
  min_dist = min(min_dist, double(y0 - yleft));
  min_dist = min(min_dist, double(yright - y0));
  min_dist = min(min_dist, double(z0 - zleft));
  min_dist = min(min_dist, double(zright - z0));

  kernel->kernel_complete = (patch_radius <= min_dist) ? 1 : 0;

  const int rcell = max(0, int(ceil(patch_radius / cell_width)));
  const int ilo = max(0, i0 - rcell);
  const int ihi = min(nx - 1, i0 + rcell);
  const int jlo = max(0, j0 - rcell);
  const int jhi = min(ny - 1, j0 + rcell);
  const int klo = max(0, k0 - rcell);
  const int khi = min(nz - 1, k0 + rcell);
  const double r2 = double(patch_radius) * double(patch_radius);

  double patch_mass = 0.0;
  double patch_metal_mass = 0.0;
  float patch_density_peak = density[index];

  /* Deterministic k-j-i loop order for MPI reproducibility.
     Ghost values may be stale by one hydro step; this is accepted for Phase 2. */
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
        const float rho = density[nindex];
        const double cell_mass = rho * cell_mass_msun_unit;
        patch_mass += cell_mass;
        if (metal_fraction != NULL)
          patch_metal_mass += cell_mass * metal_fraction[nindex];
        if (rho > patch_density_peak)
          patch_density_peak = rho;
      }
    }
  }

  kernel->patch_mass_msun = float(patch_mass);
  kernel->patch_density_peak = patch_density_peak;
  if (patch_mass > 0.0 && metal_fraction != NULL)
    kernel->patch_metallicity = float(patch_metal_mass / patch_mass);
  else
    kernel->patch_metallicity = 0.0f;
}

int mbh_maker2(grid *ThisGrid,
               int level,
               int DensNum, int Vel1Num, int Vel2Num, int Vel3Num,
               float *temperature, float *cooling_time, float *dmfield,
               float *metal_fraction, int metal_field_present,
               int dm_field_present,
               float CellWidthTemp,
               float DensityUnits, float LengthUnits, float TimeUnits,
               int MaximumRefinementLevel)
{
  (void) level;
  (void) MaximumRefinementLevel;

  if (!BHSeedingMethod || !BHSeedIsActive())
    return SUCCESS;

  int nx = ThisGrid->GridDimension[0];
  int ny = ThisGrid->GridDimension[1];
  int nz = ThisGrid->GridDimension[2];
  int size = nx*ny*nz;

  /* Not OpenMP-safe; safe because star formation/BH seeding grid loop is serial. */
  static std::vector<int> cand_index_scratch;
  static std::vector<float> cand_density_scratch;
  if (int(cand_index_scratch.size()) < size) {
    cand_index_scratch.resize(size);
    cand_density_scratch.resize(size);
  }
  int *cand_index = &cand_index_scratch[0];
  float *cand_density = &cand_density_scratch[0];
  int diag[8];
  int ncand_local = 0;
  int ibuff = NumberOfGhostZones;
  int imethod = HydroMethod;
  float *metal_ptr = (metal_field_present && metal_fraction != NULL) ?
    metal_fraction : NULL;

  if (star_maker_bh_seed(&nx, &ny, &nz, &ibuff, &imethod,
                         ThisGrid->BaryonField[DensNum], dmfield, temperature,
                         ThisGrid->BaryonField[Vel1Num], ThisGrid->BaryonField[Vel2Num],
                         ThisGrid->BaryonField[Vel3Num], cooling_time,
                         ThisGrid->BaryonField[ThisGrid->NumberOfBaryonFields], metal_ptr,
                         &CellWidthTemp, &DensityUnits, &TimeUnits,
                         &BHSeedOverdensityThreshold,
                         &BHSeedMetallicityThreshold,
                         &BHSeedTemperatureThreshold,
                         &BHSeedVelDivCrit, &BHSeedThermalCrit,
                         &BHSeedSelfBoundCrit,
                         &BHSeedRequireFinestLevel, &BHSeedRequireLocalPeak,
                         &ncand_local, cand_index,
                         cand_density, diag) == FAIL) {
    return FAIL;
  }

  BHSeedAccumulateKernelDiagnostics(ncand_local, diag);
  BHSeedNoteDMFieldAvailability(dm_field_present ? TRUE : FALSE);

  const double cell_mass_msun_unit =
    DensityUnits * pow(LengthUnits*CellWidthTemp, 3) / SolarMass;
  const float bh_mass_code =
    (cell_mass_msun_unit > 0.0) ? float(BHSeedMass / cell_mass_msun_unit) : 0.0f;

  long long dist_blocked = 0;
  long long mass_gated = 0;
  long long legacy_cellmass_would_fail = 0;
  long long nkernel_evaluated = 0;
  long long nkernel_truncated = 0;
  long long enclosed_mass_gated = 0;

  int local_has = FALSE;
  int local_index = -1;
  float local_density = -FLT_MAX;
  FLOAT local_x = 0, local_y = 0, local_z = 0;
  float local_patch_mass = -1.0f;
  float local_patch_metallicity = -1.0f;
  float local_patch_density_peak = -1.0f;
  int local_kernel_complete = -1;
  float local_host_dm_density = -1.0f;

  for (int c = 0; c < ncand_local; c++) {
    int index = cand_index[c];
    int i = index % nx;
    int j = (index / nx) % ny;
    int k = index / (nx*ny);

    FLOAT xpos = ThisGrid->CellLeftEdge[0][0] + (FLOAT(i) + 0.5f)*CellWidthTemp;
    FLOAT ypos = ThisGrid->CellLeftEdge[1][0] + (FLOAT(j) + 0.5f)*CellWidthTemp;
    FLOAT zpos = ThisGrid->CellLeftEdge[2][0] + (FLOAT(k) + 0.5f)*CellWidthTemp;

    const float cell_mass_msun = ThisGrid->BaryonField[DensNum][index] * cell_mass_msun_unit;
    if (cell_mass_msun < BHSeedMass) {
      legacy_cellmass_would_fail++;
      /* Phase 2: legacy cell-mass gate is shadow-only by default.
         See BHSeedLegacyCellMassGate. */
      if (BHSeedLegacyCellMassGate == 1) {
        mass_gated++;
        continue;
      }
    }

    if (BHSeedCandidateBlocked(xpos, ypos, zpos)) {
      dist_blocked++;
      continue;
    }

    BHSeedKernelEval kernel;
    BHSeedEvaluateKernel(ThisGrid, index, nx, ny, nz,
                         ThisGrid->BaryonField[DensNum], metal_ptr,
                         CellWidthTemp, cell_mass_msun_unit,
                         BHSeedPatchRadius, &kernel);
    nkernel_evaluated++;
    if (kernel.kernel_complete == 0)
      nkernel_truncated++;

    if (kernel.patch_mass_msun < BHSeedMinEnclosedMass) {
      enclosed_mass_gated++;
      continue;
    }

    local_has = TRUE;
    local_index = index;
    local_density = cand_density[c];
    local_x = xpos;
    local_y = ypos;
    local_z = zpos;
    local_patch_mass = kernel.patch_mass_msun;
    local_patch_metallicity = kernel.patch_metallicity;
    /* kernel maximum, not candidate-cell density */
    local_patch_density_peak = kernel.patch_density_peak;
    local_kernel_complete = kernel.kernel_complete;
    if (dm_field_present && dmfield != NULL)
      local_host_dm_density = dmfield[index];
    else
      local_host_dm_density = -1.0f;
    break;
  }

  BHSeedAddDistBlocked(dist_blocked);
  BHSeedAddMassGate(mass_gated);
  BHSeedAddLegacyCellMassWouldFail(legacy_cellmass_would_fail);
  BHSeedAddKernelEvaluationDiagnostics(nkernel_evaluated,
                                       nkernel_truncated,
                                       enclosed_mass_gated);
  BHSeedAccumulateLocalBest(local_has, local_density, local_x, local_y, local_z,
                            local_index, ThisGrid,
                            DensNum, Vel1Num, Vel2Num, Vel3Num,
                            bh_mass_code, CellWidthTemp,
                            local_patch_mass, local_patch_metallicity,
                            local_patch_density_peak, local_kernel_complete,
                            local_host_dm_density);

  return SUCCESS;
}
