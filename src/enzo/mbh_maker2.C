/***********************************************************************
/
/  MBH MAKER 2 (standalone BH seeding candidate caller)
/
/  PURPOSE:
/    Run BH candidate kernel, apply linked-cell distance exclusion, enforce
/    deterministic winner selection across MPI ranks, and create one MBH.
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
                       int *ncand, int *cand_index, float *cand_density,
                       int *diag);

int BHSeedIsActive();
void BHSeedAccumulateKernelDiagnostics(int ncand, const int diag[6]);
void BHSeedAddDistBlocked(long long nblocked);
void BHSeedAddMassGate(long long ngated);
int BHSeedCandidateBlocked(FLOAT xpos, FLOAT ypos, FLOAT zpos);
void BHSeedAccumulateLocalBest(int has, float density,
                               FLOAT xpos, FLOAT ypos, FLOAT zpos,
                               int flat_index, grid *grid_ptr,
                               int DensNum, int Vel1Num, int Vel2Num, int Vel3Num,
                               float bh_mass_code, float cell_width_code);

int mbh_maker2(grid *ThisGrid,
               int level,
               int DensNum, int Vel1Num, int Vel2Num, int Vel3Num,
               float *temperature, float *cooling_time, float *dmfield,
               float *metal_fraction, int metal_field_present,
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
  int diag[6];
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
                         &BHSeedSelfBoundCrit, &ncand_local, cand_index,
                         cand_density, diag) == FAIL) {
    return FAIL;
  }

  BHSeedAccumulateKernelDiagnostics(ncand_local, diag);

  const double cell_mass_msun =
    DensityUnits * pow(LengthUnits*CellWidthTemp, 3) / SolarMass;
  const float bh_mass_code =
    (cell_mass_msun > 0.0) ? float(BHSeedMass / cell_mass_msun) : 0.0f;

  long long dist_blocked = 0;
  long long mass_gated = 0;

  int local_has = FALSE;
  int local_index = -1;
  float local_density = -FLT_MAX;
  FLOAT local_x = 0, local_y = 0, local_z = 0;

  for (int c = 0; c < ncand_local; c++) {
    int index = cand_index[c];
    int i = index % nx;
    int j = (index / nx) % ny;
    int k = index / (nx*ny);

    FLOAT xpos = ThisGrid->CellLeftEdge[0][0] + (FLOAT(i) + 0.5f)*CellWidthTemp;
    FLOAT ypos = ThisGrid->CellLeftEdge[1][0] + (FLOAT(j) + 0.5f)*CellWidthTemp;
    FLOAT zpos = ThisGrid->CellLeftEdge[2][0] + (FLOAT(k) + 0.5f)*CellWidthTemp;

    if (BHSeedCandidateBlocked(xpos, ypos, zpos)) {
      dist_blocked++;
      continue;
    }

    if (bh_mass_code <= 0.0f || bh_mass_code > ThisGrid->BaryonField[DensNum][index]) {
      mass_gated++;
      continue;
    }

    local_has = TRUE;
    local_index = index;
    local_density = cand_density[c];
    local_x = xpos;
    local_y = ypos;
    local_z = zpos;
    break;
  }

  BHSeedAddDistBlocked(dist_blocked);
  BHSeedAddMassGate(mass_gated);
  BHSeedAccumulateLocalBest(local_has, local_density, local_x, local_y, local_z,
                            local_index, ThisGrid,
                            DensNum, Vel1Num, Vel2Num, Vel3Num,
                            bh_mass_code, CellWidthTemp);

  return SUCCESS;
}
