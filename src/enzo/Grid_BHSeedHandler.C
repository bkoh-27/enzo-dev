/***********************************************************************
/
/  GRID CLASS (BH SEED GLOBAL STATE + LEVEL GATHER/FINALIZE)
/
/  PURPOSE:
/    Maintains global MBH cache and diagnostics for BH seeding.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <vector>
#include <algorithm>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "StarParticleData.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

struct BHSeedGlobalState {
  int Active;
  int Level;
  FLOAT Time;
  FLOAT A;
  float Zred;

  FLOAT DomainLeft[MAX_DIMENSION];
  FLOAT DomainWidth[MAX_DIMENSION];
  float BoxKpcH;

  float ExclusionRadiusComKpcH;
  float ExclusionRadius2;
  float CellSizeCode;
  int NBins[MAX_DIMENSION];
  int NBinsPrev[MAX_DIMENSION];
  int NBHPrev;

  std::vector<FLOAT> X;
  std::vector<FLOAT> Y;
  std::vector<FLOAT> Z;
  std::vector<PINT> ID;

  std::vector<int> Head;
  std::vector<int> Next;

  std::vector<FLOAT> NewSeedsX;
  std::vector<FLOAT> NewSeedsY;
  std::vector<FLOAT> NewSeedsZ;

  int BestLocalHas;
  float BestLocalDensity;
  FLOAT BestLocalX;
  FLOAT BestLocalY;
  FLOAT BestLocalZ;
  int BestLocalFlatIndex;
  grid *BestLocalGridPtr;
  int BestLocalDensNum;
  int BestLocalVel1Num;
  int BestLocalVel2Num;
  int BestLocalVel3Num;
  float BestLocalBHMassCode;
  float BestLocalCellWidthCode;

  long long NCandLocal;
  long long DiagLocal[6];
  long long MassGateLocal;
  long long DistBlockedLocal;
  long long CreatedLocal;
  int PreCacheBH;
};

static BHSeedGlobalState BHSeedState = {0, INT_UNDEFINED, FLOAT_UNDEFINED, 1.0,
                                        0.0f,
                                        {0, 0, 0}, {0, 0, 0}, 0.0f,
                                        0.0f, 0.0f, 0.0f,
                                        {1, 1, 1},
                                        {0, 0, 0}, -1,
                                        std::vector<FLOAT>(),
                                        std::vector<FLOAT>(),
                                        std::vector<FLOAT>(),
                                        std::vector<PINT>(),
                                        std::vector<int>(),
                                        std::vector<int>(),
                                        std::vector<FLOAT>(),
                                        std::vector<FLOAT>(),
                                        std::vector<FLOAT>(),
                                        FALSE, -FLT_MAX,
                                        0, 0, 0,
                                        -1, NULL, -1, -1, -1, -1,
                                        0.0f, 0.0f,
                                        0, {0, 0, 0, 0, 0, 0}, 0, 0, 0, 0};

static int BHSeedStepCounter = 0;
static int BHSeedCacheInitialized = FALSE;

static int BHSeedWrapBin(int b, int n)
{
  while (b < 0)
    b += n;
  while (b >= n)
    b -= n;
  return b;
}

static int BHSeedBin1D(FLOAT x, FLOAT left, FLOAT width, float cell_size, int nbins)
{
  FLOAT dx = x - left;
  while (dx < 0)
    dx += width;
  while (dx >= width)
    dx -= width;

  int b = int(floor(dx / cell_size));
  if (b < 0)
    b = 0;
  if (b >= nbins)
    b = nbins - 1;
  return b;
}

static void BHSeedResetDiagnostics()
{
  BHSeedState.NCandLocal = 0;
  for (int n = 0; n < 6; n++)
    BHSeedState.DiagLocal[n] = 0;
  BHSeedState.MassGateLocal = 0;
  BHSeedState.DistBlockedLocal = 0;
  BHSeedState.CreatedLocal = 0;
  BHSeedState.NewSeedsX.clear();
  BHSeedState.NewSeedsY.clear();
  BHSeedState.NewSeedsZ.clear();
  BHSeedState.BestLocalHas = FALSE;
  BHSeedState.BestLocalDensity = -FLT_MAX;
  BHSeedState.BestLocalX = 0;
  BHSeedState.BestLocalY = 0;
  BHSeedState.BestLocalZ = 0;
  BHSeedState.BestLocalFlatIndex = -1;
  BHSeedState.BestLocalGridPtr = NULL;
  BHSeedState.BestLocalDensNum = -1;
  BHSeedState.BestLocalVel1Num = -1;
  BHSeedState.BestLocalVel2Num = -1;
  BHSeedState.BestLocalVel3Num = -1;
  BHSeedState.BestLocalBHMassCode = 0.0f;
  BHSeedState.BestLocalCellWidthCode = 0.0f;
}

static void BHSeedClearCache()
{
  BHSeedState.X.clear();
  BHSeedState.Y.clear();
  BHSeedState.Z.clear();
  BHSeedState.ID.clear();
  BHSeedState.Head.clear();
  BHSeedState.Next.clear();
  BHSeedState.NBinsPrev[0] = 0;
  BHSeedState.NBinsPrev[1] = 0;
  BHSeedState.NBinsPrev[2] = 0;
  BHSeedState.NBHPrev = -1;
}

static void BHSeedBuildLinkedCell()
{
  const int kMaxBinsPerDim = 1 << 20;
  const long double kMaxTotalBins = 4.0e6L;

  for (int dim = 0; dim < MAX_DIMENSION; dim++)
    BHSeedState.NBins[dim] = 1;

  if (BHSeedState.CellSizeCode > 0.0f) {
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
      const double raw_bins = floor(BHSeedState.DomainWidth[dim] / BHSeedState.CellSizeCode);
      if (!isfinite(raw_bins) || raw_bins < 1.0)
        BHSeedState.NBins[dim] = 1;
      else if (raw_bins > double(kMaxBinsPerDim))
        BHSeedState.NBins[dim] = kMaxBinsPerDim;
      else
        BHSeedState.NBins[dim] = max(1, int(raw_bins));
    }
  }

  long double total_bins_ld =
    ((long double) BHSeedState.NBins[0]) *
    ((long double) BHSeedState.NBins[1]) *
    ((long double) BHSeedState.NBins[2]);

  while (total_bins_ld > kMaxTotalBins) {
    int max_dim = 0;
    if (BHSeedState.NBins[1] > BHSeedState.NBins[max_dim]) max_dim = 1;
    if (BHSeedState.NBins[2] > BHSeedState.NBins[max_dim]) max_dim = 2;
    if (BHSeedState.NBins[max_dim] <= 1)
      break;
    BHSeedState.NBins[max_dim] = max(1, BHSeedState.NBins[max_dim] / 2);

    total_bins_ld =
      ((long double) BHSeedState.NBins[0]) *
      ((long double) BHSeedState.NBins[1]) *
      ((long double) BHSeedState.NBins[2]);
  }

  int nbh_now = int(BHSeedState.X.size());
  if (BHSeedState.NBins[0] == BHSeedState.NBinsPrev[0] &&
      BHSeedState.NBins[1] == BHSeedState.NBinsPrev[1] &&
      BHSeedState.NBins[2] == BHSeedState.NBinsPrev[2] &&
      nbh_now == BHSeedState.NBHPrev)
    return;

  BHSeedState.NBinsPrev[0] = BHSeedState.NBins[0];
  BHSeedState.NBinsPrev[1] = BHSeedState.NBins[1];
  BHSeedState.NBinsPrev[2] = BHSeedState.NBins[2];
  BHSeedState.NBHPrev = nbh_now;

  size_t total_bins =
    size_t(BHSeedState.NBins[0]) *
    size_t(BHSeedState.NBins[1]) *
    size_t(BHSeedState.NBins[2]);

  if (total_bins > size_t(INT_MAX))
    total_bins = size_t(INT_MAX);

  BHSeedState.Head.assign(total_bins, -1);
  BHSeedState.Next.assign(BHSeedState.X.size(), -1);

  if (BHSeedState.CellSizeCode <= 0.0f)
    return;

  for (int n = 0; n < int(BHSeedState.X.size()); n++) {
    int bi = BHSeedBin1D(BHSeedState.X[n], BHSeedState.DomainLeft[0],
                         BHSeedState.DomainWidth[0], BHSeedState.CellSizeCode,
                         BHSeedState.NBins[0]);
    int bj = BHSeedBin1D(BHSeedState.Y[n], BHSeedState.DomainLeft[1],
                         BHSeedState.DomainWidth[1], BHSeedState.CellSizeCode,
                         BHSeedState.NBins[1]);
    int bk = BHSeedBin1D(BHSeedState.Z[n], BHSeedState.DomainLeft[2],
                         BHSeedState.DomainWidth[2], BHSeedState.CellSizeCode,
                         BHSeedState.NBins[2]);
    int bin = (bk*BHSeedState.NBins[1] + bj)*BHSeedState.NBins[0] + bi;
    BHSeedState.Next[n] = BHSeedState.Head[bin];
    BHSeedState.Head[bin] = n;
  }
}

static void BHSeedAppendGlobalBH(FLOAT x, FLOAT y, FLOAT z)
{
  int n = BHSeedState.X.size();

  BHSeedState.X.push_back(x);
  BHSeedState.Y.push_back(y);
  BHSeedState.Z.push_back(z);
  BHSeedState.ID.push_back(INT_UNDEFINED);
  BHSeedState.Next.push_back(-1);

  if (BHSeedState.CellSizeCode <= 0.0f)
    return;

  int bi = BHSeedBin1D(x, BHSeedState.DomainLeft[0], BHSeedState.DomainWidth[0],
                       BHSeedState.CellSizeCode, BHSeedState.NBins[0]);
  int bj = BHSeedBin1D(y, BHSeedState.DomainLeft[1], BHSeedState.DomainWidth[1],
                       BHSeedState.CellSizeCode, BHSeedState.NBins[1]);
  int bk = BHSeedBin1D(z, BHSeedState.DomainLeft[2], BHSeedState.DomainWidth[2],
                       BHSeedState.CellSizeCode, BHSeedState.NBins[2]);
  int bin = (bk*BHSeedState.NBins[1] + bj)*BHSeedState.NBins[0] + bi;

  BHSeedState.Next[n] = BHSeedState.Head[bin];
  BHSeedState.Head[bin] = n;
}

static void BHSeedAppendGlobalBHUnique(FLOAT x, FLOAT y, FLOAT z)
{
  for (int n = 0; n < int(BHSeedState.X.size()); n++)
    if (BHSeedState.X[n] == x && BHSeedState.Y[n] == y && BHSeedState.Z[n] == z)
      return;

  BHSeedAppendGlobalBH(x, y, z);
}

static int BHSeedPositionBlocked(FLOAT xpos, FLOAT ypos, FLOAT zpos)
{
  if (BHSeedState.CellSizeCode <= 0.0f || BHSeedState.X.empty())
    return FALSE;

  const float max_periodic_dist2 =
    0.25f * BHSeedState.BoxKpcH * BHSeedState.BoxKpcH *
    (BHSeedState.DomainWidth[0]*BHSeedState.DomainWidth[0] +
     BHSeedState.DomainWidth[1]*BHSeedState.DomainWidth[1] +
     BHSeedState.DomainWidth[2]*BHSeedState.DomainWidth[2]);
  if (BHSeedState.ExclusionRadius2 >= max_periodic_dist2)
    return TRUE;

  int bi = BHSeedBin1D(xpos, BHSeedState.DomainLeft[0], BHSeedState.DomainWidth[0],
                       BHSeedState.CellSizeCode, BHSeedState.NBins[0]);
  int bj = BHSeedBin1D(ypos, BHSeedState.DomainLeft[1], BHSeedState.DomainWidth[1],
                       BHSeedState.CellSizeCode, BHSeedState.NBins[1]);
  int bk = BHSeedBin1D(zpos, BHSeedState.DomainLeft[2], BHSeedState.DomainWidth[2],
                       BHSeedState.CellSizeCode, BHSeedState.NBins[2]);

  for (int dk = -1; dk <= 1; dk++)
    for (int dj = -1; dj <= 1; dj++)
      for (int di = -1; di <= 1; di++) {
        int nbi = BHSeedWrapBin(bi + di, BHSeedState.NBins[0]);
        int nbj = BHSeedWrapBin(bj + dj, BHSeedState.NBins[1]);
        int nbk = BHSeedWrapBin(bk + dk, BHSeedState.NBins[2]);
        int bin = (nbk*BHSeedState.NBins[1] + nbj)*BHSeedState.NBins[0] + nbi;

        for (int n = BHSeedState.Head[bin]; n >= 0; n = BHSeedState.Next[n]) {
          float dx = xpos - BHSeedState.X[n];
          float dy = ypos - BHSeedState.Y[n];
          float dz = zpos - BHSeedState.Z[n];

          if (dx >  0.5f*BHSeedState.DomainWidth[0]) dx -= BHSeedState.DomainWidth[0];
          if (dx < -0.5f*BHSeedState.DomainWidth[0]) dx += BHSeedState.DomainWidth[0];
          if (dy >  0.5f*BHSeedState.DomainWidth[1]) dy -= BHSeedState.DomainWidth[1];
          if (dy < -0.5f*BHSeedState.DomainWidth[1]) dy += BHSeedState.DomainWidth[1];
          if (dz >  0.5f*BHSeedState.DomainWidth[2]) dz -= BHSeedState.DomainWidth[2];
          if (dz < -0.5f*BHSeedState.DomainWidth[2]) dz += BHSeedState.DomainWidth[2];

          dx *= BHSeedState.BoxKpcH;
          dy *= BHSeedState.BoxKpcH;
          dz *= BHSeedState.BoxKpcH;

          float dist2 = dx*dx + dy*dy + dz*dz;
          if (dist2 < BHSeedState.ExclusionRadius2)
            return TRUE;
        }
      }

  return FALSE;
}

struct BHSeedSortOrder {
  const std::vector<PINT> *ID;

  bool operator()(int ia, int ib) const
  {
    return (*ID)[ia] < (*ID)[ib];
  }
};

int BHSeedIsActive()
{
  return BHSeedState.Active;
}

void BHSeedAccumulateKernelDiagnostics(int ncand, const int diag[6])
{
  BHSeedState.NCandLocal += ncand;
  for (int n = 0; n < 6; n++)
    BHSeedState.DiagLocal[n] += diag[n];
}

void BHSeedAddDistBlocked(long long nblocked)
{
  BHSeedState.DistBlockedLocal += nblocked;
}

void BHSeedAddMassGate(long long ngated)
{
  BHSeedState.MassGateLocal += ngated;
}

int BHSeedCandidateBlocked(FLOAT xpos, FLOAT ypos, FLOAT zpos)
{
  return BHSeedPositionBlocked(xpos, ypos, zpos);
}

void BHSeedAccumulateLocalBest(int has, float density,
                               FLOAT xpos, FLOAT ypos, FLOAT zpos,
                               int flat_index, grid *grid_ptr,
                               int DensNum, int Vel1Num, int Vel2Num, int Vel3Num,
                               float bh_mass_code, float cell_width_code)
{
  if (!has || grid_ptr == NULL || flat_index < 0)
    return;

  if (!BHSeedState.BestLocalHas ||
      density > BHSeedState.BestLocalDensity ||
      (density == BHSeedState.BestLocalDensity &&
       (xpos < BHSeedState.BestLocalX ||
        (xpos == BHSeedState.BestLocalX && ypos < BHSeedState.BestLocalY) ||
        (xpos == BHSeedState.BestLocalX && ypos == BHSeedState.BestLocalY &&
         zpos < BHSeedState.BestLocalZ)))) {
    BHSeedState.BestLocalHas = TRUE;
    BHSeedState.BestLocalDensity = density;
    BHSeedState.BestLocalX = xpos;
    BHSeedState.BestLocalY = ypos;
    BHSeedState.BestLocalZ = zpos;
    BHSeedState.BestLocalFlatIndex = flat_index;
    BHSeedState.BestLocalGridPtr = grid_ptr;
    BHSeedState.BestLocalDensNum = DensNum;
    BHSeedState.BestLocalVel1Num = Vel1Num;
    BHSeedState.BestLocalVel2Num = Vel2Num;
    BHSeedState.BestLocalVel3Num = Vel3Num;
    BHSeedState.BestLocalBHMassCode = bh_mass_code;
    BHSeedState.BestLocalCellWidthCode = cell_width_code;
  }
}

void BHSeedAppendSeedToGlobalCache(FLOAT xpos, FLOAT ypos, FLOAT zpos)
{
  BHSeedAppendGlobalBHUnique(xpos, ypos, zpos);
}

void BHSeedRecordCreatedSeed(FLOAT xpos, FLOAT ypos, FLOAT zpos)
{
  BHSeedState.CreatedLocal++;
  BHSeedState.NewSeedsX.push_back(xpos);
  BHSeedState.NewSeedsY.push_back(ypos);
  BHSeedState.NewSeedsZ.push_back(zpos);
  BHSeedAppendGlobalBHUnique(xpos, ypos, zpos);
}

int BHSeedCreateLocalBestParticle()
{
  if (!BHSeedState.BestLocalHas || BHSeedState.BestLocalGridPtr == NULL ||
      BHSeedState.BestLocalFlatIndex < 0)
    return SUCCESS;

  grid *GridData = BHSeedState.BestLocalGridPtr;
  int index = BHSeedState.BestLocalFlatIndex;
  int DensNum = BHSeedState.BestLocalDensNum;
  int Vel1Num = BHSeedState.BestLocalVel1Num;
  int Vel2Num = BHSeedState.BestLocalVel2Num;
  int Vel3Num = BHSeedState.BestLocalVel3Num;
  float bh_mass_code = BHSeedState.BestLocalBHMassCode;
  float cell_width = BHSeedState.BestLocalCellWidthCode;

  if (bh_mass_code <= 0.0f || bh_mass_code > GridData->BaryonField[DensNum][index]) {
    BHSeedState.MassGateLocal++;
    /* Particle not created; caller must check CreatedLocal before caching. */
    return SUCCESS;
  }

  GridData->BaryonField[DensNum][index] -= bh_mass_code;

  const int nx = GridData->GridDimension[0];
  const int ny = GridData->GridDimension[1];
  const int xo = 1;
  const int yo = nx;
  const int zo = nx*ny;

  ParticleEntry p;
  p.Position[0] = BHSeedState.BestLocalX;
  p.Position[1] = BHSeedState.BestLocalY;
  p.Position[2] = BHSeedState.BestLocalZ;

  if (HydroMethod == Zeus_Hydro) {
    p.Velocity[0] = 0.5f*(GridData->BaryonField[Vel1Num][index] +
                          GridData->BaryonField[Vel1Num][index+xo]);
    p.Velocity[1] = 0.5f*(GridData->BaryonField[Vel2Num][index] +
                          GridData->BaryonField[Vel2Num][index+yo]);
    p.Velocity[2] = 0.5f*(GridData->BaryonField[Vel3Num][index] +
                          GridData->BaryonField[Vel3Num][index+zo]);
  } else {
    p.Velocity[0] = GridData->BaryonField[Vel1Num][index];
    p.Velocity[1] = GridData->BaryonField[Vel2Num][index];
    p.Velocity[2] = GridData->BaryonField[Vel3Num][index];
  }

  for (int n = 0; n < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES; n++)
    p.Attribute[n] = 0.0f;

  p.Type = PARTICLE_TYPE_MBH;
  p.Number = INT_UNDEFINED;
  /* bh_mass_code is code-density (mass/cell-volume) units.
     AddOneParticleFromList expects code-mass units, so convert by cell volume. */
  p.Mass = bh_mass_code * pow(cell_width, 3);

  if (NumberOfParticleAttributes > 0)
    p.Attribute[0] = GridData->Time;
  if (NumberOfParticleAttributes > 1)
    p.Attribute[1] = 0.0f;
  if (NumberOfParticleAttributes > 2) {
    float metal_fraction = 0.0f;
    int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum;
    int MetalIaNum, MetalIINum;
    if (GridData->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
                                       MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum)
        != FAIL) {
      float den = GridData->BaryonField[DensNum][index];
      if (den > 0.0f) {
        if (MetalNum != -1 && SNColourNum != -1)
          metal_fraction = (GridData->BaryonField[MetalNum][index] +
                            GridData->BaryonField[SNColourNum][index]) / den;
        else if (MetalNum != -1)
          metal_fraction = GridData->BaryonField[MetalNum][index] / den;
        else if (SNColourNum != -1)
          metal_fraction = GridData->BaryonField[SNColourNum][index] / den;
      }
    }
    p.Attribute[2] = metal_fraction;
  }

  if (GridData->AddOneParticleFromList(&p, 0) == FAIL)
    return FAIL;

  BHSeedRecordCreatedSeed(BHSeedState.BestLocalX, BHSeedState.BestLocalY,
                          BHSeedState.BestLocalZ);

  return SUCCESS;
}

int BHSeedBeginLevel(HierarchyEntry *Grids[], int NumberOfGrids, int level, FLOAT time,
                     LevelHierarchyEntry *LevelArray[])
{
  if (!BHSeedingMethod)
    return SUCCESS;

  if (!BHSeedRunEveryTimestep && StarFormationOncePerRootGridTimeStep) {
    int any_make_stars = 0;
    for (int g = 0; g < NumberOfGrids; g++)
      if (Grids[g]->GridData->ReturnMakeStars()) {
        any_make_stars = 1;
        break;
      }
    if (!any_make_stars)
      return SUCCESS;
  }

  BHSeedState.Active = TRUE;
  BHSeedState.Level = level;
  BHSeedState.Time = time;
  BHSeedResetDiagnostics();
  /* Cache is cleared only at level==0 below; finer levels reuse it. */

  FLOAT dadt = 0;
  BHSeedState.A = 1.0;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(time, &BHSeedState.A, &dadt);

  if (ComovingCoordinates)
    BHSeedState.Zred = float((1.0 + InitialRedshift)/BHSeedState.A - 1.0);
  else
    BHSeedState.Zred = 0.0f;

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    BHSeedState.DomainLeft[dim] = DomainLeftEdge[dim];
    BHSeedState.DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
  }

  if (ComovingCoordinates) {
    BHSeedState.BoxKpcH = ComovingBoxSize * 1000.0f;
  } else {
    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1;
    float TimeUnits = 1, VelocityUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                 &TimeUnits, &VelocityUnits, time) == FAIL)
      return FAIL;
    BHSeedState.BoxKpcH = LengthUnits / kpc_cm;
  }

  float h_param = 1.0f;
  float a_param = 1.0f;
  if (ComovingCoordinates) {
    h_param = HubbleConstantNow;
    if (h_param > 10.0f)
      h_param *= 0.01f;
    /* BHSeedState.A is Enzo internal A=a_phys*(1+InitialRedshift). */
    a_param = BHSeedState.A / (1.0f + float(InitialRedshift));
    if (h_param <= 0.0f)
      h_param = 1.0f;
    if (a_param <= 0.0f)
      a_param = 1.0f;
  }

  BHSeedState.ExclusionRadiusComKpcH = BHSeedExclusionRadius * h_param / a_param;
  BHSeedState.ExclusionRadius2 = BHSeedState.ExclusionRadiusComKpcH *
                                 BHSeedState.ExclusionRadiusComKpcH;
  BHSeedState.CellSizeCode = 0.0f;
  if (BHSeedState.ExclusionRadiusComKpcH > 0.0f && BHSeedState.BoxKpcH > 0.0f)
    BHSeedState.CellSizeCode = 0.5f * BHSeedState.ExclusionRadiusComKpcH /
                               BHSeedState.BoxKpcH;

  if (level == 0 && !BHSeedCacheInitialized) {
    BHSeedClearCache();

    std::vector<PINT> local_id;
    std::vector<FLOAT> local_xyz;

    for (int lev = 0; lev <= MaximumRefinementLevel; lev++) {
      for (LevelHierarchyEntry *lhe = LevelArray[lev]; lhe != NULL;
           lhe = lhe->NextGridThisLevel) {
        grid *GridData = lhe->GridData;

        if (GridData->ReturnProcessorNumber() != MyProcessorNumber)
          continue;

        const int np = GridData->ReturnNumberOfParticles();
        if (np <= 0)
          continue;

        ParticleEntry *plist = new ParticleEntry[np];
        int nret = GridData->ReturnParticleEntry(plist);

        for (int p = 0; p < nret; p++)
          if (plist[p].Type == PARTICLE_TYPE_MBH ||
              plist[p].Type == PARTICLE_TYPE_BLACK_HOLE) {
            local_id.push_back(plist[p].Number);
            local_xyz.push_back(plist[p].Position[0]);
            local_xyz.push_back(plist[p].Position[1]);
            local_xyz.push_back(plist[p].Position[2]);
          }

        delete [] plist;
      }
    }

    int local_n = local_id.size();
    int global_n = local_n;

#ifdef USE_MPI
    std::vector<Eint32> counts, displs;
    std::vector<Eint32> counts3, displs3;
    if (NumberOfProcessors > 1) {
      counts.assign(NumberOfProcessors, 0);
      displs.assign(NumberOfProcessors, 0);

      MPI_Allgather(&local_n, 1, MPI_INT, &counts[0], 1, MPI_INT, MPI_COMM_WORLD);

      global_n = 0;
      for (int p = 0; p < NumberOfProcessors; p++) {
        displs[p] = global_n;
        global_n += counts[p];
      }

      counts3.assign(NumberOfProcessors, 0);
      displs3.assign(NumberOfProcessors, 0);
      for (int p = 0; p < NumberOfProcessors; p++) {
        counts3[p] = 3 * counts[p];
        displs3[p] = 3 * displs[p];
      }

      BHSeedState.ID.resize(global_n);
      std::vector<FLOAT> global_xyz(3*global_n);

      MPI_Allgatherv((local_n > 0) ? &local_id[0] : NULL, local_n,
                     PINTDataType,
                     (global_n > 0) ? &BHSeedState.ID[0] : NULL,
                     &counts[0], &displs[0], PINTDataType,
                     MPI_COMM_WORLD);

      MPI_Allgatherv((local_n > 0) ? &local_xyz[0] : NULL, 3*local_n,
                     FLOATDataType,
                     (global_n > 0) ? &global_xyz[0] : NULL,
                     &counts3[0], &displs3[0], FLOATDataType,
                     MPI_COMM_WORLD);

      BHSeedState.X.resize(global_n);
      BHSeedState.Y.resize(global_n);
      BHSeedState.Z.resize(global_n);

      if (global_n > 0) {
        std::vector<int> order(global_n);
        for (int n = 0; n < global_n; n++)
          order[n] = n;

        BHSeedSortOrder sorter;
        sorter.ID = &BHSeedState.ID;
        std::sort(order.begin(), order.end(), sorter);

        std::vector<PINT> sorted_id(global_n);
        for (int n = 0; n < global_n; n++) {
          int old = order[n];
          BHSeedState.X[n] = global_xyz[3*old + 0];
          BHSeedState.Y[n] = global_xyz[3*old + 1];
          BHSeedState.Z[n] = global_xyz[3*old + 2];
          sorted_id[n] = BHSeedState.ID[old];
        }
        BHSeedState.ID.swap(sorted_id);
      }
    } else
#endif
    {
      BHSeedState.ID = local_id;
      BHSeedState.X.resize(local_n);
      BHSeedState.Y.resize(local_n);
      BHSeedState.Z.resize(local_n);
      for (int n = 0; n < local_n; n++) {
        BHSeedState.X[n] = local_xyz[3*n + 0];
        BHSeedState.Y[n] = local_xyz[3*n + 1];
        BHSeedState.Z[n] = local_xyz[3*n + 2];
      }
    }

    BHSeedCacheInitialized = TRUE;
  }

  BHSeedBuildLinkedCell();
  BHSeedState.PreCacheBH = int(BHSeedState.X.size());

  return SUCCESS;
}

int BHSeedFinalizeLevel()
{
  if (!BHSeedingMethod || !BHSeedState.Active)
    return SUCCESS;

  double local5[5];
  if (BHSeedState.BestLocalHas) {
    local5[0] = BHSeedState.BestLocalDensity;
    local5[1] = BHSeedState.BestLocalX;
    local5[2] = BHSeedState.BestLocalY;
    local5[3] = BHSeedState.BestLocalZ;
    local5[4] = MyProcessorNumber;
  } else {
    local5[0] = -DBL_MAX;
    local5[1] = DBL_MAX;
    local5[2] = DBL_MAX;
    local5[3] = DBL_MAX;
    local5[4] = MyProcessorNumber;
  }

  double best_density = local5[0];
  double best_x = local5[1];
  double best_y = local5[2];
  double best_z = local5[3];
  Eint32 best_rank = Eint32(local5[4]);

#ifdef USE_MPI
  if (NumberOfProcessors > 1) {
    std::vector<double> all5(5*NumberOfProcessors, 0.0);
    MPI_Allgather(local5, 5, MPI_DOUBLE, &all5[0], 5, MPI_DOUBLE, MPI_COMM_WORLD);

    best_density = -DBL_MAX;
    best_x = DBL_MAX;
    best_y = DBL_MAX;
    best_z = DBL_MAX;
    best_rank = ROOT_PROCESSOR;

    for (int p = 0; p < NumberOfProcessors; p++) {
      double d = all5[5*p + 0];
      double x = all5[5*p + 1];
      double y = all5[5*p + 2];
      double z = all5[5*p + 3];
      Eint32 r = Eint32(all5[5*p + 4]);

      if (d > best_density ||
          (d == best_density &&
           (x < best_x ||
            (x == best_x && y < best_y) ||
            (x == best_x && y == best_y && z < best_z) ||
            (x == best_x && y == best_y && z == best_z && r < best_rank)))) {
        best_density = d;
        best_x = x;
        best_y = y;
        best_z = z;
        best_rank = r;
      }
    }
  }
#endif

  int seed_created = FALSE;
  if (best_density > -0.5*DBL_MAX) {
    if (MyProcessorNumber == best_rank) {
      long long created_before = BHSeedState.CreatedLocal;
      if (BHSeedCreateLocalBestParticle() == FAIL)
        return FAIL;
      seed_created = (BHSeedState.CreatedLocal > created_before) ? TRUE : FALSE;
    }

#ifdef USE_MPI
    if (NumberOfProcessors > 1)
      MPI_Bcast(&seed_created, 1, MPI_INT, best_rank, MPI_COMM_WORLD);
#endif

    if (seed_created)
      BHSeedAppendSeedToGlobalCache(best_x, best_y, best_z);
  }

  long long local[10], global[10];
  local[0] = BHSeedState.NCandLocal;
  for (int n = 0; n < 6; n++)
    local[n+1] = BHSeedState.DiagLocal[n];
  local[7] = BHSeedState.MassGateLocal;
  local[8] = BHSeedState.DistBlockedLocal;
  local[9] = BHSeedState.CreatedLocal;

#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Reduce(local, global, 10, MPI_LONG_LONG_INT, MPI_SUM,
               ROOT_PROCESSOR, MPI_COMM_WORLD);
  else
#endif
    for (int i = 0; i < 10; i++)
      global[i] = local[i];

  long long ncand_min = BHSeedState.NCandLocal;
  long long ncand_max = BHSeedState.NCandLocal;
#ifdef USE_MPI
  if (NumberOfProcessors > 1) {
    MPI_Reduce(&BHSeedState.NCandLocal, &ncand_min, 1, MPI_LONG_LONG_INT, MPI_MIN,
               ROOT_PROCESSOR, MPI_COMM_WORLD);
    MPI_Reduce(&BHSeedState.NCandLocal, &ncand_max, 1, MPI_LONG_LONG_INT, MPI_MAX,
               ROOT_PROCESSOR, MPI_COMM_WORLD);
  }
#endif

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    BHSeedStepCounter++;
    FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;
    fprintf(logptr,
            "[BHSEED] step=%d level=%d z=%.4f a_phys=%.6f "
            "excl_phys_kpc=%.1f excl_com_kpch=%.1f "
            "cell_code=%.6e nbins=%d,%d,%d "
            "ncand_local_min=%lld ncand_local_max=%lld ncand_global=%lld "
            "ngates_density=%lld ngates_temp=%lld ngates_metal=%lld "
            "ngates_conv=%lld ngates_cool=%lld ngates_bound=%lld "
            "ngates_mass=%lld dist_blocked=%lld created=%lld total_mbh=%d pre_cache_bh=%d\n",
            BHSeedStepCounter, BHSeedState.Level, BHSeedState.Zred,
            double(BHSeedState.A)/(1.0 + InitialRedshift),
            BHSeedExclusionRadius, BHSeedState.ExclusionRadiusComKpcH,
            BHSeedState.CellSizeCode,
            BHSeedState.NBins[0], BHSeedState.NBins[1], BHSeedState.NBins[2],
            ncand_min, ncand_max, global[0],
            global[1], global[2], global[3], global[4], global[5], global[6],
            global[7], global[8], global[9], int(BHSeedState.X.size()),
            BHSeedState.PreCacheBH);
  }

  BHSeedState.Active = FALSE;

  return SUCCESS;
}
