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
#include "BHSeedCandidate.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
double ReturnWallTime();

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
  float BestLocalPatchMass;
  float BestLocalPatchMetallicity;
  float BestLocalPatchDensityPeak;
  int BestLocalKernelComplete;
  float BestLocalHostDMDensity;

  long long NCandLocal;
  long long DiagLocal[8];
  long long MassGateLocal;
  long long DistBlockedLocal;
  long long CreatedLocal;
  long long LegacyCellMassWouldFailLocal;
  long long KernelEvaluatedLocal;
  long long KernelTruncatedLocal;
  long long EnclosedMassGateLocal;
  int DMFieldAvailableLocal;
  float LevelCellWidthCode;
  double SeedingWallStart;
  int PreCacheBH;
  float HubbleParam;
  float APhys;
  float CandidateSeparationComKpcH;

  std::vector<BHSeedCandidate> LocalCandidates;
  std::vector<int> LocalGridIDs;
  std::vector<grid*> LocalGridPtrs;
};

static BHSeedGlobalState BHSeedState;

static int BHSeedStepCounter = 0;
static int BHSeedCacheInitialized = FALSE;
static int BHSeedDMAvailabilityLogged = FALSE;

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
  for (int n = 0; n < 8; n++)
    BHSeedState.DiagLocal[n] = 0;
  BHSeedState.MassGateLocal = 0;
  BHSeedState.DistBlockedLocal = 0;
  BHSeedState.CreatedLocal = 0;
  BHSeedState.LegacyCellMassWouldFailLocal = 0;
  BHSeedState.KernelEvaluatedLocal = 0;
  BHSeedState.KernelTruncatedLocal = 0;
  BHSeedState.EnclosedMassGateLocal = 0;
  BHSeedState.DMFieldAvailableLocal = FALSE;
  BHSeedState.NewSeedsX.clear();
  BHSeedState.NewSeedsY.clear();
  BHSeedState.NewSeedsZ.clear();
  BHSeedState.LocalCandidates.clear();
  BHSeedState.LocalGridIDs.clear();
  BHSeedState.LocalGridPtrs.clear();
  BHSeedState.CandidateSeparationComKpcH = 0.0f;
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
  BHSeedState.BestLocalPatchMass = -1.0f;
  BHSeedState.BestLocalPatchMetallicity = -1.0f;
  BHSeedState.BestLocalPatchDensityPeak = -1.0f;
  BHSeedState.BestLocalKernelComplete = -1;
  BHSeedState.BestLocalHostDMDensity = -1.0f;
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

static float BHSeedPeriodicDistance2ComKpcH(FLOAT x1, FLOAT y1, FLOAT z1,
                                             FLOAT x2, FLOAT y2, FLOAT z2)
{
  float dx = x1 - x2;
  float dy = y1 - y2;
  float dz = z1 - z2;

  if (dx >  0.5f*BHSeedState.DomainWidth[0]) dx -= BHSeedState.DomainWidth[0];
  if (dx < -0.5f*BHSeedState.DomainWidth[0]) dx += BHSeedState.DomainWidth[0];
  if (dy >  0.5f*BHSeedState.DomainWidth[1]) dy -= BHSeedState.DomainWidth[1];
  if (dy < -0.5f*BHSeedState.DomainWidth[1]) dy += BHSeedState.DomainWidth[1];
  if (dz >  0.5f*BHSeedState.DomainWidth[2]) dz -= BHSeedState.DomainWidth[2];
  if (dz < -0.5f*BHSeedState.DomainWidth[2]) dz += BHSeedState.DomainWidth[2];

  dx *= BHSeedState.BoxKpcH;
  dy *= BHSeedState.BoxKpcH;
  dz *= BHSeedState.BoxKpcH;

  return dx*dx + dy*dy + dz*dz;
}

static float BHSeedExclusionRadiusComKpcHForCandidate(float dx_local)
{
  float radius_com = 0.0f;

  if (BHSeedExclusionMode == 1) {
    /* Mode 1: parameter is already comoving kpc/h. */
    radius_com = BHSeedExclusionRadius;
  } else if (BHSeedExclusionMode == 2) {
    /* Mode 2: candidate-centric resolution-scaled exclusion radius. */
    radius_com = BHSeedExclusionCells * dx_local * BHSeedState.BoxKpcH;
  } else {
    /* Mode 0: fixed physical kpc converted to comoving kpc/h in BeginLevel. */
    radius_com = BHSeedState.ExclusionRadiusComKpcH;
  }

  if (radius_com < 0.0f)
    radius_com = 0.0f;

  return radius_com;
}

static int BHSeedPositionBlockedWithRadius(FLOAT xpos, FLOAT ypos, FLOAT zpos,
                                           float radius_com_kpch)
{
  if (radius_com_kpch <= 0.0f || BHSeedState.X.empty())
    return FALSE;

  const float radius2 = radius_com_kpch * radius_com_kpch;

  const float max_periodic_dist2 =
    0.25f * BHSeedState.BoxKpcH * BHSeedState.BoxKpcH *
    (BHSeedState.DomainWidth[0]*BHSeedState.DomainWidth[0] +
     BHSeedState.DomainWidth[1]*BHSeedState.DomainWidth[1] +
     BHSeedState.DomainWidth[2]*BHSeedState.DomainWidth[2]);
  if (radius2 >= max_periodic_dist2)
    return TRUE;

  if (BHSeedState.CellSizeCode <= 0.0f || BHSeedState.Head.empty()) {
    for (int n = 0; n < int(BHSeedState.X.size()); n++) {
      float dist2 = BHSeedPeriodicDistance2ComKpcH(xpos, ypos, zpos,
                                                   BHSeedState.X[n],
                                                   BHSeedState.Y[n],
                                                   BHSeedState.Z[n]);
      if (dist2 < radius2)
        return TRUE;
    }
    return FALSE;
  }

  int bi = BHSeedBin1D(xpos, BHSeedState.DomainLeft[0], BHSeedState.DomainWidth[0],
                       BHSeedState.CellSizeCode, BHSeedState.NBins[0]);
  int bj = BHSeedBin1D(ypos, BHSeedState.DomainLeft[1], BHSeedState.DomainWidth[1],
                       BHSeedState.CellSizeCode, BHSeedState.NBins[1]);
  int bk = BHSeedBin1D(zpos, BHSeedState.DomainLeft[2], BHSeedState.DomainWidth[2],
                       BHSeedState.CellSizeCode, BHSeedState.NBins[2]);

  const float bin_width_kpch = BHSeedState.CellSizeCode * BHSeedState.BoxKpcH;
  int nsearch = 0;
  if (bin_width_kpch > 0.0f)
    nsearch = int(ceil(radius_com_kpch / bin_width_kpch));

  for (int dk = -nsearch; dk <= nsearch; dk++)
    for (int dj = -nsearch; dj <= nsearch; dj++)
      for (int di = -nsearch; di <= nsearch; di++) {
        int nbi = BHSeedWrapBin(bi + di, BHSeedState.NBins[0]);
        int nbj = BHSeedWrapBin(bj + dj, BHSeedState.NBins[1]);
        int nbk = BHSeedWrapBin(bk + dk, BHSeedState.NBins[2]);
        int bin = (nbk*BHSeedState.NBins[1] + nbj)*BHSeedState.NBins[0] + nbi;

        for (int n = BHSeedState.Head[bin]; n >= 0; n = BHSeedState.Next[n]) {
          float dist2 = BHSeedPeriodicDistance2ComKpcH(xpos, ypos, zpos,
                                                       BHSeedState.X[n],
                                                       BHSeedState.Y[n],
                                                       BHSeedState.Z[n]);
          if (dist2 < radius2)
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

struct BHSeedCandidateComparator {
  int RankingOrder;
  int UseDeterministicTiebreak;

  bool operator()(const BHSeedCandidate &a, const BHSeedCandidate &b) const
  {
    if (RankingOrder == 1) {
      if (a.patch_density_peak != b.patch_density_peak)
        return a.patch_density_peak > b.patch_density_peak;
      if (a.patch_mass != b.patch_mass)
        return a.patch_mass > b.patch_mass;
    } else {
      if (a.patch_mass != b.patch_mass)
        return a.patch_mass > b.patch_mass;
      if (a.patch_density_peak != b.patch_density_peak)
        return a.patch_density_peak > b.patch_density_peak;
    }

    if (a.patch_metallicity != b.patch_metallicity)
      return a.patch_metallicity < b.patch_metallicity;

    if (UseDeterministicTiebreak) {
      if (a.pos[0] != b.pos[0])
        return a.pos[0] < b.pos[0];
      if (a.pos[1] != b.pos[1])
        return a.pos[1] < b.pos[1];
      if (a.pos[2] != b.pos[2])
        return a.pos[2] < b.pos[2];
    }

    if (a.owning_rank != b.owning_rank)
      return a.owning_rank < b.owning_rank;
    if (a.owning_grid_id != b.owning_grid_id)
      return a.owning_grid_id < b.owning_grid_id;
    if (a.cell_k != b.cell_k)
      return a.cell_k < b.cell_k;
    if (a.cell_j != b.cell_j)
      return a.cell_j < b.cell_j;
    if (a.cell_i != b.cell_i)
      return a.cell_i < b.cell_i;

    return a.cell_index < b.cell_index;
  }
};

int BHSeedIsActive()
{
  return BHSeedState.Active;
}

void BHSeedAccumulateKernelDiagnostics(int ncand, const int diag[8])
{
  BHSeedState.NCandLocal += ncand;
  for (int n = 0; n < 8; n++)
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

void BHSeedAddLegacyCellMassWouldFail(long long nshadow)
{
  BHSeedState.LegacyCellMassWouldFailLocal += nshadow;
}

void BHSeedAddKernelEvaluationDiagnostics(long long nevaluated,
                                          long long ntruncated,
                                          long long ngenclosedmass)
{
  BHSeedState.KernelEvaluatedLocal += nevaluated;
  BHSeedState.KernelTruncatedLocal += ntruncated;
  BHSeedState.EnclosedMassGateLocal += ngenclosedmass;
}

void BHSeedNoteDMFieldAvailability(int available)
{
  if (available)
    BHSeedState.DMFieldAvailableLocal = TRUE;
}

float BHSeedCurrentRedshift()
{
  return BHSeedState.Zred;
}

void BHSeedAddLocalCandidate(const BHSeedCandidate &cand)
{
  BHSeedState.LocalCandidates.push_back(cand);
}

int BHSeedCandidateBlocked(FLOAT xpos, FLOAT ypos, FLOAT zpos, float dx_local)
{
  const float radius_com = BHSeedExclusionRadiusComKpcHForCandidate(dx_local);
  return BHSeedPositionBlockedWithRadius(xpos, ypos, zpos, radius_com);
}

void BHSeedAccumulateLocalBest(int has, float density,
                               FLOAT xpos, FLOAT ypos, FLOAT zpos,
                               int flat_index, grid *grid_ptr,
                               int DensNum, int Vel1Num, int Vel2Num, int Vel3Num,
                               float bh_mass_code, float cell_width_code,
                               float patch_mass_msun,
                               float patch_metallicity,
                               float patch_density_peak,
                               int kernel_complete,
                               float host_dm_density)
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
    BHSeedState.BestLocalPatchMass = patch_mass_msun;
    BHSeedState.BestLocalPatchMetallicity = patch_metallicity;
    BHSeedState.BestLocalPatchDensityPeak = patch_density_peak;
    BHSeedState.BestLocalKernelComplete = kernel_complete;
    BHSeedState.BestLocalHostDMDensity = host_dm_density;
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
  const float seed_density_before = GridData->BaryonField[DensNum][index];

  if (bh_mass_code <= 0.0f) {
    BHSeedState.MassGateLocal++;
    /* Particle not created; caller must check CreatedLocal before caching. */
    return SUCCESS;
  }

  if (bh_mass_code > seed_density_before && BHSeedLegacyCellMassGate == 1) {
    BHSeedState.MassGateLocal++;
    return SUCCESS;
  }

  /* Phase 2: legacy cell-mass gate is shadow-only by default, so creation can
     proceed even when bh_mass_code > seed_density_before. Clamp residual density
     to avoid negative cell values. */
  GridData->BaryonField[DensNum][index] = max(0.0f, seed_density_before - bh_mass_code);

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

  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_CREATION_TIME)
    p.Attribute[PARTICLE_ATTRIBUTE_CREATION_TIME] = GridData->Time;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_DYNAMICAL_TIME)
    p.Attribute[PARTICLE_ATTRIBUTE_DYNAMICAL_TIME] = 0.0f;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_METALLICITY_FRACTION) {
    float metal_fraction = 0.0f;
    int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum;
    int MetalIaNum, MetalIINum;
    if (GridData->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
                                       MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum)
        != FAIL) {
      float den = seed_density_before;
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
    p.Attribute[PARTICLE_ATTRIBUTE_METALLICITY_FRACTION] = metal_fraction;
  }

  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_CHANNEL) {
    /* Stored in float slot; explicit int->float cast for restart-safe roundtrip. */
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_CHANNEL] = float(BHSeedChannel);
  }
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_REDSHIFT)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_REDSHIFT] = BHSeedState.Zred;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_PATCH_MASS)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_PATCH_MASS] = BHSeedState.BestLocalPatchMass;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_PATCH_METALLICITY)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_PATCH_METALLICITY] =
      BHSeedState.BestLocalPatchMetallicity;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_PATCH_DENSITY_PEAK)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_PATCH_DENSITY_PEAK] =
      BHSeedState.BestLocalPatchDensityPeak;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_KERNEL_COMPLETE) {
    /* Int semantic field stored in float slot; explicit int->float cast. */
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_KERNEL_COMPLETE] =
      float(BHSeedState.BestLocalKernelComplete);
  }
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_HOST_DM_DENSITY)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_HOST_DM_DENSITY] =
      BHSeedState.BestLocalHostDMDensity;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_ACCEPT_RANK) {
    /* Int semantic field stored in float slot; -1 means unresolved in Phase 1. */
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_ACCEPT_RANK] = -1.0f; // Phase 3+
  }

  if (GridData->AddOneParticleFromList(&p, 0) == FAIL)
    return FAIL;

  BHSeedRecordCreatedSeed(BHSeedState.BestLocalX, BHSeedState.BestLocalY,
                          BHSeedState.BestLocalZ);
  /* transient flag, not checkpointed, zeroed fresh every pass */
  if (GridData->BaryonField[GridData->NumberOfBaryonFields] != NULL)
    GridData->BaryonField[GridData->NumberOfBaryonFields][index] = 1.0f;

  if (BHSeedVerbose >= 1) {
    FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;
    int seed_channel = (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_CHANNEL) ?
      int(p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_CHANNEL]) : -1;
    int seed_kernel_complete =
      (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_KERNEL_COMPLETE) ?
      int(p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_KERNEL_COMPLETE]) : -1;
    int seed_accept_rank =
      (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_ACCEPT_RANK) ?
      int(p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_ACCEPT_RANK]) : -1;

    fprintf(logptr,
            "[BHSEED_SEED] level=%d x=%.8g y=%.8g z=%.8g channel=%d redshift=%.8g "
            "patch_mass=%.8g patch_metal=%.8g patch_density_peak=%.8g "
            "kernel_complete=%d host_dm_density=%.8g accept_rank=%d\n",
            BHSeedState.Level,
            double(BHSeedState.BestLocalX), double(BHSeedState.BestLocalY),
            double(BHSeedState.BestLocalZ),
            seed_channel,
            (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_REDSHIFT) ?
            double(p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_REDSHIFT]) : -1.0,
            (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_PATCH_MASS) ?
            double(p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_PATCH_MASS]) : -1.0,
            (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_PATCH_METALLICITY) ?
            double(p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_PATCH_METALLICITY]) : -1.0,
            (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_PATCH_DENSITY_PEAK) ?
            double(p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_PATCH_DENSITY_PEAK]) : -1.0,
            seed_kernel_complete,
            (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_HOST_DM_DENSITY) ?
            double(p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_HOST_DM_DENSITY]) : -1.0,
            seed_accept_rank);
  }

  return SUCCESS;
}

static grid* BHSeedFindLocalGridByID(int grid_id)
{
  for (int n = 0; n < int(BHSeedState.LocalGridIDs.size()); n++)
    if (BHSeedState.LocalGridIDs[n] == grid_id)
      return BHSeedState.LocalGridPtrs[n];
  return NULL;
}

int BHSeedCreateAcceptedCandidate(const BHSeedCandidate &cand,
                                         long long *nskipped_insufficient_gas)
{
  if (cand.owning_rank != MyProcessorNumber)
    return SUCCESS;

  grid *GridData = BHSeedFindLocalGridByID(cand.owning_grid_id);
  if (GridData == NULL)
    ENZO_FAIL("BHSeedCreateAcceptedCandidate: owning grid not found on owning rank.");

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  if (GridData->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                           Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL)
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");

  const int nx = GridData->GridDimension[0];
  const int ny = GridData->GridDimension[1];
  const int nz = GridData->GridDimension[2];
  const int xo = 1;
  const int yo = nx;
  const int zo = nx*ny;

  const int i0 = cand.cell_i;
  const int j0 = cand.cell_j;
  const int k0 = cand.cell_k;
  const int index0 = cand.cell_index;
  const float cell_width = cand.dx_local;
  const double cell_volume_code = pow(cell_width, 3.0);
  const double bh_mass_code = double(cand.seed_mass_code_density) * cell_volume_code;

  if (bh_mass_code <= 0.0) {
    (*nskipped_insufficient_gas)++;
    return SUCCESS;
  }

  const int isx = GridData->GetGridStartIndex(0);
  const int isy = GridData->GetGridStartIndex(1);
  const int isz = GridData->GetGridStartIndex(2);
  const int iex = GridData->GetGridEndIndex(0);
  const int iey = GridData->GetGridEndIndex(1);
  const int iez = GridData->GetGridEndIndex(2);

  const int rcell = max(0, int(ceil(BHSeedPatchRadius / cell_width)));
  const int ilo = max(0, i0 - rcell);
  const int ihi = min(nx - 1, i0 + rcell);
  const int jlo = max(0, j0 - rcell);
  const int jhi = min(ny - 1, j0 + rcell);
  const int klo = max(0, k0 - rcell);
  const int khi = min(nz - 1, k0 + rcell);
  const double r2 = double(BHSeedPatchRadius) * double(BHSeedPatchRadius);

  std::vector<int> active_index;
  std::vector<double> active_mass;
  std::vector<double> active_velx;
  std::vector<double> active_vely;
  std::vector<double> active_velz;
  active_index.reserve((2*rcell+1)*(2*rcell+1)*(2*rcell+1));
  active_mass.reserve((2*rcell+1)*(2*rcell+1)*(2*rcell+1));
  active_velx.reserve((2*rcell+1)*(2*rcell+1)*(2*rcell+1));
  active_vely.reserve((2*rcell+1)*(2*rcell+1)*(2*rcell+1));
  active_velz.reserve((2*rcell+1)*(2*rcell+1)*(2*rcell+1));

  double active_zone_mass = 0.0;

  /* Deterministic k-j-i loop order for MPI reproducibility. */
  for (int k = klo; k <= khi; k++) {
    const double dz = (double(k - k0)) * double(cell_width);
    for (int j = jlo; j <= jhi; j++) {
      const double dy = (double(j - j0)) * double(cell_width);
      for (int i = ilo; i <= ihi; i++) {
        const double dx = (double(i - i0)) * double(cell_width);
        const double dist2 = dx*dx + dy*dy + dz*dz;
        if (dist2 > r2)
          continue;

        if (i < isx || i > iex || j < isy || j > iey || k < isz || k > iez)
          continue;

        const int nindex = (k*ny + j)*nx + i;
        const double rho = GridData->BaryonField[DensNum][nindex];
        if (rho <= 0.0)
          continue;

        const double mcell = rho * cell_volume_code;
        if (mcell <= 0.0)
          continue;

        double vx = 0.0, vy = 0.0, vz = 0.0;
        if (HydroMethod == Zeus_Hydro) {
          vx = 0.5*(GridData->BaryonField[Vel1Num][nindex] +
                    GridData->BaryonField[Vel1Num][nindex+xo]);
          vy = 0.5*(GridData->BaryonField[Vel2Num][nindex] +
                    GridData->BaryonField[Vel2Num][nindex+yo]);
          vz = 0.5*(GridData->BaryonField[Vel3Num][nindex] +
                    GridData->BaryonField[Vel3Num][nindex+zo]);
        } else {
          vx = GridData->BaryonField[Vel1Num][nindex];
          vy = GridData->BaryonField[Vel2Num][nindex];
          vz = GridData->BaryonField[Vel3Num][nindex];
        }

        active_index.push_back(nindex);
        active_mass.push_back(mcell);
        active_velx.push_back(vx);
        active_vely.push_back(vy);
        active_velz.push_back(vz);
        active_zone_mass += mcell;
      }
    }
  }

  if (active_zone_mass < bh_mass_code || active_index.empty()) {
    (*nskipped_insufficient_gas)++;
    return SUCCESS;
  }

  const float seed_density_before = GridData->BaryonField[DensNum][index0];

  std::vector<double> dm_removed(active_index.size(), 0.0);

  double removed_mass = 0.0;
  double momentum_x = 0.0;
  double momentum_y = 0.0;
  double momentum_z = 0.0;
  double remaining_to_remove = bh_mass_code;
  double remaining_available = active_zone_mass;

  for (int n = 0; n < int(active_index.size()); n++) {
    const double mcell = active_mass[n];

    double dm = 0.0;
    if (n == int(active_index.size()) - 1)
      dm = remaining_to_remove;
    else if (remaining_available > 0.0)
      dm = remaining_to_remove * (mcell / remaining_available);

    dm = max(0.0, dm);
    dm = min(dm, mcell);

    dm_removed[n] = dm;
    remaining_to_remove -= dm;
    remaining_available -= mcell;
  }

  if (remaining_to_remove != 0.0) {
    if (remaining_to_remove > 0.0) {
      for (int n = int(active_index.size()) - 1;
           n >= 0 && remaining_to_remove > 0.0; n--) {
        const double cap = active_mass[n] - dm_removed[n];
        if (cap <= 0.0)
          continue;
        const double add = min(cap, remaining_to_remove);
        dm_removed[n] += add;
        remaining_to_remove -= add;
      }
    } else {
      for (int n = int(active_index.size()) - 1;
           n >= 0 && remaining_to_remove < 0.0; n--) {
        if (dm_removed[n] <= 0.0)
          continue;
        const double sub = min(dm_removed[n], -remaining_to_remove);
        dm_removed[n] -= sub;
        remaining_to_remove += sub;
      }
    }
  }

  for (int n = 0; n < int(active_index.size()); n++) {
    const int nindex = active_index[n];
    const double dm = dm_removed[n];
    if (dm <= 0.0)
      continue;

    const double rho_old = GridData->BaryonField[DensNum][nindex];
    const double mold = rho_old * cell_volume_code;
    const double mnew = max(0.0, mold - dm);
    double rho_new = mnew / cell_volume_code;
    if (rho_new < 0.0)
      rho_new = 0.0;

    if (HydroMethod == PPM_DirectEuler) {
      if (mold > 0.0 && mnew > 0.0) {
        const double vx = active_velx[n];
        const double vy = active_vely[n];
        const double vz = active_velz[n];
        const double v2 = vx*vx + vy*vy + vz*vz;

        const double etot = mold * GridData->BaryonField[TENum][nindex];
        double eint = 0.0;
        if (GENum >= 0 && DualEnergyFormalism)
          eint = mold * GridData->BaryonField[GENum][nindex];
        else
          eint = etot - 0.5*mold*v2;

        const double ke = 0.5*mold*v2;
        const double frac = mnew / mold;
        const double eint_new = eint * frac;
        const double ke_new = ke * frac;
        const double etot_new = eint_new + ke_new;
        GridData->BaryonField[TENum][nindex] = float(etot_new / mnew);
      }
    } else if (HydroMethod == Zeus_Hydro) {
      // Total energy is internal specific energy for Zeus; unchanged.
      ;
    } else {
      ENZO_FAIL("BHSeed kernel removal does not support RK Hydro or RK MHD.");
    }

    GridData->BaryonField[DensNum][nindex] = float(rho_new);

    removed_mass += dm;
    momentum_x += dm * active_velx[n];
    momentum_y += dm * active_vely[n];
    momentum_z += dm * active_velz[n];
  }

  const double rel_mass_err =
    (bh_mass_code > 0.0) ? fabs(removed_mass - bh_mass_code) / bh_mass_code : 0.0;

  ParticleEntry p;
  p.Position[0] = cand.pos[0];
  p.Position[1] = cand.pos[1];
  p.Position[2] = cand.pos[2];
  p.Velocity[0] = (bh_mass_code > 0.0) ? float(momentum_x / bh_mass_code) : 0.0f;
  p.Velocity[1] = (bh_mass_code > 0.0) ? float(momentum_y / bh_mass_code) : 0.0f;
  p.Velocity[2] = (bh_mass_code > 0.0) ? float(momentum_z / bh_mass_code) : 0.0f;

  for (int n = 0; n < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES; n++)
    p.Attribute[n] = 0.0f;

  p.Type = PARTICLE_TYPE_MBH;
  p.Number = INT_UNDEFINED;
  p.Mass = bh_mass_code;

  if (BHSeedVerbose >= 2) {
    const double mom_x_particle = double(p.Mass) * double(p.Velocity[0]);
    const double mom_y_particle = double(p.Mass) * double(p.Velocity[1]);
    const double mom_z_particle = double(p.Mass) * double(p.Velocity[2]);
    const double rel_mom_x_err =
      (fabs(momentum_x) > 0.0) ? fabs(momentum_x - mom_x_particle) / fabs(momentum_x) : 0.0;
    const double rel_mom_y_err =
      (fabs(momentum_y) > 0.0) ? fabs(momentum_y - mom_y_particle) / fabs(momentum_y) : 0.0;
    const double rel_mom_z_err =
      (fabs(momentum_z) > 0.0) ? fabs(momentum_z - mom_z_particle) / fabs(momentum_z) : 0.0;
    FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;
    fprintf(logptr,
            "[BHSEED_DEBUG] active_zone_mass_code=%.15e sum_dm_removed=%.15e "
            "bh_mass_code=%.15e p_mass_code=%.15e rel_mass_err=%.3e "
            "mom_x_removed=%.15e mom_y_removed=%.15e mom_z_removed=%.15e "
            "mom_x_particle=%.15e mom_y_particle=%.15e mom_z_particle=%.15e "
            "rel_mom_x_err=%.3e rel_mom_y_err=%.3e rel_mom_z_err=%.3e\n",
            active_zone_mass, removed_mass, bh_mass_code, double(p.Mass), rel_mass_err,
            momentum_x, momentum_y, momentum_z,
            mom_x_particle, mom_y_particle, mom_z_particle,
            rel_mom_x_err, rel_mom_y_err, rel_mom_z_err);
  }

  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_CREATION_TIME)
    p.Attribute[PARTICLE_ATTRIBUTE_CREATION_TIME] = GridData->Time;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_DYNAMICAL_TIME)
    p.Attribute[PARTICLE_ATTRIBUTE_DYNAMICAL_TIME] = 0.0f;

  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_METALLICITY_FRACTION) {
    float metal_fraction = 0.0f;
    int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum;
    int MetalIaNum, MetalIINum;
    if (GridData->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
                                       MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum)
        != FAIL) {
      float den = seed_density_before;
      if (den > 0.0f) {
        if (MetalNum != -1 && SNColourNum != -1)
          metal_fraction = (GridData->BaryonField[MetalNum][index0] +
                            GridData->BaryonField[SNColourNum][index0]) / den;
        else if (MetalNum != -1)
          metal_fraction = GridData->BaryonField[MetalNum][index0] / den;
        else if (SNColourNum != -1)
          metal_fraction = GridData->BaryonField[SNColourNum][index0] / den;
      }
    }
    p.Attribute[PARTICLE_ATTRIBUTE_METALLICITY_FRACTION] = metal_fraction;
  }

  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_CHANNEL)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_CHANNEL] = float(cand.seed_channel);
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_REDSHIFT)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_REDSHIFT] = cand.seed_redshift;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_PATCH_MASS)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_PATCH_MASS] = cand.patch_mass;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_PATCH_METALLICITY)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_PATCH_METALLICITY] = cand.patch_metallicity;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_PATCH_DENSITY_PEAK)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_PATCH_DENSITY_PEAK] = cand.patch_density_peak;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_KERNEL_COMPLETE)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_KERNEL_COMPLETE] = float(cand.kernel_complete);
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_HOST_DM_DENSITY)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_HOST_DM_DENSITY] = cand.host_dm_density;
  if (NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHSEED_ACCEPT_RANK)
    p.Attribute[PARTICLE_ATTRIBUTE_BHSEED_ACCEPT_RANK] = float(cand.accept_rank);

  if (GridData->AddOneParticleFromList(&p, 0) == FAIL)
    return FAIL;

  BHSeedRecordCreatedSeed(cand.pos[0], cand.pos[1], cand.pos[2]);
  if (GridData->BaryonField[GridData->NumberOfBaryonFields] != NULL)
    GridData->BaryonField[GridData->NumberOfBaryonFields][index0] = 1.0f;

  if (BHSeedVerbose >= 1) {
    FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;
    fprintf(logptr,
            "[BHSEED_SEED] level=%d x=%.8g y=%.8g z=%.8g channel=%d redshift=%.8g "
            "patch_mass=%.8g patch_metal=%.8g patch_density_peak=%.8g "
            "kernel_complete=%d host_dm_density=%.8g accept_rank=%d\n",
            BHSeedState.Level,
            cand.pos[0], cand.pos[1], cand.pos[2], cand.seed_channel,
            cand.seed_redshift, cand.patch_mass, cand.patch_metallicity,
            cand.patch_density_peak, cand.kernel_complete, cand.host_dm_density,
            cand.accept_rank);
  }

  return SUCCESS;
}

static int BHSeedSyncCreatedSeedsToGlobalCache()
{
  int local_n = int(BHSeedState.NewSeedsX.size());
  Eint32 local_n32 = Eint32(local_n);
  Eint32 any32 = local_n32;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Allreduce(&local_n32, &any32, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (any32 == 0)
    return SUCCESS;


  std::vector<FLOAT> local_xyz(3*local_n);
  for (int n = 0; n < local_n; n++) {
    local_xyz[3*n + 0] = BHSeedState.NewSeedsX[n];
    local_xyz[3*n + 1] = BHSeedState.NewSeedsY[n];
    local_xyz[3*n + 2] = BHSeedState.NewSeedsZ[n];
  }

  int global_n = local_n;
  std::vector<Eint32> counts, displs;
  std::vector<Eint32> counts3, displs3;

#ifdef USE_MPI
  if (NumberOfProcessors > 1) {
    counts.assign(NumberOfProcessors, 0);
    displs.assign(NumberOfProcessors, 0);

    MPI_Allgather(&local_n32, 1, MPI_INT, &counts[0], 1, MPI_INT, MPI_COMM_WORLD);

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

    std::vector<FLOAT> global_xyz(3*global_n, 0.0f);
    MPI_Allgatherv((local_n > 0) ? &local_xyz[0] : NULL, 3*local_n,
                   FLOATDataType,
                   (global_n > 0) ? &global_xyz[0] : NULL,
                   &counts3[0], &displs3[0], FLOATDataType,
                   MPI_COMM_WORLD);

    for (int n = 0; n < global_n; n++)
      BHSeedAppendGlobalBHUnique(global_xyz[3*n + 0],
                                 global_xyz[3*n + 1],
                                 global_xyz[3*n + 2]);
  } else
#endif
  {
    for (int n = 0; n < local_n; n++)
      BHSeedAppendGlobalBHUnique(local_xyz[3*n + 0],
                                 local_xyz[3*n + 1],
                                 local_xyz[3*n + 2]);
  }

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
  BHSeedState.SeedingWallStart = ReturnWallTime();
  BHSeedState.LevelCellWidthCode = 0.0f;

  double local_dx = 1.0e99;
  BHSeedState.LocalGridIDs.reserve(NumberOfGrids);
  BHSeedState.LocalGridPtrs.reserve(NumberOfGrids);
  for (int g = 0; g < NumberOfGrids; g++) {
    grid *GridData = Grids[g]->GridData;
    if (GridData->ReturnProcessorNumber() != MyProcessorNumber)
      continue;
    local_dx = min(local_dx, double(GridData->GetCellWidth(0, 0)));
    BHSeedState.LocalGridIDs.push_back(GridData->GetGridID());
    BHSeedState.LocalGridPtrs.push_back(GridData);
  }

  double global_dx = local_dx;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Allreduce(&local_dx, &global_dx, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
  if (global_dx < 1.0e90)
    BHSeedState.LevelCellWidthCode = float(global_dx);

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

  BHSeedState.HubbleParam = h_param;
  BHSeedState.APhys = a_param;

  BHSeedState.CandidateSeparationComKpcH =
    BHSeedMinCandidateSeparation * h_param / a_param;

  float representative_excl_com = 0.0f;
  if (BHSeedExclusionMode == 1) {
    representative_excl_com = BHSeedExclusionRadius;
  } else if (BHSeedExclusionMode == 2) {
    representative_excl_com =
      BHSeedExclusionCells * BHSeedState.LevelCellWidthCode * BHSeedState.BoxKpcH;
  } else {
    representative_excl_com = BHSeedExclusionRadius * h_param / a_param;
  }
  if (representative_excl_com < 0.0f)
    representative_excl_com = 0.0f;

  BHSeedState.ExclusionRadiusComKpcH = representative_excl_com;
  BHSeedState.ExclusionRadius2 = representative_excl_com * representative_excl_com;
  BHSeedState.CellSizeCode = 0.0f;
  if (representative_excl_com > 0.0f && BHSeedState.BoxKpcH > 0.0f)
    BHSeedState.CellSizeCode = 0.5f * representative_excl_com /
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

  const int packed_bytes = BHSeedCandidatePackedBytes();
  const int expected_packed_bytes =
    3*int(sizeof(double)) + 7*int(sizeof(float)) + 9*int(sizeof(int));
  if (packed_bytes != expected_packed_bytes)
    ENZO_FAIL("BHSeedCandidate packed size mismatch.");

  const int local_count = int(BHSeedState.LocalCandidates.size());
  int global_count = local_count;

  std::vector<BHSeedCandidate> gathered;

#ifdef USE_MPI
  std::vector<Eint32> counts;
  if (NumberOfProcessors > 1) {
    counts.assign(NumberOfProcessors, 0);
    Eint32 local_count32 = Eint32(local_count);
    MPI_Allgather((void*) &local_count32, 1, MPI_INT,
                  (void*) &counts[0], 1, MPI_INT, MPI_COMM_WORLD);

    global_count = 0;
    for (int p = 0; p < NumberOfProcessors; p++)
      global_count += counts[p];

    if (global_count > 0) {
      std::vector<Eint32> byte_counts(NumberOfProcessors, 0);
      std::vector<Eint32> byte_displs(NumberOfProcessors, 0);
      int total_bytes = 0;
      for (int p = 0; p < NumberOfProcessors; p++) {
        byte_counts[p] = counts[p] * packed_bytes;
        byte_displs[p] = total_bytes;
        total_bytes += byte_counts[p];
      }

      std::vector<char> send_buffer(local_count * packed_bytes);
      for (int n = 0; n < local_count; n++)
        BHSeedCandidatePack(BHSeedState.LocalCandidates[n],
                            &send_buffer[n*packed_bytes]);

      std::vector<char> recv_buffer(total_bytes);
      MPI_Allgatherv((local_count > 0) ? &send_buffer[0] : NULL,
                     local_count * packed_bytes, MPI_PACKED,
                     (total_bytes > 0) ? &recv_buffer[0] : NULL,
                     &byte_counts[0], &byte_displs[0], MPI_PACKED,
                     MPI_COMM_WORLD);

      gathered.resize(global_count);
      int out = 0;
      for (int p = 0; p < NumberOfProcessors; p++)
        for (int n = 0; n < counts[p]; n++, out++)
          BHSeedCandidateUnpack(&recv_buffer[byte_displs[p] + n*packed_bytes],
                                gathered[out]);
    }
  } else
#endif
  {
    gathered = BHSeedState.LocalCandidates;
    global_count = local_count;
  }

  BHSeedCandidateComparator cmp;
  cmp.RankingOrder = BHSeedRankingOrder;
  cmp.UseDeterministicTiebreak = BHSeedDeterministicTiebreak;
  std::sort(gathered.begin(), gathered.end(), cmp);

  long long ncandidates_exclusion_rejected = 0;
  long long ncandidates_dedup_rejected = 0;
  int walk_stopped_at_max = 0;

  std::vector<BHSeedCandidate> accepted;
  accepted.reserve(gathered.size());

  int max_per_pass = BHSeedMaxPerPass;
  if (max_per_pass < 0)
    max_per_pass = 0;

  const float min_sep_com = max(0.0f, BHSeedState.CandidateSeparationComKpcH);
  const float min_sep2 = min_sep_com * min_sep_com;

  if (max_per_pass == 0 && global_count > 0) {
    walk_stopped_at_max = 1;
  } else {
    for (int n = 0; n < int(gathered.size()); n++) {
      BHSeedCandidate cand = gathered[n];

      const float exclusion_radius_com =
        BHSeedExclusionRadiusComKpcHForCandidate(cand.dx_local);
      if (BHSeedPositionBlockedWithRadius(cand.pos[0], cand.pos[1], cand.pos[2],
                                          exclusion_radius_com)) {
        ncandidates_exclusion_rejected++;
        continue;
      }

      int reject_dedup = FALSE;
      if (min_sep_com > 0.0f) {
        for (int a = 0; a < int(accepted.size()); a++) {
          float d2 = BHSeedPeriodicDistance2ComKpcH(cand.pos[0], cand.pos[1], cand.pos[2],
                                                    accepted[a].pos[0], accepted[a].pos[1], accepted[a].pos[2]);
          if (d2 < min_sep2) {
            reject_dedup = TRUE;
            break;
          }
        }
      }

      if (reject_dedup) {
        ncandidates_dedup_rejected++;
        continue;
      }

      cand.accept_rank = int(accepted.size()) + 1;
      accepted.push_back(cand);

      if (int(accepted.size()) >= max_per_pass) {
        walk_stopped_at_max = 1;
        break;
      }
    }
  }

  long long nseeds_skipped_insufficient_gas_local = 0;
  for (int n = 0; n < int(accepted.size()); n++)
    if (BHSeedCreateAcceptedCandidate(accepted[n],
                                      &nseeds_skipped_insufficient_gas_local) == FAIL)
      return FAIL;

  if (BHSeedSyncCreatedSeedsToGlobalCache() == FAIL)
    return FAIL;

  long long local[16], global[16];
  local[0] = BHSeedState.NCandLocal;
  for (int n = 0; n < 6; n++)
    local[n+1] = BHSeedState.DiagLocal[n];
  local[7] = BHSeedState.MassGateLocal;
  local[8] = BHSeedState.DistBlockedLocal;
  local[9] = BHSeedState.CreatedLocal;
  local[10] = BHSeedState.DiagLocal[6];
  local[11] = BHSeedState.DiagLocal[7];
  local[12] = BHSeedState.KernelEvaluatedLocal;
  local[13] = BHSeedState.KernelTruncatedLocal;
  local[14] = BHSeedState.LegacyCellMassWouldFailLocal;
  local[15] = BHSeedState.EnclosedMassGateLocal;

#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Reduce(local, global, 16, MPI_LONG_LONG_INT, MPI_SUM,
               ROOT_PROCESSOR, MPI_COMM_WORLD);
  else
#endif
    for (int i = 0; i < 16; i++)
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

  long long ncandidates_gathered_local = local_count;
  long long ncandidates_gathered_global = ncandidates_gathered_local;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Reduce(&ncandidates_gathered_local, &ncandidates_gathered_global, 1,
               MPI_LONG_LONG_INT, MPI_SUM, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  long long ncandidates_dedup_rejected_local = ncandidates_dedup_rejected;
  long long ncandidates_dedup_rejected_global = ncandidates_dedup_rejected_local;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Reduce(&ncandidates_dedup_rejected_local, &ncandidates_dedup_rejected_global,
               1, MPI_LONG_LONG_INT, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  long long ncandidates_exclusion_rejected_local = ncandidates_exclusion_rejected;
  long long ncandidates_exclusion_rejected_global = ncandidates_exclusion_rejected_local;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Reduce(&ncandidates_exclusion_rejected_local, &ncandidates_exclusion_rejected_global,
               1, MPI_LONG_LONG_INT, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  long long nseeds_skipped_insufficient_gas_global = nseeds_skipped_insufficient_gas_local;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Reduce(&nseeds_skipped_insufficient_gas_local,
               &nseeds_skipped_insufficient_gas_global,
               1, MPI_LONG_LONG_INT, MPI_SUM, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  int walk_stopped_at_max_global = walk_stopped_at_max;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Reduce(&walk_stopped_at_max, &walk_stopped_at_max_global,
               1, MPI_INT, MPI_MAX, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  int dm_available_local = BHSeedState.DMFieldAvailableLocal;
  int dm_available_global = dm_available_local;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Reduce(&dm_available_local, &dm_available_global, 1, MPI_INT, MPI_MAX,
               ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  double seeding_wall_ms_local =
    1000.0 * (ReturnWallTime() - BHSeedState.SeedingWallStart);
  double seeding_wall_ms = seeding_wall_ms_local;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Reduce(&seeding_wall_ms_local, &seeding_wall_ms, 1, MPI_DOUBLE, MPI_MAX,
               ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  long long accepted_count_local = accepted.size();
  long long accepted_count_global = accepted_count_local;
  long long created_count_global_all = BHSeedState.CreatedLocal;
  long long skipped_count_global_all = nseeds_skipped_insufficient_gas_local;
#ifdef USE_MPI
  if (NumberOfProcessors > 1) {
    MPI_Allreduce(&accepted_count_local, &accepted_count_global, 1,
                  MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&BHSeedState.CreatedLocal, &created_count_global_all, 1,
                  MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nseeds_skipped_insufficient_gas_local,
                  &skipped_count_global_all, 1,
                  MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
  }
#endif

  int invariant_bad = 0;
  if (created_count_global_all + skipped_count_global_all != accepted_count_global)
    invariant_bad = 1;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    MPI_Bcast(&invariant_bad, 1, MPI_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif
  if (invariant_bad)
    ENZO_FAIL("BHSeed invariant failed: created + skipped must equal accepted candidates.");

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    BHSeedStepCounter++;
    FILE *logptr = (Outfptr != NULL) ? Outfptr : stdout;

    if (!BHSeedDMAvailabilityLogged) {
      fprintf(logptr, "[BHSEED_INFO] dm_density_field_available=%d\n",
              dm_available_global);
      BHSeedDMAvailabilityLogged = TRUE;
    }

    if (BHSeedState.LevelCellWidthCode > 0.0f) {
      if (BHSeedPatchRadius > 3.0f * BHSeedState.LevelCellWidthCode)
        fprintf(logptr, "WARNING: BHSeedPatchRadius exceeds ghost-zone support for this grid resolution.\n");
      if (BHSeedPatchRadius / BHSeedState.LevelCellWidthCode < 1.5f)
        fprintf(logptr, "WARNING: BHSeedPatchRadius is under-resolved for this grid resolution.\n");
    }

    if (global[12] > 100)
      fprintf(logptr, "WARNING: large number of candidates reached kernel evaluation. Consider tightening gate thresholds.\n");

    if (global[12] > 0 && global[13] > (long long) (0.25 * double(global[12])))
      fprintf(logptr, "WARNING: >25%% of kernel evaluations were truncated at grid boundaries.\n");

    if (ncandidates_gathered_global > 500)
      fprintf(logptr, "WARNING: large number of candidates gathered. Consider tightening gate thresholds.\n");

    float excl_phys_kpc = BHSeedExclusionRadius;
    if (BHSeedExclusionMode == 1)
      excl_phys_kpc = BHSeedExclusionRadius * BHSeedState.APhys / BHSeedState.HubbleParam;
    else if (BHSeedExclusionMode == 2)
      excl_phys_kpc = BHSeedState.ExclusionRadiusComKpcH * BHSeedState.APhys / BHSeedState.HubbleParam;

    fprintf(logptr,
            "[BHSEED] step=%d level=%d z=%.4f a_phys=%.6f "
            "excl_phys_kpc=%.1f excl_com_kpch=%.1f "
            "cell_code=%.6e nbins=%d,%d,%d "
            "ncand_local_min=%lld ncand_local_max=%lld ncand_global=%lld "
            "ngates_density=%lld ngates_temp=%lld ngates_metal=%lld "
            "ngates_conv=%lld ngates_cool=%lld ngates_bound=%lld "
            "ngates_mass=%lld dist_blocked=%lld created=%lld total_mbh=%d pre_cache_bh=%d "
            "ngates_finestlevel=%lld ngates_peak=%lld "
            "nkernel_evaluated=%lld nkernel_truncated=%lld "
            "nlegacy_cellmass_would_fail=%lld ngates_enclosedmass=%lld "
            "seeding_wall_ms=%.3f "
            "nseeds_created=%lld ncandidates_gathered=%lld "
            "ncandidates_dedup_rejected=%lld ncandidates_exclusion_rejected=%lld "
            "walk_stopped_at_max=%d nseeds_skipped_insufficient_gas=%lld\n",
            BHSeedStepCounter, BHSeedState.Level, BHSeedState.Zred,
            double(BHSeedState.A)/(1.0 + InitialRedshift),
            excl_phys_kpc, BHSeedState.ExclusionRadiusComKpcH,
            BHSeedState.CellSizeCode,
            BHSeedState.NBins[0], BHSeedState.NBins[1], BHSeedState.NBins[2],
            ncand_min, ncand_max, global[0],
            global[1], global[2], global[3], global[4], global[5], global[6],
            global[7], global[8], global[9], int(BHSeedState.X.size()),
            BHSeedState.PreCacheBH, global[10], global[11],
            global[12], global[13], global[14], global[15], seeding_wall_ms,
            global[9], ncandidates_gathered_global,
            ncandidates_dedup_rejected_global, ncandidates_exclusion_rejected_global,
            walk_stopped_at_max_global, nseeds_skipped_insufficient_gas_global);
  }

  BHSeedState.Active = FALSE;

  return SUCCESS;
}
