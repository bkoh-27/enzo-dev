/***********************************************************************
/
/  BH SEED CANDIDATE KERNEL (GATES 1-7)
/
/  PURPOSE:
/    Select candidate cells for BH seeding and collect gate diagnostics.
/
/  NOTES:
/    This mirrors star_maker2-style gating, with temperature and metallicity
/    ceilings, and no Jeans-mass gate.
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"

struct BHSeedCandidate {
  int index;
  float density;
};

static bool CompareBHSeedCandidates(const BHSeedCandidate &ca,
                                    const BHSeedCandidate &cb)
{
  if (ca.density != cb.density)
    return ca.density > cb.density;
  return ca.index < cb.index;
}

int star_maker_bh_seed(int *nx, int *ny, int *nz, int *ibuff, int *imethod,
                       float *d, float *dm, float *temp, float *u, float *v, float *w,
                       float *cooltime, float *r, float *metal,
                       float *dx, float *d1, float *t1,
                       float *odthresh, float *metalthresh, float *tempthresh,
                       int *veldivcrit, int *thermalcrit, int *selfboundcrit,
                       int *ncand, int *cand_index, float *cand_density,
                       int *diag)
{
  *ncand = 0;
  for (int n = 0; n < 6; n++)
    diag[n] = 0;

  const int xo = 1;
  const int yo = *nx;
  const int zo = (*nx) * (*ny);
  const int max_cells = (*nx) * (*ny) * (*nz);
  /* Not OpenMP-safe; safe because star formation/BH seeding grid loop is serial. */
  static std::vector<BHSeedCandidate> tmp_scratch;
  if (int(tmp_scratch.size()) < max_cells)
    tmp_scratch.resize(max_cells);
  BHSeedCandidate *tmp = &tmp_scratch[0];

  int count = 0;

  for (int k = *ibuff; k < *nz - *ibuff; k++) {
    for (int j = *ibuff; j < *ny - *ibuff; j++) {
      int index = (k * (*ny) + j) * (*nx) + *ibuff;
      for (int i = *ibuff; i < *nx - *ibuff; i++, index++) {

        if (r[index] != 0.0f)
          continue;

        if (d[index] < *odthresh) {
          diag[0]++;
          continue;
        }

        if (temp[index] > *tempthresh) {
          diag[1]++;
          continue;
        }

        if (metal != NULL && metal[index] > *metalthresh) {
          diag[2]++;
          continue;
        }

        float div = 0.0f;
        if (*veldivcrit == 1) {
          if (*imethod == Zeus_Hydro) {
            div = u[index+xo] - u[index] +
                  v[index+yo] - v[index] +
                  w[index+zo] - w[index];
          } else {
            div = u[index+xo] - u[index-xo] +
                  v[index+yo] - v[index-yo] +
                  w[index+zo] - w[index-zo];
          }
          if (div >= 0.0f) {
            diag[3]++;
            continue;
          }
        }

        const double dtot = (d[index] + dm[index]) * (*d1);
        if (dtot <= 0.0)
          continue;

        const double tdyn = sqrt(3.0*pi/(32.0*GravConst*dtot))/(*t1);

        if (*thermalcrit == 1 && tdyn < cooltime[index]) {
          diag[4]++;
          continue;
        }

        if (*selfboundcrit == 1) {
          const float dvx = (w[index+yo] - w[index-yo]) -
                            (v[index+zo] - v[index-zo]);
          const float dvy = (u[index+zo] - u[index-zo]) -
                            (w[index+xo] - w[index-xo]);
          const float dvz = (v[index+xo] - v[index-xo]) -
                            (u[index+yo] - u[index-yo]);
          const float divvel2 = div*div / ((*dx) * (*dx));
          const float curlvel2 = (dvx*dvx + dvy*dvy + dvz*dvz) /
                                 ((*dx) * (*dx));
          const float alpha = 0.5f * (divvel2 + curlvel2) /
                              (GravConst * d[index]);
          if (alpha >= 1.0f) {
            diag[5]++;
            continue;
          }
        }

        tmp[count].index = index;
        tmp[count].density = d[index];
        count++;
      }
    }
  }

  std::sort(tmp, tmp + count, CompareBHSeedCandidates);

  for (int n = 0; n < count; n++) {
    cand_index[n] = tmp[n].index;
    cand_density[n] = tmp[n].density;
  }

  *ncand = count;

  return SUCCESS;
}
