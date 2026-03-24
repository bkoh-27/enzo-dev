/***********************************************************************
/
/  BH SEED CANDIDATE KERNEL (PHASE 2 GATE ORDER)
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

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"

int star_maker_bh_seed(int *nx, int *ny, int *nz, int *ibuff, int *imethod,
                       float *d, float *dm, float *temp, float *u, float *v, float *w,
                       float *cooltime, float *r, float *metal,
                       float *dx, float *d1, float *t1,
                       float *odthresh, float *metalthresh, float *tempthresh,
                       int *veldivcrit, int *thermalcrit, int *selfboundcrit,
                       int *requirefinestlevel, int *requirelocalpeak,
                       int *ncand, int *maxcand, int *cand_index, float *cand_density,
                       int *diag)
{
  *ncand = 0;
  for (int n = 0; n < 8; n++)
    diag[n] = 0;

  const int xo = 1;
  const int yo = *nx;
  const int zo = (*nx) * (*ny);
  const int max_output = (*maxcand > 0) ? *maxcand : 0;

  int count = 0;

  for (int k = *ibuff; k < *nz - *ibuff; k++) {
    for (int j = *ibuff; j < *ny - *ibuff; j++) {
      int index = (k * (*ny) + j) * (*nx) + *ibuff;
      for (int i = *ibuff; i < *nx - *ibuff; i++, index++) {

        if (*requirefinestlevel == 1 && r[index] != 0.0f) {
          diag[6]++;
          continue;
        }

        if (!isfinite(d[index])) {
          diag[0]++;
          continue;
        }

        if (d[index] < *odthresh) {
          diag[0]++;
          continue;
        }

        /* Phase 2 gate order: local-peak check is immediately after density gate. */
        if (*requirelocalpeak == 1) {
          int is_local_peak = TRUE;
          for (int kk = -1; kk <= 1 && is_local_peak; kk++)
            for (int jj = -1; jj <= 1 && is_local_peak; jj++)
              for (int ii = -1; ii <= 1; ii++) {
                if (ii == 0 && jj == 0 && kk == 0)
                  continue;

                int ni = i + ii;
                int nj = j + jj;
                int nk = k + kk;

                /* Out-of-bounds neighbors are treated as rho=0 by skipping compare. */
                if (ni < 0 || ni >= *nx ||
                    nj < 0 || nj >= *ny ||
                    nk < 0 || nk >= *nz)
                  continue;

                int nindex = (nk * (*ny) + nj) * (*nx) + ni;
                if (d[index] < d[nindex]) {
                  is_local_peak = FALSE;
                  break;
                }
              }

          if (!is_local_peak) {
            diag[7]++;
            continue;
          }
        }

        if (temp[index] > *tempthresh) {
          diag[1]++;
          continue;
        }

        if (metal != NULL && metal[index] > *metalthresh) {
          diag[2]++;
          continue;
        }

        const int has_forward = (i+1 < *nx && j+1 < *ny && k+1 < *nz);
        const int has_centered = (i-1 >= 0 && i+1 < *nx &&
                                  j-1 >= 0 && j+1 < *ny &&
                                  k-1 >= 0 && k+1 < *nz);

        float div = 0.0f;
        if (*veldivcrit == 1) {
          if ((*imethod == Zeus_Hydro && !has_forward) ||
              (*imethod != Zeus_Hydro && !has_centered)) {
            diag[3]++;
            continue;
          }

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
          if (!has_centered) {
            diag[5]++;
            continue;
          }

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

        if (count < max_output) {
          cand_index[count] = index;
          cand_density[count] = d[index];
          count++;
        }
      }
    }
  }

  *ncand = count;

  return SUCCESS;
}
