#ifndef BH_SEED_CANDIDATE_H
#define BH_SEED_CANDIDATE_H

#include <stddef.h>
#include <string.h>

struct BHSeedCandidate {
  double pos[3];

  float patch_mass;
  float patch_density_peak;
  float patch_metallicity;

  int kernel_complete;
  float host_dm_density;
  float seed_redshift;
  int seed_channel;

  float seed_mass_code_density;

  int owning_rank;
  int owning_grid_id;
  int cell_i;
  int cell_j;
  int cell_k;
  int cell_index;

  float dx_local;

  int accept_rank;
};

inline int BHSeedCandidatePackedBytes()
{
  return 3*int(sizeof(double)) +
         7*int(sizeof(float)) +
         9*int(sizeof(int));
}

inline void BHSeedCandidatePack(const BHSeedCandidate &cand, char *buffer)
{
  size_t offset = 0;

#define BHSEED_PACK_SCALAR(value) \
  do { \
    memcpy(buffer + offset, &(value), sizeof(value)); \
    offset += sizeof(value); \
  } while (0)

  BHSEED_PACK_SCALAR(cand.pos[0]);
  BHSEED_PACK_SCALAR(cand.pos[1]);
  BHSEED_PACK_SCALAR(cand.pos[2]);

  BHSEED_PACK_SCALAR(cand.patch_mass);
  BHSEED_PACK_SCALAR(cand.patch_density_peak);
  BHSEED_PACK_SCALAR(cand.patch_metallicity);

  BHSEED_PACK_SCALAR(cand.kernel_complete);
  BHSEED_PACK_SCALAR(cand.host_dm_density);
  BHSEED_PACK_SCALAR(cand.seed_redshift);
  BHSEED_PACK_SCALAR(cand.seed_channel);

  BHSEED_PACK_SCALAR(cand.seed_mass_code_density);

  BHSEED_PACK_SCALAR(cand.owning_rank);
  BHSEED_PACK_SCALAR(cand.owning_grid_id);
  BHSEED_PACK_SCALAR(cand.cell_i);
  BHSEED_PACK_SCALAR(cand.cell_j);
  BHSEED_PACK_SCALAR(cand.cell_k);
  BHSEED_PACK_SCALAR(cand.cell_index);

  BHSEED_PACK_SCALAR(cand.dx_local);

  BHSEED_PACK_SCALAR(cand.accept_rank);

#undef BHSEED_PACK_SCALAR
}

inline void BHSeedCandidateUnpack(const char *buffer, BHSeedCandidate &cand)
{
  size_t offset = 0;

#define BHSEED_UNPACK_SCALAR(value) \
  do { \
    memcpy(&(value), buffer + offset, sizeof(value)); \
    offset += sizeof(value); \
  } while (0)

  BHSEED_UNPACK_SCALAR(cand.pos[0]);
  BHSEED_UNPACK_SCALAR(cand.pos[1]);
  BHSEED_UNPACK_SCALAR(cand.pos[2]);

  BHSEED_UNPACK_SCALAR(cand.patch_mass);
  BHSEED_UNPACK_SCALAR(cand.patch_density_peak);
  BHSEED_UNPACK_SCALAR(cand.patch_metallicity);

  BHSEED_UNPACK_SCALAR(cand.kernel_complete);
  BHSEED_UNPACK_SCALAR(cand.host_dm_density);
  BHSEED_UNPACK_SCALAR(cand.seed_redshift);
  BHSEED_UNPACK_SCALAR(cand.seed_channel);

  BHSEED_UNPACK_SCALAR(cand.seed_mass_code_density);

  BHSEED_UNPACK_SCALAR(cand.owning_rank);
  BHSEED_UNPACK_SCALAR(cand.owning_grid_id);
  BHSEED_UNPACK_SCALAR(cand.cell_i);
  BHSEED_UNPACK_SCALAR(cand.cell_j);
  BHSEED_UNPACK_SCALAR(cand.cell_k);
  BHSEED_UNPACK_SCALAR(cand.cell_index);

  BHSEED_UNPACK_SCALAR(cand.dx_local);

  BHSEED_UNPACK_SCALAR(cand.accept_rank);

#undef BHSEED_UNPACK_SCALAR
}

#endif
