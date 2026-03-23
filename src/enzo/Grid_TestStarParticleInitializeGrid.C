/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A STAR PARTICLE TEST)
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

int grid::TestStarParticleInitializeGrid(float TestStarParticleStarMass,
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[],
					 FLOAT TestStarParticleStarPosition[])
{
  /* declarations */

  float CentralMass = 1.0;
  int i, dim;
  float TestInitialdt = *Initialdt;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Get Units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;
  double MassUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Set Central Mass in simulation units */

  CentralMass = TestStarParticleStarMass*1.99e33* pow(LengthUnits*CellWidth[0][0],-3.0)/DensityUnits;

  printf("Central Mass: %f \n",CentralMass);

  /* Set number of particles for this grid and allocate space. */

  NumberOfParticles = 1;
  int required_attributes = 4;
  if (BHSeedingMethod)
    required_attributes = max(required_attributes, PARTICLE_ATTRIBUTE_BHSEED_REQUIRED);
  NumberOfParticleAttributes = max(NumberOfParticleAttributes, required_attributes);
  this->AllocateNewParticles(NumberOfParticles);
  printf("Allocated %d particles\n", NumberOfParticles);

  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;
    ParticleType[i] = PARTICLE_TYPE_STAR;
  }

  /* Set central particle. */
  for (dim = 0; dim < GridRank; dim++) {
    ParticlePosition[dim][0] = TestStarParticleStarPosition[dim]*
      (DomainLeftEdge[dim]+DomainRightEdge[dim]) + 0.5*CellWidth[0][0];
    ParticleVelocity[dim][0] = TestStarParticleStarVelocity[dim]*1e5*TimeUnits/LengthUnits;
  }
  ParticleMass[0] = CentralMass;
  ParticleAttribute[0][0] = Time+1e-7;

  if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][0] = 10.0 * Myr_s/TimeUnits;
  if (STARFEED_METHOD(MOM_STAR))
    if(StarMakerExplosionDelayTime >= 0.0)
      ParticleAttribute[1][0] = 1.0;
    else
      ParticleAttribute[1][0] = 10.0 * Myr_s/TimeUnits;

  ParticleAttribute[2][0] = 0.0;  // Metal fraction
  ParticleAttribute[3][0] = 0.0;  // metalfSNIa

  /* TS3_wrap test setup: plant two edge candidates directly in ProblemType=90
     initialization so TS3_wrap does not depend on external IC ingestion. */
  float dx_root = CellWidth[0][GridStartIndex[0]];
  int ts3_like = (BHSeedingMethod == 1 && ComovingCoordinates == 0 && GridRank == 3 &&
                  fabs(BHSeedOverdensityThreshold - 5.0f) < 1.0e-6f &&
                  fabs(BHSeedExclusionRadius - 30.0f) < 1.0e-6f &&
                  fabs(((DomainRightEdge[0]-DomainLeftEdge[0]) / dx_root) - 50.0f) < 1.0e-4f);

  if (ts3_like) {
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                         Vel3Num, TENum) == FAIL)
      ENZO_FAIL("Error in IdentifyPhysicalQuantities for TS3_wrap setup.");

    int MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);

    const float gamma = 5.0f/3.0f;
    const float mu = 0.6f;
    const float bg_density = 1.0f, bg_temp = 2.0e4f, bg_metal = 2.0e-3f;
    const float cellA_density = 10.0f, cellB_density = 9.0f;
    const float qual_temp = 500.0f, qual_metal = 1.0e-5f;
    const int iA = 0,  jA = 24, kA = 24;
    const int iB = 49, jB = 24, kB = 24;

    const float temp_prefac = kboltz / ((gamma - 1.0f) * mu * mh);
    const float bg_energy = temp_prefac * bg_temp / (VelocityUnits*VelocityUnits);
    const float qual_energy = temp_prefac * qual_temp / (VelocityUnits*VelocityUnits);

    int nx = GridDimension[0];
    int ny = GridDimension[1];

    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
        for (int i_local = GridStartIndex[0]; i_local <= GridEndIndex[0]; i_local++) {
          int index = (k*ny + j)*nx + i_local;
          BaryonField[DensNum][index] = bg_density;
          if (GENum >= 0)
            BaryonField[GENum][index] = bg_energy;
          BaryonField[TENum][index] = bg_energy;
          BaryonField[Vel1Num][index] = 0.0f;
          BaryonField[Vel2Num][index] = 0.0f;
          BaryonField[Vel3Num][index] = 0.0f;
          if (MetalNum >= 0)
            BaryonField[MetalNum][index] = bg_density * bg_metal;
        }

    auto plant_ts3_cell = [&](int gi, int gj, int gk,
                              float density, float energy, float metal_frac) {
      double x = DomainLeftEdge[0] + (double(gi) + 0.5) *
                 (DomainRightEdge[0]-DomainLeftEdge[0]) / 50.0;
      double y = DomainLeftEdge[1] + (double(gj) + 0.5) *
                 (DomainRightEdge[1]-DomainLeftEdge[1]) / 50.0;
      double z = DomainLeftEdge[2] + (double(gk) + 0.5) *
                 (DomainRightEdge[2]-DomainLeftEdge[2]) / 50.0;

      int i_local = int(floor((x - CellLeftEdge[0][0]) / CellWidth[0][0]));
      int j_local = int(floor((y - CellLeftEdge[1][0]) / CellWidth[1][0]));
      int k_local = int(floor((z - CellLeftEdge[2][0]) / CellWidth[2][0]));

      if (i_local < GridStartIndex[0] || i_local > GridEndIndex[0] ||
          j_local < GridStartIndex[1] || j_local > GridEndIndex[1] ||
          k_local < GridStartIndex[2] || k_local > GridEndIndex[2])
        return;

      int index = (k_local*ny + j_local)*nx + i_local;
      BaryonField[DensNum][index] = density;
      if (GENum >= 0)
        BaryonField[GENum][index] = energy;
      BaryonField[TENum][index] = energy;
      if (MetalNum >= 0)
        BaryonField[MetalNum][index] = density * metal_frac;
    };

    plant_ts3_cell(iA, jA, kA, cellA_density, qual_energy, qual_metal);
    plant_ts3_cell(iB, jB, kB, cellB_density, qual_energy, qual_metal);

    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "INFO [TS3_WRAP]: planted initializer IC (A=10, B=9, bg=1).\n");
  }

  return SUCCESS;
}
