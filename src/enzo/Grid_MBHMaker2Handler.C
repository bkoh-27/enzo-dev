/***********************************************************************
/
/  GRID CLASS (RUN MBH MAKER 2 AS STANDALONE BH SEEDING MODULE)
/
************************************************************************/

#include <stdio.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "phys_constants.h"
#include "Star.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                       int *sdim1, int *sdim2, int *sdim3,
                                       int *ddim1, int *ddim2, int *ddim3,
                                       int *sstart1, int *sstart2, int *sstart3,
                                       int *dstart1, int *dstart2, int *dststart3);

int mbh_maker2(grid *ThisGrid,
               int level,
               int DensNum, int Vel1Num, int Vel2Num, int Vel3Num,
               float *temperature, float *cooling_time, float *dmfield,
               float *metal_fraction, int metal_field_present,
               int dm_field_present,
               float CellWidthTemp,
               float DensityUnits, float LengthUnits, float TimeUnits,
               int MaximumRefinementLevel);

int grid::MBHMaker2Handler(HierarchyEntry* SubgridPointer, int level,
                           float dtLevelAbove)
{
  (void) dtLevelAbove;

  if (!BHSeedingMethod)
    return SUCCESS;

  if (!BHSeedRunEveryTimestep &&
      StarFormationOncePerRootGridTimeStep &&
      !this->MakeStars)
    return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
  for (HierarchyEntry *sg = SubgridPointer; sg != NULL; sg = sg->NextGridThisLevel)
    this->ZeroSolutionUnderSubgrid(sg->GridData, ZERO_UNDER_SUBGRID_FIELD);

  int dim, i;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL)
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  float *temperature = new float[size];
  this->ComputeTemperatureField(temperature);

  float *cooling_time = new float[size];
  this->ComputeCoolingTime(cooling_time);

  float *dmfield = new float[size];
  int StartIndex[MAX_DIMENSION], Zero[] = {0, 0, 0};
  int dm_field_present = (level <= MaximumGravityRefinementLevel &&
                          GravitatingMassFieldParticles != NULL);
  if (dm_field_present) {
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      StartIndex[dim] =
        nint((CellLeftEdge[dim][0] - GravitatingMassFieldParticlesLeftEdge[dim])/
             GravitatingMassFieldParticlesCellSize);

    FORTRAN_NAME(copy3d)(GravitatingMassFieldParticles, dmfield,
                         GravitatingMassFieldParticlesDimension,
                         GravitatingMassFieldParticlesDimension+1,
                         GravitatingMassFieldParticlesDimension+2,
                         GridDimension, GridDimension+1, GridDimension+2,
                         Zero, Zero+1, Zero+2,
                         StartIndex, StartIndex+1, StartIndex+2);
  } else {
    for (i = 0; i < size; i++)
      dmfield[i] = 0.0f;
  }

  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum;
  int MetalIaNum, MetalIINum;
  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
                                 MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  int metal_field_present = (MetalNum != -1 || SNColourNum != -1);
  float *metal_fraction = NULL;

  if (metal_field_present) {
    metal_fraction = new float[size];

    if (MetalNum != -1 && SNColourNum != -1) {
      for (i = 0; i < size; i++) {
        float den = BaryonField[DensNum][i];
        if (den > 0.0f)
          metal_fraction[i] = (BaryonField[MetalNum][i] + BaryonField[SNColourNum][i]) / den;
        else
          metal_fraction[i] = 0.0f;
      }
    } else {
      int MetalField = (MetalNum != -1) ? MetalNum : SNColourNum;
      for (i = 0; i < size; i++) {
        float den = BaryonField[DensNum][i];
        if (den > 0.0f)
          metal_fraction[i] = BaryonField[MetalField][i] / den;
        else
          metal_fraction[i] = 0.0f;
      }
    }
  }

  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1;
  float TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL)
    ENZO_FAIL("Error in GetUnits.");

  float CellWidthTemp = float(CellWidth[0][0]);

  int result = mbh_maker2(this, level,
                          DensNum, Vel1Num, Vel2Num, Vel3Num,
                          temperature, cooling_time, dmfield,
                          metal_fraction, metal_field_present,
                          dm_field_present,
                          CellWidthTemp,
                          DensityUnits, LengthUnits, TimeUnits,
                          MaximumRefinementLevel);

  delete [] temperature;
  delete [] cooling_time;
  delete [] dmfield;
  if (metal_fraction != NULL)
    delete [] metal_fraction;

  if (result == FAIL)
    return FAIL;

  return SUCCESS;
}
