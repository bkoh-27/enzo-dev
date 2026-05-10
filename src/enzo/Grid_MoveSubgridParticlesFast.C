/***********************************************************************
/
/  GRID CLASS (MOVE APPROPRIATE PARTICLES FROM SPECIFIED GRID TO THIS GRID)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
//
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
 
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

static int ParticleHasAbsoluteBHMass(int particle_type)
{
  return particle_type == PARTICLE_TYPE_MBH ||
         particle_type == PARTICLE_TYPE_BLACK_HOLE;
}

static void ReportBHParticleAMRMassTransfer(const char *path,
                                            const char *direction,
                                            int from_level, int to_level,
                                            const char *factor_name,
                                            double factor,
                                            int particle_type,
                                            PINT particle_id,
                                            double mass_before,
                                            double mass_after,
                                            float **attributes,
                                            int particle_index,
                                            const char *action)
{
  if (!ParticleHasAbsoluteBHMass(particle_type))
    return;

  const int has_mass_attrs =
    attributes != NULL &&
    NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BH_FORMATION_MASS &&
    NumberOfParticleAttributes > PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS;

  double formation_mass = 0.0;
  double accreted_mass = 0.0;
  double expected_mass = 0.0;
  double ratio_before = -1.0;
  double ratio_after = -1.0;
  int suspicious = FALSE;

  if (has_mass_attrs) {
    formation_mass =
      attributes[PARTICLE_ATTRIBUTE_BH_FORMATION_MASS][particle_index];
    accreted_mass =
      attributes[PARTICLE_ATTRIBUTE_BHACCR_ACCRETED_MASS][particle_index];
    expected_mass = formation_mass + accreted_mass;

    if (expected_mass > 0.0 && mass_before > 0.0)
      ratio_before = mass_before / expected_mass;
    if (expected_mass > 0.0 && mass_after > 0.0)
      ratio_after = mass_after / expected_mass;

    suspicious =
      ratio_before > 1.5 || ratio_after > 1.5 ||
      (ratio_before > 0.0 && ratio_before < 1.0/1.5) ||
      (ratio_after > 0.0 && ratio_after < 1.0/1.5);
  }

  if (!suspicious && BHAccretionVerbose < 2)
    return;

  FILE *log_fptr = suspicious ? stderr : (Outfptr != NULL ? Outfptr : stdout);
  fprintf(log_fptr,
          "[BH_PARTICLE_AMR_MASS_TRANSFER] path=%s direction=%s "
          "from_level=%d to_level=%d factor_name=%s factor=%.8e "
          "particle_type=%d particle_id=%lld ParticleMass_before=%.16e "
          "ParticleMass_after=%.16e BHFormationMass=%.16e "
          "BHAccretedMass=%.16e expected_mass=%.16e "
          "ratio_before_to_expected=%.8e ratio_after_to_expected=%.8e "
          "processor=%d action=%s\n",
          path, direction, from_level, to_level, factor_name, factor,
          particle_type, (long long) particle_id, mass_before, mass_after,
          formation_mass, accreted_mass, expected_mass, ratio_before,
          ratio_after, MyProcessorNumber, action);
  fflush(log_fptr);
}

int grid::MoveSubgridParticlesFast(int NumberOfSubgrids, grid* ToGrids[],
				   int AllLocal)
{

  if (debug1) 
    printf("MoveSubgridParticlesFast: %"ISYM"\n", NumberOfParticles);
 
  /* If there are no particles to move, we're done. */
 
  if (NumberOfParticles == 0 || NumberOfSubgrids == 0)
    return SUCCESS;
 
  int i, j, dim, index, subgrid, n;
 
  /* Initialize. */
 
  int *ParticlesToMove = new int[NumberOfSubgrids];
  for (i = 0; i < NumberOfSubgrids; i++)
    ParticlesToMove[i] = 0;
 
  /* Error check. */
 
  if (BaryonField[NumberOfBaryonFields] == NULL &&
      MyProcessorNumber == ProcessorNumber) {
    ENZO_FAIL("Subgrid field not present.\n");
  }
 
  /* Loop over particles and count the number in each subgrid. */
 
  int i0 = 0, j0 = 0, k0 = 0;
  if (MyProcessorNumber == ProcessorNumber) {
    for (i = 0; i < NumberOfParticles; i++) {
 
      /* Compute index of particle position. */
 
      i0 = int((ParticlePosition[0][i] - CellLeftEdge[0][0])/CellWidth[0][0]);
      if (GridRank > 0)
       j0 = int((ParticlePosition[1][i] - CellLeftEdge[1][0])/CellWidth[1][0]);
      if (GridRank > 1)
       k0 = int((ParticlePosition[2][i] - CellLeftEdge[2][0])/CellWidth[2][0]);
 
      i0 = max(min(GridEndIndex[0], i0), GridStartIndex[0]);
      j0 = max(min(GridEndIndex[1], j0), GridStartIndex[1]);
      k0 = max(min(GridEndIndex[2], k0), GridStartIndex[2]);
 
      index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;
 
      /* Find subgrid number of this particle, and add to count. */
 
      subgrid = nint(BaryonField[NumberOfBaryonFields][index])-1;
      if (subgrid >= 0)
	ParticlesToMove[subgrid]++;
      if (subgrid < -1 || subgrid > NumberOfSubgrids-1) {
	ENZO_VFAIL("particle subgrid (%"ISYM"/%"ISYM") out of range\n", subgrid,
		NumberOfSubgrids)
      }
 
    }  // end: loop over particles
  } // end: if (MyProcessorNumber)
 
  /* Communicate number of send particles to subgrids */
 
  if (AllLocal == FALSE)
    for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
      if (CommunicationBroadcastValue(&ParticlesToMove[subgrid],
				      ProcessorNumber) == FAIL) {
	ENZO_FAIL("Error in CommunicationBroadcastValue.\n");
      }
/*
    if ((MyProcessorNumber == ProcessorNumber ||
	 MyProcessorNumber == ToGrids[subgrid]->ProcessorNumber) &&
	ProcessorNumber != ToGrids[subgrid]->ProcessorNumber) {
      ENZO_FAIL("this routine not parallelized.\n");
      if (CommunicationSendInt(MyProcessorNumber,
			       ToGrids[subgrid]->ProcessorNumber,
			       &ParticlesToMove[subgrid]) == FAIL) {
        ENZO_FAIL("Error in CommunicationSendInt.\n");
      }
    }
*/
 
  /* Allocate space on all the subgrids with particles. */
 
  if (MyProcessorNumber == ProcessorNumber)
    for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
 
      if (ParticlesToMove[subgrid] > 0) {
 
	if (ToGrids[subgrid]->ParticlePosition[0] != NULL ||
	    ToGrids[subgrid]->NumberOfParticles != 0) {
	  ENZO_VFAIL("Particles already in subgrid %"ISYM" (n=%"ISYM", nm=%"ISYM")\n",
		  subgrid, ToGrids[subgrid]->NumberOfParticles,
		  ParticlesToMove[subgrid])
	}
 
	ToGrids[subgrid]->AllocateNewParticles(ParticlesToMove[subgrid]);
 
	if (debug1) printf("MoveSubgridParticles: subgrid[%"ISYM"] = %"ISYM"\n",
			  subgrid, ParticlesToMove[subgrid]);
 
      } // end: if (ParticlesToMove > 0)
 
  /* Compute the increase in mass for particles moving to the subgrid. */
 
  float RefinementFactors[MAX_DIMENSION], MassIncrease = 1.0;
  this->ComputeRefinementFactorsFloat(ToGrids[0], RefinementFactors);
  for (dim = 0; dim < GridRank; dim++)
    MassIncrease *= RefinementFactors[dim];
 
  if (MyProcessorNumber == ProcessorNumber) {
 
    /* Loop over particles and move them to the appropriate ToGrid, depending
       on their position. */
 
    for (i = 0; i < NumberOfParticles; i++) {
 
      /* Compute index of particle position. */
 
      i0 = int((ParticlePosition[0][i] - CellLeftEdge[0][0])/CellWidth[0][0]);
      if (GridRank > 0)
       j0 = int((ParticlePosition[1][i] - CellLeftEdge[1][0])/CellWidth[1][0]);
      if (GridRank > 1)
       k0 = int((ParticlePosition[2][i] - CellLeftEdge[2][0])/CellWidth[2][0]);
 
      i0 = max(min(GridEndIndex[0], i0), GridStartIndex[0]);
      j0 = max(min(GridEndIndex[1], j0), GridStartIndex[1]);
      k0 = max(min(GridEndIndex[2], k0), GridStartIndex[2]);
 
      index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;
 
      /* Find subgrid number of this particle, and move it. */
 
      subgrid = nint(BaryonField[NumberOfBaryonFields][index])-1;
 
      if (subgrid >= 0) {
	n = ToGrids[subgrid]->NumberOfParticles;
	const int is_bh_mass_absolute =
	  ParticleHasAbsoluteBHMass(ParticleType[i]);
	const float mass_before = ParticleMass[i];
	const float mass_after = is_bh_mass_absolute ?
	  mass_before : mass_before * MassIncrease;
	ToGrids[subgrid]->ParticleMass[n] = mass_after;
	ToGrids[subgrid]->ParticleNumber[n] = ParticleNumber[i];
	ToGrids[subgrid]->ParticleType[n] = ParticleType[i];
	for (dim = 0; dim < GridRank; dim++) {
	 ToGrids[subgrid]->ParticlePosition[dim][n] = ParticlePosition[dim][i];
	 ToGrids[subgrid]->ParticleVelocity[dim][n] = ParticleVelocity[dim][i];
	}
	for (j = 0; j < NumberOfParticleAttributes; j++)
	  ToGrids[subgrid]->ParticleAttribute[j][n] = ParticleAttribute[j][i];
	ReportBHParticleAMRMassTransfer("Grid_MoveSubgridParticlesFast",
					 "parent_to_child",
					 this->GetLevel(),
					 ToGrids[subgrid]->GetLevel(),
					 "MassIncrease", MassIncrease,
					 ParticleType[i], ParticleNumber[i],
					 mass_before, mass_after,
					 ParticleAttribute, i,
					 "bh_preserved_absolute_mass");
 
	ToGrids[subgrid]->NumberOfParticles++;
 
	/* Mark this particle removed. */
 
	ParticleMass[i] = FLOAT_UNDEFINED;
 
      } // end: if (subgrid >= 0)
 
    } // end: loop over particles
 
    /* Clean up the moved particles. */
 
    this->CleanUpMovedParticles();
 
    delete [] BaryonField[NumberOfBaryonFields];
    BaryonField[NumberOfBaryonFields] = NULL;
 
  } // end: if (MyProcessorNumber)
 
  /* Transfer particles from fake to real grids (and clean up). */
 
  for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
    if ((MyProcessorNumber == ProcessorNumber ||
         MyProcessorNumber == ToGrids[subgrid]->ProcessorNumber) &&
	ProcessorNumber != ToGrids[subgrid]->ProcessorNumber)
      if (ParticlesToMove[subgrid] != 0) {
	if (this->CommunicationSendParticles(ToGrids[subgrid],
             ToGrids[subgrid]->ProcessorNumber, 0, ParticlesToMove[subgrid], 0)
	    == FAIL) {
	  ENZO_FAIL("Error in grid->CommunicationSendParticles.\n");
	}
	if (MyProcessorNumber == ProcessorNumber)

	  ToGrids[subgrid]->DeleteAllFields();
      }
 
  delete [] ParticlesToMove;
 
  return SUCCESS;
}
