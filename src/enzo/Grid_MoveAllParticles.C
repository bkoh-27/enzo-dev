/***********************************************************************
/
/  GRID CLASS (MOVE ALL PARTICLES FROM SPECIFIED GRID TO THIS GRID)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  May, 2009 by John Wise
/                Move particles to "empty" grid on the processor local
/                to the subgrid.  This distributes memory usage during
/                RebuildHierarchy in simulations with nested grids.
/
/  PURPOSE:
/
/    NOTE: We assume all the from grids are at the same level!
/
************************************************************************/
 
//
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "ActiveParticle.h"

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

int grid::MoveAllParticles(int NumberOfGrids, grid* FromGrid[])
{

  if (NumberOfGrids < 1) {
    ENZO_VFAIL("NumberOfGrids(%"ISYM") must be > 0.\n", NumberOfGrids)
  }
 
  /* Determine total number of local particles. */

  int NumberOfSubgridParticles = 0;
  int NumberOfSubgridActiveParticles = 0;
  int TotalNumberOfParticles = NumberOfParticles;
  int TotalNumberOfActiveParticles = NumberOfActiveParticles;
  int i, j, grid, dim, *Type;
  PINT *Number;
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (MyProcessorNumber == FromGrid[grid]->ProcessorNumber) {
      NumberOfSubgridParticles += FromGrid[grid]->NumberOfParticles;
      NumberOfSubgridActiveParticles += FromGrid[grid]->NumberOfActiveParticles;
    }
  if (NumberOfSubgridParticles + NumberOfSubgridActiveParticles == 0) 
    return SUCCESS;
  
  TotalNumberOfParticles += NumberOfSubgridParticles;
  TotalNumberOfActiveParticles += NumberOfSubgridActiveParticles;
 
  /* Debugging info. */

  if (debug1)
    printf("MoveAllParticles: %"ISYM",%"ISYM
           " (before: ThisGrid = %"ISYM",%"ISYM").\n",
           TotalNumberOfParticles, TotalNumberOfActiveParticles,
           NumberOfParticles, NumberOfActiveParticles);
 
  /* Allocate space for the particles. */
 
  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION], *Mass,
        *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
 
  Mass = new float[TotalNumberOfParticles];
  Number = new PINT[TotalNumberOfParticles];
  Type = new int[TotalNumberOfParticles];
  for (int dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[TotalNumberOfParticles];
    Velocity[dim] = new float[TotalNumberOfParticles];
  }
  for (int i = 0; i < NumberOfParticleAttributes; i++)
    Attribute[i] = new float[TotalNumberOfParticles];
  
  if (Velocity[GridRank-1] == NULL) {
    ENZO_FAIL("malloc error (out of memory?)\n");
  }
 
  /* Compute the decrease in mass for particles moving to this grid
     (We assume here all grids are from the same level). */
 
  float RefinementFactors[MAX_DIMENSION];
  this->ComputeRefinementFactorsFloat(FromGrid[0], RefinementFactors);
  float MassDecrease = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    MassDecrease *= RefinementFactors[dim];
  MassDecrease = 1.0/MassDecrease;
 
  /* Copy this grid's particles to the new space. */
 
  for (i = 0; i < NumberOfParticles; i++) {
    Mass[i]   = ParticleMass[i];
    Number[i] = ParticleNumber[i];
    Type[i]   = ParticleType[i];
  }
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < NumberOfParticles; i++) {
      Position[dim][i] = ParticlePosition[dim][i];
      Velocity[dim][i] = ParticleVelocity[dim][i];
    }
  for (j = 0; j < NumberOfParticleAttributes; j++)
    for (i = 0; i < NumberOfParticles; i++)
      Attribute[j][i] = ParticleAttribute[j][i];
 
  /* Delete this grid's particles (now copied). */
 
  this->DeleteParticles();
 
  /* Copy new pointers into their correct position. */
 
  this->SetParticlePointers(Mass, Number, Type, Position, Velocity,
			    Attribute);
 
  /* Copy FromGrids' particles to new space on local "fake" grid. */
 
  int Index = NumberOfParticles;
  for (grid = 0; grid < NumberOfGrids; grid++) {

    for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++) {
      const int is_bh_mass_absolute =
	ParticleHasAbsoluteBHMass(FromGrid[grid]->ParticleType[i]);
      const float mass_before = FromGrid[grid]->ParticleMass[i];
      const float mass_after = is_bh_mass_absolute ?
	mass_before : mass_before * MassDecrease;
      Mass[Index+i] = mass_after;
      Number[Index+i] = FromGrid[grid]->ParticleNumber[i];
      Type[Index+i] = FromGrid[grid]->ParticleType[i];
      ReportBHParticleAMRMassTransfer("Grid_MoveAllParticles",
				      "child_to_parent",
				      FromGrid[grid]->GetLevel(),
				      this->GetLevel(),
				      "MassDecrease", MassDecrease,
				      FromGrid[grid]->ParticleType[i],
				      FromGrid[grid]->ParticleNumber[i],
				      mass_before, mass_after,
				      FromGrid[grid]->ParticleAttribute, i,
				      "bh_preserved_absolute_mass");
    }
    
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++) {
	Position[dim][Index+i] = FromGrid[grid]->ParticlePosition[dim][i];
	Velocity[dim][Index+i] = FromGrid[grid]->ParticleVelocity[dim][i];
      }

    for (j = 0; j < NumberOfParticleAttributes; j++)
      for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++)
	Attribute[j][Index+i] = FromGrid[grid]->ParticleAttribute[j][i];
    
    Index += FromGrid[grid]->NumberOfParticles;

  } // ENDFOR grids 
  
  /* Set new number of particles in this grid. */
 
  NumberOfParticles = TotalNumberOfParticles;
 
  /* Delete FromGrid's particles (and set number of particles to zero). */
 
  for (grid = 0; grid < NumberOfGrids; grid++) {
    FromGrid[grid]->NumberOfParticles = 0;
    FromGrid[grid]->DeleteParticles();
  }

  /******************** ACTIVE PARTICLES ********************/

  ActiveParticleList<ActiveParticleType> MoveParticles(NumberOfSubgridActiveParticles);

  int dlevel = nint(logf(RefinementFactors[0]) / logf(RefineBy));  /* Fix by JHW */
  //int dlevel = logf(RefinementFactors[0]) / logf(RefineBy);

  for (grid = 0; grid < NumberOfGrids; grid++) {
    for (i = 0; i < FromGrid[grid]->NumberOfActiveParticles; i++) {
      FromGrid[grid]->ActiveParticles[i]->AdjustMassByFactor(MassDecrease);
      FromGrid[grid]->ActiveParticles[i]->ReduceLevel(dlevel);
      FromGrid[grid]->ActiveParticles[i]->AssignCurrentGrid(this);
      FromGrid[grid]->ActiveParticles[i]->SetGridID(this->ID);
      MoveParticles.copy_and_insert(*(FromGrid[grid]->ActiveParticles[i]));
    }
    FromGrid[grid]->DeleteActiveParticles();
  }

  this->AddActiveParticles(MoveParticles, 0, NumberOfSubgridActiveParticles);
 
  return SUCCESS;
}











/* Old version of MoveAllParticles : 
   might just work for less computationally intensive runs */

int grid::MoveAllParticlesOld(int NumberOfGrids, grid* FromGrid[])  
{

  if (NumberOfGrids < 1) {
    fprintf(stderr, "NumberOfGrids(%d) must be > 0.\n", NumberOfGrids);
    return FAIL;
  }

  /* Determine total number of particles. */

  int TotalNumberOfParticles = NumberOfParticles;
  int i, j, grid, dim, *Type;
  PINT *Number;

  for (grid = 0; grid < NumberOfGrids; grid++)
    TotalNumberOfParticles += FromGrid[grid]->NumberOfParticles;
  if (TotalNumberOfParticles == 0)
    return SUCCESS;

  /* Debugging info. */

  if (debug) printf("MoveAllParticles: %d (before: ThisGrid = %d).\n",
		    TotalNumberOfParticles, NumberOfParticles);

  /* Allocate space for the particles. */

  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION], *Mass,
        *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];

  if (MyProcessorNumber == ProcessorNumber) {
     Mass = new float[TotalNumberOfParticles];
     Number = new PINT[TotalNumberOfParticles]; 
     Type = new int[TotalNumberOfParticles];
     for (int dim = 0; dim < GridRank; dim++) {
       Position[dim] = new FLOAT[TotalNumberOfParticles];
       Velocity[dim] = new float[TotalNumberOfParticles];
     }
     for (int i = 0; i < NumberOfParticleAttributes; i++)
       Attribute[i] = new float[TotalNumberOfParticles];

     if (Velocity[GridRank-1] == NULL) {
       fprintf(stderr, "malloc error (out of memory?)\n");
       return FAIL;
     }
  }

  /* Compute the decrease in mass for particles moving to this grid
     (We assume here all grids are from the same level). */

  float RefinementFactors[MAX_DIMENSION];
  this->ComputeRefinementFactorsFloat(FromGrid[0], RefinementFactors);
  float MassDecrease = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    MassDecrease *= RefinementFactors[dim];
  MassDecrease = 1.0/MassDecrease;

  /* Copy this grid's particles to the new space. */

   if (MyProcessorNumber == ProcessorNumber) {
     for (i = 0; i < NumberOfParticles; i++) {
       Mass[i]   = ParticleMass[i];
       Number[i] = ParticleNumber[i];
       Type[i]   = ParticleType[i];
     }
     for (dim = 0; dim < GridRank; dim++)
       for (i = 0; i < NumberOfParticles; i++) {
	 Position[dim][i] = ParticlePosition[dim][i];
	 Velocity[dim][i] = ParticleVelocity[dim][i];
       }
     for (j = 0; j < NumberOfParticleAttributes; j++)
       for (i = 0; i < NumberOfParticles; i++)
	 Attribute[j][i] = ParticleAttribute[j][i];
   }

  /* Delete this grid's particles (now copied). */  

  if (MyProcessorNumber == ProcessorNumber) {
    this->DeleteParticles();

    /* Copy new pointers into their correct position. */

    this->SetParticlePointers(Mass, Number, Type, Position, Velocity,  
			      Attribute);
  }

  /* Copy FromGrids' particles to new space (starting at NumberOfParticles). */

  int Index = NumberOfParticles;
  for (grid = 0; grid < NumberOfGrids; grid++) {

   /* If on the same processor, just copy. */

    if (MyProcessorNumber == ProcessorNumber &&
        MyProcessorNumber == FromGrid[grid]->ProcessorNumber) {

      //      fprintf(stderr, "P(%d) copying %d particles\n", MyProcessorNumber,
      //	     FromGrid[grid]->NumberOfParticles);

      for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++) {
	const int is_bh_mass_absolute =
	  ParticleHasAbsoluteBHMass(FromGrid[grid]->ParticleType[i]);
	const float mass_before = FromGrid[grid]->ParticleMass[i];
	const float mass_after = is_bh_mass_absolute ?
	  mass_before : mass_before * MassDecrease;
	Mass[Index+i] = mass_after;
	Number[Index+i] = FromGrid[grid]->ParticleNumber[i];
	Type[Index+i] = FromGrid[grid]->ParticleType[i];
	ReportBHParticleAMRMassTransfer("Grid_MoveAllParticlesOld",
					"child_to_parent",
					FromGrid[grid]->GetLevel(),
					this->GetLevel(),
					"MassDecrease", MassDecrease,
					FromGrid[grid]->ParticleType[i],
					FromGrid[grid]->ParticleNumber[i],
					mass_before, mass_after,
					FromGrid[grid]->ParticleAttribute, i,
					"bh_preserved_absolute_mass");
      }

      for (dim = 0; dim < GridRank; dim++)
	for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++) {
	  Position[dim][Index+i] = FromGrid[grid]->ParticlePosition[dim][i];
	  Velocity[dim][Index+i] = FromGrid[grid]->ParticleVelocity[dim][i];
	}
      for (j = 0; j < NumberOfParticleAttributes; j++)
	for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++)
	  Attribute[j][Index+i] = FromGrid[grid]->ParticleAttribute[j][i];
    }

    /* Otherwise, communicate. */

    else {
      if (MyProcessorNumber == ProcessorNumber ||
          MyProcessorNumber == FromGrid[grid]->ProcessorNumber)
	if (FromGrid[grid]->CommunicationSendParticles(this, ProcessorNumber,
              0, FromGrid[grid]->NumberOfParticles, Index) == FAIL) {
	  fprintf(stderr, "Error in grid->CommunicationSendParticles.\n");
	  return FAIL;
        }

      /* Change mass, as required. */

      if (MyProcessorNumber == ProcessorNumber)
	for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++) {
	  const int is_bh_mass_absolute =
	    ParticleHasAbsoluteBHMass(Type[Index+i]);
	  const float mass_before = Mass[Index+i];
	  const float mass_after = is_bh_mass_absolute ?
	    mass_before : mass_before * MassDecrease;
	  Mass[Index+i] = mass_after;
	  ReportBHParticleAMRMassTransfer("Grid_MoveAllParticlesOld",
					  "child_to_parent",
					  FromGrid[grid]->GetLevel(),
					  this->GetLevel(),
					  "MassDecrease", MassDecrease,
					  Type[Index+i], Number[Index+i],
					  mass_before, mass_after,
					  Attribute, Index+i,
					  "bh_preserved_absolute_mass");
	}

    }

    Index += FromGrid[grid]->NumberOfParticles;

  } // end: loop over grids.

  NumberOfParticles = TotalNumberOfParticles; 

  /* Delete FromGrid's particles (and set number of particles to zero). */

  for (grid = 0; grid < NumberOfGrids; grid++) {
    FromGrid[grid]->NumberOfParticles = 0;
    if (MyProcessorNumber == FromGrid[grid]->ProcessorNumber)

      FromGrid[grid]->DeleteParticles();
  }

  return SUCCESS;
}
