#ifndef BH_ACCRETION_SAFETY_H
#define BH_ACCRETION_SAFETY_H

#include <math.h>

static const double BHAccretionMaxCellRemovalFraction = 0.9;
static const double BHAccretionSmallRhoFloorFactor = 1.0e-10;

static inline double BHAccretionMaxDouble(double a, double b)
{
  return (a > b) ? a : b;
}

static inline double BHAccretionDensityFloorCode(double tiny_number_code,
                                                 double small_rho_code)
{
  double floor = 0.0;
  if (isfinite(tiny_number_code) && tiny_number_code > floor)
    floor = tiny_number_code;
  if (isfinite(small_rho_code) && small_rho_code > 0.0)
    floor = BHAccretionMaxDouble(floor,
                                 BHAccretionSmallRhoFloorFactor * small_rho_code);
  return floor;
}

static inline double BHAccretionSafeRemovableMass(double cell_mass,
                                                  double density_floor_mass)
{
  if (!isfinite(cell_mass) || cell_mass <= 0.0)
    return 0.0;

  const double fraction_floor =
    (1.0 - BHAccretionMaxCellRemovalFraction) * cell_mass;
  const double mass_floor =
    BHAccretionMaxDouble(BHAccretionMaxDouble(0.0, density_floor_mass),
                         fraction_floor);

  if (mass_floor >= cell_mass)
    return 0.0;

  return cell_mass - mass_floor;
}

static inline double BHAccretionInternalEnergySpecific(double total_energy,
                                                       double gas_energy,
                                                       double vx,
                                                       double vy,
                                                       double vz,
                                                       bool use_dual_energy)
{
  if (use_dual_energy)
    return gas_energy;

  return total_energy - 0.5*(vx*vx + vy*vy + vz*vz);
}

static inline double BHAccretionPressureFromState(double density,
                                                  double total_energy,
                                                  double gas_energy,
                                                  double vx,
                                                  double vy,
                                                  double vz,
                                                  bool use_dual_energy,
                                                  double gamma)
{
  if (!isfinite(density) || density <= 0.0)
    return -1.0;

  const double eint =
    BHAccretionInternalEnergySpecific(total_energy, gas_energy,
                                      vx, vy, vz, use_dual_energy);
  if (!isfinite(eint))
    return -1.0;

  return (gamma - 1.0) * density * eint;
}

static inline double BHAccretionGasEnergyFloor(double density,
                                               double pressure_floor,
                                               double gamma)
{
  if (!isfinite(density) || density <= 0.0 || gamma <= 1.0)
    return 0.0;

  double floor = 0.0;
  if (isfinite(pressure_floor) && pressure_floor > 0.0)
    floor = pressure_floor / ((gamma - 1.0) * density);

  return BHAccretionMaxDouble(floor, 0.0);
}

#endif
