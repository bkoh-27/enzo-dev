#include <cmath>
#include <cstdlib>
#include <iostream>

#include "../../../src/enzo/BHAccretionSafety.h"

static void require_true(bool ok, const char *message)
{
  if (!ok) {
    std::cerr << message << std::endl;
    std::exit(1);
  }
}

static void require_close(double got, double expected, double tol,
                          const char *message)
{
  if (std::fabs(got - expected) > tol) {
    std::cerr << message << ": got=" << got
              << " expected=" << expected << std::endl;
    std::exit(1);
  }
}

int main()
{
  const double mcell = 1.0;
  const double density_floor_mass = 1.0e-20;
  const double safe =
    BHAccretionSafeRemovableMass(mcell, density_floor_mass);

  require_close(safe, 0.9, 1.0e-14,
                "safe removable mass must leave the configured cell fraction");
  require_true(safe < mcell,
               "safe removable mass must never consume a whole positive cell");

  const double large_floor_safe =
    BHAccretionSafeRemovableMass(mcell, 0.25);
  require_close(large_floor_safe, 0.75, 1.0e-14,
                "density floor must further reduce removable mass");

  const double rejected =
    BHAccretionSafeRemovableMass(1.0e-30, 1.0e-20);
  require_close(rejected, 0.0, 0.0,
                "cells at or below the density floor must be rejected");

  const double pressure =
    BHAccretionPressureFromState(1.0, 10.0, 1.0, 2.0, 3.0, -1.0,
                                 true, 5.0/3.0);
  require_true(pressure > 0.0,
               "dual-energy pressure must come from positive gas energy");

  const double nondual_pressure =
    BHAccretionPressureFromState(1.0, 10.0, -1.0, 1.0, 2.0, 3.0,
                                 false, 5.0/3.0);
  require_true(nondual_pressure > 0.0,
               "non-dual pressure must come from total minus kinetic energy");

  return 0;
}
