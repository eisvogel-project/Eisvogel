#include "units.h"

namespace constants {
  const double speed_of_light = 299792458 * units::meter / units::second;
  const double c = speed_of_light;

  const double vacuum_permittivity = 55.26349406 * units::e_charge *
                                     units::e_charge / units::eV /
                                     units::micrometer;
  const double epsilon_0 = vacuum_permittivity;
} // namespace constants
