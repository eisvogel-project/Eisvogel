#include "ice_profile.h"
#include "units.h"

double environment::IceProfile::get_density(float x, float y, float z)
{
  return 917 * units::kg / units::cubic_meter;
}
double environment::IceProfile::get_index_of_refraction(float x, float y,
                                                        float z)
{
  return 1.77;
}

double environment::IceProfile::get_maximum_index_of_refraction()
{
  return 1.77;
}
