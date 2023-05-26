#include <meep.hpp>
#include <cmath>
using namespace meep;

double eps(const vec &p) {
  if (p.z() < 5)
    return 12.0;
  return 1.0;
}

std::complex<double> srcfunc(double t, void*) {
  return std::exp(-50 * std::pow(t - 10, 2));
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  double resolution = 20;
  grid_volume v = volcyl(5, 10, resolution); // (r, z, resolution)
  
  structure s(v, eps, pml(1.0));
  fields f(&s);

  f.output_hdf5(Dielectric, v.surroundings());

  custom_src_time src(srcfunc, NULL);
  f.add_point_source(Ez, src, veccyl(0.0, 4.0));
  
  while (f.time() < 13) {
     f.step();
  }

  f.output_hdf5(Ez, v.surroundings());

  return 0;
}
