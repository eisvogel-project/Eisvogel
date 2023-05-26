#include <meep.hpp>
#include <cmath>
using namespace meep;

double eps(const vec &p) {
  if (p.z() < 10)
    return 2.0;
  return 1.0;
}

std::complex<double> srcfunc(double t, void*) {
  if(t > 20) {
    return 0.0;
  }
  else {
    return std::exp(-50 * std::pow(t - 10, 2));
  }
  return 0.0;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  double resolution = 20;
  grid_volume v = volcyl(20, 20, resolution); // (r, z, resolution)
  
  structure s(v, eps, pml(1.0), identity(), 0, 0.5);
  fields f(&s);

  f.output_hdf5(Dielectric, v.surroundings());

  custom_src_time src(srcfunc, NULL);
  f.add_point_source(Ez, src, veccyl(0.0, 8.0));
  
  while (f.time() < 20) {
     f.step();
  }

  f.output_hdf5(Ez, v.surroundings());

  return 0;
}
