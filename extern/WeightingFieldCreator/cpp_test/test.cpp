#include <meep.hpp>
#include <cmath>
using namespace meep;

double eps(const vec &p) {
  if (p.x() < 2 && p.y() < 3)
    return 12.0;
  return 1.0;
}

std::complex<double> srcfunc(double t, void*) {
  return std::exp(-50 * std::pow(t - 10, 2));
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv); // do this even for non-MPI Meep
  double resolution = 20;     // pixels per distance
  grid_volume v = vol2d(5, 10, resolution); // 5x10 2d cell
  structure s(v, eps, pml(1.0));
  fields f(&s);

  f.output_hdf5(Dielectric, v.surroundings());

  custom_src_time src(srcfunc, NULL);
  f.add_point_source(Ez, src, vec(1.1, 2.3));
  while (f.time() < 20) {
     f.step();
  }

  f.output_hdf5(Ez, v.surroundings());

  return 0;
}
