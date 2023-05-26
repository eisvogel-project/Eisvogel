#include <meep.hpp>
using namespace meep;

double eps(const vec &p) {
  if (p.x() < 2 && p.y() < 3)
    return 12.0;
  return 1.0;
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv); // do this even for non-MPI Meep
  double resolution = 20;     // pixels per distance
  grid_volume v = vol2d(5, 10, resolution); // 5x10 2d cell
  structure s(v, eps, pml(1.0));
  fields f(&s);

  f.output_hdf5(Dielectric, v.surroundings());

  double freq = 0.3, fwidth = 0.1;
  gaussian_src_time src(freq, fwidth);
  f.add_point_source(Ey, src, vec(1.1, 2.3));
  while (f.time() < f.last_source_time()) {
     f.step();
  }

  f.output_hdf5(Hz, v.surroundings());

  return 0;
}
