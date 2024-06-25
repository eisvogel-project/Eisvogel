#ifndef DENSITYPROFILE_CLASS
#define DENSITYPROFILE_CLASS

namespace environment {
  class IceProfile {
  public:
    double get_density(float x, float y, float z);
    double get_index_of_refraction(float x, float y, float z);
    double get_maximum_index_of_refraction();
  };
} // namespace environment

#endif
