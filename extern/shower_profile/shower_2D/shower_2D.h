#ifndef SHOWER2D_CLASS
#define SHOWER2D_CLASS
#include <array>
#include "ice_profile.h"
#include "charge_excess_profile.cpp"
#include "Eisvogel/SparseCurrentDensity3D.hh"
#include "Eisvogel/Common.hh"
#include "Eisvogel/CoordUtils.hh"
using namespace CoordUtils;

namespace showers {
    class Shower2D {
        public:
            void get_shower(
                scalar_t delta_t,
                std::vector<scalar_t> *t,
                std::vector<scalar_t> *x,
                std::vector<scalar_t> *y,
                std::vector<scalar_t> *z,
                std::vector<scalar_t> *ce
            );

            SparseCurrentDensity3D get_current(
                DeltaVector voxel_size
            );
            SparseCurrentDensity3D get_current_1d(
                DeltaVector voxel_size
            );

            void dump_profile();

            Shower2D(
                std::array<scalar_t, 3> pos,
                scalar_t en,
                scalar_t ze,
                scalar_t az,
                ChargeExcessProfile2D ce,
                scalar_t ce_scaling,
                environment::IceProfile & ice
            );
        private:
            std::array<scalar_t, 3> starting_position;
            scalar_t energy;
            scalar_t zenith;
            scalar_t azimuth;
            ChargeExcessProfile2D charge_excess_profile;
            scalar_t charge_excess_profile_scaling;
            environment::IceProfile ice_profile;
    };

}

#endif