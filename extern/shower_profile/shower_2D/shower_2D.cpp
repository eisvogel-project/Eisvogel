#include "shower_2D.h"
#include <iostream>
#include <fstream>
#include <array>
#include "charge_excess_profile.cpp"
#include "Eisvogel/CoordUtils.hh"
#include "Eisvogel/SparseCurrentDensity3D.hh"
#include "constants.h"
#include "units.h"

using namespace CoordUtils;

showers::Shower2D::Shower2D(
                std::array<scalar_t, 3> pos,
                scalar_t en,
                scalar_t ze,
                scalar_t az,
                ChargeExcessProfile2D ce,
                scalar_t ce_scaling,
                environment::IceProfile &ice
) : 
charge_excess_profile(ce),
starting_position(pos),
energy(en),
zenith(ze),
azimuth(az),
charge_excess_profile_scaling(ce_scaling),
ice_profile(ice)
{}

SparseCurrentDensity3D showers::Shower2D::get_current(DeltaVector voxel_size) {
    int n_points = 0;
    scalar_t delta_t = getT(voxel_size);
    double max_grammage = charge_excess_profile.grammage[charge_excess_profile.grammage.size() - 1];
    double max_radius = charge_excess_profile.radius[charge_excess_profile.radius.size() - 1];
    std::cout << "Max Radius: " << max_radius << "\n";
    double integrated_grammage = 0;
    double delta_s = delta_t * constants::c;
    double shower_t = 0;
     SparseCurrentDensity3D current_density(voxel_size);    
    std::array<double, 3> shower_pos = {starting_position[0], starting_position[1], starting_position[2]};
    std::array<double, 3> shower_axis = {sin(zenith) * cos(azimuth), sin(zenith) * sin(azimuth), cos(zenith)};
    std::array<double, 3> axis_normal_hor = {sin(azimuth), cos(azimuth), 0};
    std::array<double, 3> axis_normal_ver = {
        -cos(zenith) * cos(azimuth),
        cos(zenith) * sin(azimuth),
        sin(zenith) * cos(azimuth) * cos(azimuth) - sin(zenith) * sin(azimuth) * sin(azimuth)
        };
    double current_radius;
    double current_ce;
    while (integrated_grammage < max_grammage / 2) {
        integrated_grammage = integrated_grammage + ice_profile.get_density(
            shower_pos[0], shower_pos[1], shower_pos[2]
            ) * delta_s;
        shower_t = shower_t + delta_t;
        shower_pos[0] = shower_pos[0] + shower_axis[0] * delta_s;
        shower_pos[1] = shower_pos[1] + shower_axis[1] * delta_s;
        shower_pos[2] = shower_pos[2] + shower_axis[2] * delta_s;
        for (float hor_offset = -max_radius; hor_offset < max_radius; hor_offset = hor_offset + delta_s) {
            for (float ver_offset = -max_radius; ver_offset < max_radius; ver_offset = ver_offset + delta_s) {
                current_radius = sqrt(hor_offset * hor_offset + ver_offset * ver_offset);

                if (current_radius < max_radius) {
                    current_ce = charge_excess_profile.get_charge_excess(integrated_grammage, current_radius);
                    if (current_ce > 0) {
                        CoordVector current_pos = CoordUtils::MakeCoordVectorTXYZ(
                            shower_t,
                            (shower_pos[0] + hor_offset * axis_normal_hor[0] + ver_offset * axis_normal_ver[0]) / constants::c,
                            (shower_pos[1] + hor_offset * axis_normal_hor[1] + ver_offset * axis_normal_ver[1]) / constants::c,
                            (shower_pos[2] + hor_offset * axis_normal_hor[2] + ver_offset * axis_normal_ver[2]) / constants::c
                        );
                        FieldVector current_current = CoordUtils::MakeFieldVectorXYZ(
                            shower_axis[0] * current_ce,
                            shower_axis[1] * current_ce,
                            shower_axis[2] * current_ce
                        );
                        current_density.addCurrentElement(current_pos, current_current);
                    }
                }
            }
        }
    }
    return current_density;

}
void showers::Shower2D::dump_profile(){
    std::ofstream dump_file;
    std::cout << charge_excess_profile.radius.size() << "! \n";
    dump_file.open("profile2D.csv");
    /*
    dump_file << "-1, ";
    for (int j=0;j < charge_excess_profile.grammage.size();j++) {
        dump_file << charge_excess_profile.grammage[j] << ", ";
        }
    dump_file << "\n";
    */
    for (int i=0;i< charge_excess_profile.radius.size();i++) {
        //dump_file << charge_excess_profile.radius[i] << ", ";
        for (int j=0;j < charge_excess_profile.grammage.size();j++) {
            dump_file << charge_excess_profile.charge_excess(j, i) << ", ";
        }
        dump_file << "\n";
    }
    dump_file.close();
}