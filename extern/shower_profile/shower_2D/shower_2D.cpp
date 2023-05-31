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

): charge_excess_profile(ce){
    starting_position = pos;
    energy = en;
    zenith = ze;
    azimuth = az;
    charge_excess_profile_scaling = ce_scaling;
    ice_profile = ice;
}

SparseCurrentDensity3D showers::Shower2D::get_current(scalar_t delta_t) {
    

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