#ifndef SHOWERCREATOR_CLASS
#define SHOWERCREATOR_CLASS
#include "shower_1D.h"
// #include "shower_2D.h"
#include "ice_profile.h"
#include "charge_excess_profile.cpp"
// #include "NDVecArray.hh"
#include <array>
#include <string>
#include <map>
#include <random>

namespace showers {
  class ShowerCreator {
  public:
    ShowerCreator(std::string file_path);
    Shower1D create_shower(std::array<float, 3> pos, double en, double zen,
                           double az, int had, int i_shower = -1);
    int read_shower(FILE* f, int* N, int* hadronic, double* E,
                    std::vector<double>* grammage, std::vector<double>* q);

  private:
    std::string shower_file;
    environment::IceProfile density_profile;
    std::map<int, std::map<double, std::vector<showers::ChargeExcessProfile>>>
        ce_profiles;
    std::map<int, std::vector<double>> stored_energies;
  };

  // class ShowerCreator2D {
  // public:
  // 	ShowerCreator2D(std::string file_path);
  // 	Shower2D create_shower(
  // 			std::array<float,3> pos,
  // 			double en,
  // 			double zen,
  // 			double az,
  // 			int had
  // 	);
  // 	ChargeExcessProfile2D read_shower(
  // 		FILE *f,
  // 		int *N
  // 		);
  // private:
  // 	std::string shower_file;
  // 	environment::IceProfile density_profile;
  // 	std::map<int, std::map<double,
  // std::vector<showers::ChargeExcessProfile2D>>>ce_profiles; 	std::map<int,
  // std::vector<double>> stored_energies;
  // };

} // namespace showers

#endif
