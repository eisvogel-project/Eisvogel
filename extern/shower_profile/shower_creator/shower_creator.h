#ifndef SHOWERCREATOR_CLASS
#define SHOWERCREATOR_CLASS
#include "shower_1D.h"
#include "ice_profile.h"
#include "charge_excess_profile.cpp"
#include "Eisvogel/NDArray.hh"
#include <array>
#include <string>
#include <map>
#include <random>

namespace showers {
class ShowerCreator {
public:
	ShowerCreator(std::string file_path, bool is_2 = false);
	Shower1D create_shower(
			std::array<float,3> pos,
			double en,
			double zen,
			double az,
			int had
	);
	int read_shower(
		FILE *f,
		int *N,
		int *hadronic,
		double *E,
		std::vector<double> *grammage,
		std::vector<double> *q
		);
	ChargeExcessProfile2D read_shower(
		FILE *f,
		unsigned int &N
		);
private:
	std::string shower_file;
	environment::IceProfile density_profile;
	std::map<int, std::map<double, std::vector<showers::ChargeExcessProfile>>>ce_profiles;
	std::map<int, std::map<double, std::vector<showers::ChargeExcessProfile2D>>>ce_profiles_2D;
	std::map<int, std::vector<double>> stored_energies;
	bool is_2d_shower;
};
}

#endif
