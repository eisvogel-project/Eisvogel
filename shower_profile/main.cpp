#include "shower_1D.h"
#include <array>
#include <iostream>
#include "shower_creator.h"
int main () {
	std::array<float, 3> pos = {1.0, 2.0, 3.0};

	showers::ShowerCreator shower_creator("./shower_file");
	showers::Shower1D new_shower = shower_creator.create_shower(
			pos,
			2.e18,
			1.,
			1.,
			1
	);
	std::vector<double> shower_x;
	std::vector<double> shower_y;
	std::vector<double> shower_z;
	std::vector<double> shower_t;
	std::vector<double> shower_ce;
	new_shower.get_shower(
		1.e-8,
		&shower_x,
		&shower_y,
		&shower_z,
		&shower_t,
		&shower_ce
	);

}
