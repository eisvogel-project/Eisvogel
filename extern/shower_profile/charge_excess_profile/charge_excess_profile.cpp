#ifndef CHARGEEXCESSPROFILE_CLASS
#define CHARGEEXCESSPROFILE_CLASS
#include <vector>
#include "Eisvogel/Common.hh"
#include "Eisvogel/NDArray.hh"

namespace showers {
class ChargeExcessProfile {
	public:
	std::vector<double> grammage;
	std::vector<double> charge_excess;
	int hadronic;
	int N;
	double energy;
};

class ChargeExcessProfile2D {
	public:
	ChargeExcessProfile2D(std::array<std::size_t, 2> size, scalar_t initial_val = 0):charge_excess_2d(size, initial_val), charge_excess_1d(size[0], initial_val){
	}
	DenseNDArray<double, 2> charge_excess_2d;
	DenseNDArray<double, 1> charge_excess_1d;
	std::vector<double> grammage;
	std::vector<double> radius;
	int hadronic;
	int N;
	double energy;
	double pi = 3.14159265;
	

	void set_charge_excess(int i_grammage, int i_radius, double value) {
		charge_excess_2d(i_grammage, i_radius) = value;
		double delta_grammage;
		double delta_radius;
		if (i_grammage > 0) {
			delta_grammage = grammage[i_grammage] - grammage[i_grammage - 1];
		} else {
			delta_grammage = grammage[1] - grammage[0];
		}
		if (i_radius > 0) {
			delta_radius = radius[i_radius] - radius[i_radius - 1];
		} else {
			delta_radius = radius[1] - radius[0];
		}
		charge_excess_1d(i_grammage) += value * pi * (2 * radius[i_radius] * delta_radius + delta_radius * delta_radius);
	}

	double get_charge_excess(double gram, double rad) {
		int radius_index = 0;
		//std::cout << "---> " << gram << ", " << grammage[2] << "\n";
		for (int i_radius=1;i_radius < radius.size(); i_radius++) {
			radius_index = i_radius - 1;
			if (rad < radius[i_radius]) {
				break;
			}
		}
		int grammage_index = 0;
		for (int i_grammage = 1; i_grammage < grammage.size() - 1; i_grammage++) {
			grammage_index = i_grammage - 1;
			if (gram < grammage[i_grammage]) {
				break;
			}
		}
		double delta_grammage = grammage[grammage_index + 1] - grammage[grammage_index];
		double delta_ce_1 = charge_excess_2d(grammage_index + 1, radius_index) - charge_excess_2d(grammage_index, radius_index);
		double grammage_interpolation_1 = charge_excess_2d(grammage_index, radius_index) + delta_ce_1 / delta_grammage * (gram - grammage[grammage_index]);

		double delta_ce_2 = charge_excess_2d(grammage_index + 1, radius_index + 1) - charge_excess_2d(grammage_index, radius_index + 1);
		double grammage_interpolation_2 = charge_excess_2d(grammage_index, radius_index + 1) + delta_ce_2 / delta_grammage * (gram - grammage[grammage_index]);
		double delta_r = radius[radius_index + 1] - radius[radius_index];
		return grammage_interpolation_1 + (grammage_interpolation_2 - grammage_interpolation_1) / delta_r * (rad - radius[radius_index]);
	}

	double get_charge_excess(double gram) {
		int grammage_index = 0;
		for (int i_grammage = 1; i_grammage < grammage.size() - 1; i_grammage++) {
			grammage_index = i_grammage - 1;
			if (gram < grammage[i_grammage]) {
				break;
			}
		}
		double delta_grammage = grammage[grammage_index + 1] - grammage[grammage_index];
		return charge_excess_1d(grammage_index) + (charge_excess_1d(grammage_index + 1) - charge_excess_1d(grammage_index)) * (gram - grammage[grammage_index]) / delta_grammage;
	}

};
}

#endif
