#ifndef CHARGEEXCESSPROFILE_CLASS
#define CHARGEEXCESSPROFILE_CLASS
#include <vector>

namespace showers {
class ChargeExcessProfile {
	public:
	std::vector<double> grammage;
	std::vector<double> charge_excess;
	int hadronic;
	int N;
	double energy;
};

// class ChargeExcessProfile2D {
// 	public:
//   ChargeExcessProfile2D(std::array<std::size_t, 2> size, scalar_t initial_val = 0) : charge_excess(Vector<std::size_t,2>(size), initial_val){
// 	}
//   NDVecArray<double, 2, 1> charge_excess;
// 	std::vector<double> grammage;
// 	std::vector<double> radius;
// 	int hadronic;
// 	int N;
// 	double energy;

// 	double get_charge_excess(double gram, double rad) {
// 		int radius_index = 0;
// 		for (int i_radius=1;i_radius < radius.size(); i_radius++) {
// 			radius_index = i_radius - 1;
// 			if (rad > radius[i_radius]) {
// 				break;
// 			}
// 		}
// 		int grammage_index = 0;
// 		for (int i_grammage = 1; i_grammage < grammage.size() - 1; i_grammage++) {
// 			grammage_index = i_grammage - 1;
// 			if (gram > grammage[i_grammage]) {
// 				break;
// 			}
// 		}
// 		double delta_grammage = grammage[grammage_index + 1] - grammage[grammage_index];
// 		double delta_ce_1 = charge_excess[{grammage_index + 1, radius_index}][0] - charge_excess[{grammage_index, radius_index}][0];
// 		double grammage_interpolation_1 = charge_excess[{grammage_index, radius_index}][0] + delta_ce_1 / delta_grammage * (gram - grammage[grammage_index]);

// 		double delta_ce_2 = charge_excess[{grammage_index + 1, radius_index + 1}][0] - charge_excess[{grammage_index, radius_index + 1}][0];
// 		double grammage_interpolation_2 = charge_excess[{grammage_index, radius_index + 1}][0] + delta_ce_2 / delta_grammage * (gram - grammage[grammage_index]);

// 		double delta_r = radius[radius_index + 1] - radius[radius_index];
// 		return grammage_interpolation_1 + (grammage_interpolation_2 - grammage_interpolation_1) / delta_r * (rad - radius[radius_index]);
// 	}
// };
  
}

#endif
