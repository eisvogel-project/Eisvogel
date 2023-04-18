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
}

#endif
