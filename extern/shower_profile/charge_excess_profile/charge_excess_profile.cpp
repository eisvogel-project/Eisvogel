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
	ChargeExcessProfile2D(std::array<std::size_t, 2> size, scalar_t initial_val = 0):charge_excess(size, initial_val){
	}
	DenseNDArray<scalar_t, 2> charge_excess;
	std::vector<scalar_t> grammage;
	std::vector<scalar_t> radius;
	int hadronic;
	int N;
	double energy;
};
}

#endif
