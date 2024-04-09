#ifndef SHOWER1D_CLASS
#define SHOWER1D_CLASS
#include <array>
#include "ice_profile.h"
#include "charge_excess_profile.cpp"
#include "Eisvogel/Current0DOld.hh"

namespace showers {
	class Shower1D {
	public:
		void get_shower(
				double delta_t,
				std::vector<double> *t,
				std::vector<double> *x,
				std::vector<double> *y,
				std::vector<double> *z,
				std::vector<double> *ce
		);
		Current0D get_current(
				double delta_t
		);
		Shower1D(
				std::array<float, 3> pos,
				double en,
				double ze,
				double az,
				ChargeExcessProfile ce,
				double ce_scaling,
				environment::IceProfile &ice
				);
		void print_dimensions();
	private:
		std::array<float, 3> start_position;
		double energy;
		double charge_excess_profile_scaling;
		ChargeExcessProfile charge_excess_profile;
		double zenith;
		double azimuth;
		float start_time;
		environment::IceProfile ice_profile;
	};
}

#endif
