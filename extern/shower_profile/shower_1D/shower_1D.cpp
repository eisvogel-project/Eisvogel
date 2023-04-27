#include "shower_1D.h"
#include <iostream>
#include <array>
#include <math.h>
#include "Eisvogel/CoordUtils.hh"
#include "Eisvogel/Current0D.hh"
#include "constants.h"
#include "units.h"
#include <vector>


showers::Shower1D::Shower1D(
	std::array<float, 3> pos,
	double en,
	double ze,
	double az,
	ChargeExcessProfile ce,
	double ce_scaling,
	environment::IceProfile &ice
) {
	start_position = pos;
	energy = en;
	azimuth = az;
	zenith = ze;
	charge_excess_profile = ce;
	charge_excess_profile_scaling = ce_scaling;
	ice_profile = ice;
}

void showers::Shower1D::get_shower(
	double delta_t,
	std::vector<double> *t,
	std::vector<double> *x,
	std::vector<double> *y,
	std::vector<double> *z,
	std::vector<double> *ce
) {
	int n_points = 0;
	double delta_s = delta_t * constants::c;
	double max_grammage = charge_excess_profile.grammage[charge_excess_profile.grammage.size() - 1];
	double integrated_grammage = 0;
	double xx = start_position[0];
	double yy = start_position[1];
	double zz = start_position[2];
	while (integrated_grammage < max_grammage) {
		integrated_grammage = integrated_grammage + ice_profile.get_density(
				xx,
				yy,
				zz
		) * delta_s;
		n_points ++;
		xx = xx + delta_s * sin(zenith) * cos(azimuth);
		yy = yy+ delta_s * sin(zenith) * sin(azimuth);
		zz = zz + delta_s * cos(zenith);
	}
	t -> resize(n_points);
	x -> resize(n_points);
	y -> resize(n_points);
	z -> resize(n_points);
	ce -> resize(n_points);
	(*t)[0] = 0;
	(*x)[0] = start_position[0];
	(*y)[0] = start_position[1];
	(*z)[0] = start_position[2];
	(*ce)[0] = charge_excess_profile.charge_excess[0];
	integrated_grammage = 0;
	int grammage_i = 0;
	double delta_ce;
	double delta_grammage;
	for (int i=1; i < n_points; i++) {
		(*x)[i] = (*x)[i - 1] + delta_s * sin(zenith) * cos(azimuth);
		(*y)[i] = (*y)[i - 1] + delta_s * sin(zenith) * sin(azimuth);
		(*z)[i] = (*z)[i - 1] + delta_s * cos(zenith);
		(*t)[i] = (*t)[i - 1] + delta_t;
		integrated_grammage = integrated_grammage + delta_s * ice_profile.get_density(
			(*x)[i-1],
			(*y)[i-1],
			(*z)[i-1]
		);
		while(grammage_i < charge_excess_profile.grammage.size() && charge_excess_profile.grammage[grammage_i] < integrated_grammage) {
			grammage_i ++;
		}
		if (grammage_i == 0) {
			(*ce)[i] = charge_excess_profile.charge_excess[0] * charge_excess_profile_scaling;
		} else {
			delta_ce = charge_excess_profile.charge_excess[grammage_i] - charge_excess_profile.charge_excess[grammage_i - 1];
			delta_grammage = (integrated_grammage - charge_excess_profile.grammage[grammage_i - 1]) / (charge_excess_profile.grammage[grammage_i] - charge_excess_profile.grammage[grammage_i - 1]);
			(*ce)[i] = (charge_excess_profile.charge_excess[grammage_i - 1]+ delta_ce * delta_grammage) * charge_excess_profile_scaling;
		}
	}
}

Current0D showers::Shower1D::get_current(
		double delta_t
){
        std::vector<double> t;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
	std::vector<double> ce;
	get_shower(
			delta_t,
			&t,
			&x,
			&y,
			&z,
			&ce
	);
        std::vector<CoordVector> positions;
        std::vector<scalar_t> charge_excess;
        for (int i = 0; i < t.size(); i++) {
            positions.push_back(
                    CoordUtils::MakeCoordVectorTXYZ(
                    t[i],
                    x[i] / constants::c,
                    y[i] / constants::c,
                    z[i] / constants::c
                    )
                );
            if (i < t.size() -1) {
                charge_excess.push_back(ce[i]);
            }
        }
        Current0D current(positions, charge_excess);
        return current;
}




