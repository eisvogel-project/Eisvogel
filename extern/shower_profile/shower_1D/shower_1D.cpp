#include "shower_1D.h"
#include <iostream>
#include <fstream>
#include <array>
#include <math.h>
#include "constants.h"
#include "units.h"
#include <vector>
#include "Current.hh"

showers::Shower1D::Shower1D(std::array<float, 3> pos, double en, double ze,
                            double az, ChargeExcessProfile ce,
                            double ce_scaling, environment::IceProfile& ice)
{
  start_position = pos;
  energy = en;
  azimuth = az;
  zenith = ze;
  charge_excess_profile = ce;
  charge_excess_profile_scaling = ce_scaling;
  ice_profile = ice;
}

void showers::Shower1D::get_shower(double delta_t, std::vector<double>* t,
                                   std::vector<double>* x,
                                   std::vector<double>* y,
                                   std::vector<double>* z,
                                   std::vector<double>* ce)
{
  int n_points = 0;
  double delta_s = delta_t * constants::c;
  double max_grammage =
      charge_excess_profile.grammage[charge_excess_profile.grammage.size() - 1];
  double integrated_grammage = 0;
  double xx = start_position[0];
  double yy = start_position[1];
  double zz = start_position[2];
  while (integrated_grammage < max_grammage) {
    integrated_grammage =
        integrated_grammage + ice_profile.get_density(xx, yy, zz) * delta_s;
    n_points++;
    xx = xx + delta_s * sin(zenith) * cos(azimuth);
    yy = yy + delta_s * sin(zenith) * sin(azimuth);
    zz = zz + delta_s * cos(zenith);
  }
  t->resize(n_points);
  x->resize(n_points);
  y->resize(n_points);
  z->resize(n_points);
  ce->resize(n_points);
  (*t)[0] = 0;
  (*x)[0] = start_position[0];
  (*y)[0] = start_position[1];
  (*z)[0] = start_position[2];
  (*ce)[0] = charge_excess_profile.charge_excess[0];
  integrated_grammage = 0;
  int grammage_i = 0;
  double delta_ce;
  double delta_grammage;

  for (int i = 1; i < n_points; i++) {
    (*x)[i] = (*x)[i - 1] + delta_s * sin(zenith) * cos(azimuth);
    (*y)[i] = (*y)[i - 1] + delta_s * sin(zenith) * sin(azimuth);
    (*z)[i] = (*z)[i - 1] + delta_s * cos(zenith);
    (*t)[i] = (*t)[i - 1] + delta_t;
    integrated_grammage = integrated_grammage +
                          delta_s * ice_profile.get_density(
                                        (*x)[i - 1], (*y)[i - 1], (*z)[i - 1]);
    while (grammage_i < charge_excess_profile.grammage.size() &&
           charge_excess_profile.grammage[grammage_i] < integrated_grammage) {
      grammage_i++;
    }
    if (grammage_i == 0) {
      (*ce)[i] = charge_excess_profile.charge_excess[0] *
                 charge_excess_profile_scaling;
    }
    else {
      delta_ce = charge_excess_profile.charge_excess[grammage_i] -
                 charge_excess_profile.charge_excess[grammage_i - 1];
      delta_grammage = (integrated_grammage -
                        charge_excess_profile.grammage[grammage_i - 1]) /
                       (charge_excess_profile.grammage[grammage_i] -
                        charge_excess_profile.grammage[grammage_i - 1]);
      (*ce)[i] = (charge_excess_profile.charge_excess[grammage_i - 1] +
                  delta_ce * delta_grammage) *
                 charge_excess_profile_scaling;
    }
  }
}

void showers::Shower1D::fill_tracks(double delta_t,
                                    std::vector<LineCurrentSegment>& out)
{
  std::vector<double> t;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> ce;

  get_shower(delta_t, &t, &x, &y, &z, &ce);

  for (int i = 0; i < t.size() - 1; i++) {

    scalar_t c = constants::c;
    LineCurrentSegment cur_track(
        XYZCoordVector{(scalar_t)x[i] / c, (scalar_t)y[i] / c,
                       (scalar_t)z[i] / c}, // track start position
        XYZCoordVector{(scalar_t)x[i + 1] / c, (scalar_t)y[i + 1] / c,
                       (scalar_t)z[i + 1] / c}, // track end position
        (scalar_t)t[i], (scalar_t)t[i + 1], (scalar_t)ce[i]);
    out.push_back(cur_track);
  }
}

void showers::Shower1D::print_dimensions()
{
  std::vector<double> t;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> ce;
  get_shower(1., &t, &x, &y, &z, &ce);
  double r_max = sqrt(x[0] * x[0] + y[0] * y[0] + z[0] * z[0]);
  double r_min = sqrt(x[0] * x[0] + y[0] * y[0] + z[0] * z[0]);
  double z_max = z[0];
  double z_min = z[0];
  double r;
  for (int i = 0; i < t.size(); i++) {
    r = sqrt(x[i] * x[i] + y[i] * y[i]);
    r_min = std::min(r_min, r);
    r_max = std::max(r_max, r);
    z_min = std::min(z_min, z[i]);
    z_max = std::max(z_max, z[i]);
  }
  r_max /= constants::c;
  r_min /= constants::c;
  z_max /= constants::c;
  z_min /= constants::c;

  std::cout << "Required weighting field size:\n";
  std::cout << r_min << "< r < " << r_max << "\n";
  std::cout << z_min << "< z < " << z_max << "\n";
}
