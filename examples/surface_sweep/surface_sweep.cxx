#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <sstream>
//#include <format>

#include "Common.hh"
#include "SignalCalculator.hh"
#include "Current.hh"



//Want this function for picking zenith angles using inverse CDF (no analytical solution for cos^2 pdf)
double findZero(std::function<double(double)> f, double low, double high, double tolerance = 1e-6) {
    double mid;
    while (high - low > tolerance) {
        mid = low + (high - low) / 2;
        if (f(mid) == 0) {
            return mid;
        } else if (f(mid) * f(low) < 0) {
            high = mid;
        } else {
            low = mid;
        }
    }
    return low;
}

int main(int argc, char* argv[]) {
//first argument is greens function
  if(argc < 2) {
    throw;
  }

  //Coordinate origin is at phased array location, distances in natural feet

  std::string gf_path = argv[1];
  SignalCalculator calc(gf_path, 20);

  //Make a vertical Gaussian charge excess with sigma=sigma centered at z0

  //These parameters for testing a shower in inhomogeneous ice:
  scalar_t sigma = 15.0f;

  //Set this parameter for a shower peaking at the surface
  scalar_t d_peak = 0;
  scalar_t tstart = -70, tend = 60;
  scalar_t beta = 0.999f;

  //d is distance from surface
  auto charge = [&](scalar_t d){
    if(d > 0) {
      return 1.0;
    } else {
      return exp(-(d-d_peak)*(d-d_peak) / (sigma*sigma));
    }
  };

  scalar_t bmax = 450.0;
  scalar_t bmin = 100.0;

  scalar_t bstep = 5.0;
  //scalar_t bstep = 10.0;
  std::vector<float> b_arr;
  for (float i = bmin; i < bmax; i += bstep) {
        b_arr.push_back(i);
        //std::cout << "bval " << i << std::endl;
    }

  //scalar_t phi_step = 90;
  scalar_t phi_step = 30;
  std::vector<float> phi_arr;
  for (int i = 0; i < 360; i += phi_step) {
        phi_arr.push_back(i * M_PI / 180.0);
    } 
  
  scalar_t zen_val = 0;
  scalar_t zen_cdf_step = 0.05;
  std::vector<float> zen_arr;
  for (float cdfval = 0.01; cdfval < 0.9; cdfval += zen_cdf_step) {
        std::function<double(double)> cdf = [cdfval](double x) {return std::sin(2*x) + 2*x - M_PI*cdfval; };
        zen_val = findZero(cdf, 0, M_PI);
        zen_arr.push_back(zen_val);
        std::cout << "zenith: " << zen_val*180.0/M_PI << std::endl;
  }

  scalar_t t_sig_samp = 0.5;
  std::size_t num_samples = 400;
  scalar_t delta_t = 0.5;
  scalar_t t_sig_start;
  //Loop over all radial impact locations

  std::vector<scalar_t> signal_values(num_samples);
  for (float b : b_arr) {
    
    t_sig_start = 1.7*sqrt(580.0*580.0 + b*b);
    std::cout << "Throwing showers at b = " << b << std::endl;

    for (float zen : zen_arr) {
      for (float phi : phi_arr) {
        
        
        for(scalar_t cur_t = tstart; cur_t < tend; cur_t += delta_t) {
          LineCurrentSegment cur_track(XYZCoordVector{b + -beta * cur_t * 
              std::sin(zen)*std::cos(phi), beta * cur_t * std::sin(zen)*std::sin(phi),
              -beta * cur_t*std::cos(zen)},               // track start position
            XYZCoordVector{b + -beta * (cur_t + delta_t) * 
              std::sin(zen)*std::cos(phi), beta * (cur_t + delta_t) * std::sin(zen)*std::sin(phi),
              -beta * (cur_t + delta_t) *std::cos(zen)},   // track end position
              cur_t, cur_t + delta_t, charge(-beta * cur_t));

          //std::cout << charge(-beta * cur_t) << std::endl;
          calc.AccumulateSignal(cur_track, t_sig_start, t_sig_samp, num_samples, signal_values);       
        }

        std::stringstream ss;
        ss << "rad_" << std::round(b) << "_zen" << std::round(zen*180/M_PI) << "_phi" << std::round(phi*180/M_PI) << ".csv";
        std::string fname = ss.str();

        std::string directory = "/home/nalden/Eisvogel/build/surface_sweep_files/";
        std::string filepath = directory + fname;

        //std::cout << filepath << std::endl;

        std::ofstream myfile;
        myfile.open(filepath);

        scalar_t i = 0;
        for(scalar_t& cur: signal_values) {
          //std::cout << cur << std::endl;
          myfile << t_sig_start + i*t_sig_samp << ", " << cur<< "\n";
          i++;
        }
        myfile.close();
        
      }
    }
  }

  return 0;
}
