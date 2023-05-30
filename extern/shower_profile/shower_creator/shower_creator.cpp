#include "Eisvogel/Common.hh"
#include "shower_creator.h"
#include "shower_1D.h"
#include "charge_excess_profile.cpp"
#include "units.h"
#include <iostream>
#include <array>
#include <stdio.h>
#include <time.h>
#include <stdexcept>
#include <random>


showers::ShowerCreator::ShowerCreator(
		std::string file_path,
		bool is_2d
		) {
	shower_file = file_path;
	density_profile = environment::IceProfile();
	int reading_complete = 0;
	unsigned int N;
	FILE * pFile = fopen(shower_file.c_str(), "rb");
	is_2d_shower = is_2d;
	ce_profiles[0] = std::map<double, std::vector<showers::ChargeExcessProfile>>();
	ce_profiles[1] = std::map<double, std::vector<showers::ChargeExcessProfile>>();
	ce_profiles_2D[0] = std::map<double, std::vector<showers::ChargeExcessProfile2D>>();
	ce_profiles_2D[1] = std::map<double, std::vector<showers::ChargeExcessProfile2D>>();
	while (reading_complete == 0) {
		if (is_2d) {

			fread(&N, sizeof(uint32_t),1,pFile);
			if (N > 0) {
			showers::ChargeExcessProfile2D profile = read_shower(pFile, N);
			int find_shower = ce_profiles[profile.hadronic].count(profile.energy);
			if (find_shower == 0) {
				ce_profiles_2D[profile.hadronic][profile.energy] = std::vector<showers::ChargeExcessProfile2D>();
				stored_energies[profile.hadronic].push_back(profile.energy);
			}
			ce_profiles_2D[profile.hadronic][profile.energy].push_back(profile);		
			} else {
				reading_complete = 1;
			}
		} else {
			showers::ChargeExcessProfile profile;
			reading_complete = read_shower(
					pFile,
					&profile.N,
					&profile.hadronic,
					&profile.energy,
					&profile.grammage,
					&profile.charge_excess
					);
			int find_shower = ce_profiles[profile.hadronic].count(profile.energy);
			if (find_shower == 0) {
				ce_profiles[profile.hadronic][profile.energy] = std::vector<showers::ChargeExcessProfile>();
				stored_energies[profile.hadronic].push_back(profile.energy);
			}
			ce_profiles[profile.hadronic][profile.energy].push_back(profile);
		}
	}
	fclose(pFile);
};

showers::Shower1D showers::ShowerCreator::create_shower(
		std::array<float, 3> pos,
		double en,
		double zen,
		double az,
		int had
		) {
	double closest_energy = stored_energies[had][0];
	double energy_diff = std::abs(en - stored_energies[had][0]);
	for (int i=1; i < stored_energies[had].size(); i++) {
		if (std::abs(en - stored_energies[had][i]) < energy_diff) {
			closest_energy = stored_energies[had][i];
			energy_diff = std::abs(en - stored_energies[had][i]);
		}
	}
	std::default_random_engine generator{static_cast<long unsigned int>(time(NULL))};;
	std::uniform_int_distribution<int> dist(0, ce_profiles[had][closest_energy].size());
        int i_shower = 0;
        std::cout << "Pick shower " << i_shower << " out of " << ce_profiles[had][closest_energy].size() << " showers. \n";
	return showers::Shower1D(
			pos,
			en,
			zen,
			az,
			ce_profiles[had][closest_energy][i_shower],
			en / closest_energy,
			density_profile
	);
};

int showers::ShowerCreator::read_shower(FILE *f, int *N, int *hadronic, double *E, std::vector<double> *grammage, std::vector<double> *q)
{
	if (fread(N, sizeof(uint32_t),1,f) == 0) return 1;
	fread(hadronic, sizeof(uint32_t), 1, f);
	fread(E, sizeof(double),1,f);
	grammage -> resize(*N);
	q -> resize(*N);
	fread(&((*grammage)[0]), sizeof(double),*N,f);
	for (int i=0; i < (*grammage).size(); i++) {
		(*grammage)[i] = (*grammage)[i] * units::kilogram / units::square_meter;
		
	}
	fread(&((*q)[0]), sizeof(double), *N, f);
        
	return 0;
}
showers::ChargeExcessProfile2D showers::ShowerCreator::read_shower(
	FILE *f,
	unsigned int &N
) {
	std::vector<scalar_t> radius = {};
	radius.resize(10);
	for (int i=1;i<radius.size();i++) {
		radius[i] = radius[i - 1] + 1 * units::cm;
	}
	showers::ChargeExcessProfile2D profile({N, radius.size()});
	
	fread(&profile.hadronic, sizeof(uint32_t), 1, f);
	fread(&profile.energy, sizeof(scalar_t),1,f);
	profile.grammage.resize(N);

	scalar_t r_max = radius[radius.size() - 1];
	fread(&profile.grammage, sizeof(double),N,f);
	std::vector<double> charge_io = {};
	charge_io.resize(N);
	fread(&charge_io[0], sizeof(double), N, f);
	for (int i=0; i < profile.grammage.size(); i++) {
		profile.grammage[i] = profile.grammage[i] * units::kilogram / units::square_meter;
		for (int j=0; j < radius.size(); j ++) {
			profile.charge_excess(i, j) = charge_io[i] * (r_max - radius[j]);
		}
	}
	return profile;
}

