#include "shower_creator.h"
#include "shower_1D.h"
#include "charge_excess_profile.cpp"
#include <iostream>
#include <array>
#include <stdio.h>
#include <stdexcept>
#include <random>


showers::ShowerCreator::ShowerCreator(
		std::string file_path
		) {
	shower_file = file_path;
	density_profile = environment::IceProfile();
	int reading_complete = 0;
	FILE * pFile = fopen(shower_file.c_str(), "rb");
	ce_profiles[0] = std::map<double, std::vector<showers::ChargeExcessProfile>>();
	ce_profiles[1] = std::map<double, std::vector<showers::ChargeExcessProfile>>();
	while (reading_complete == 0) {
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
	std::default_random_engine generator;
	std::uniform_int_distribution<int> dist(0, ce_profiles[had][closest_energy].size());

	return showers::Shower1D(
			pos,
			en,
			zen,
			az,
			ce_profiles[had][closest_energy][dist(generator)],
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
	fread(&((*q)[0]), sizeof(double), *N, f);
	return 0;
}

