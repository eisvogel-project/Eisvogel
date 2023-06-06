#include "Eisvogel/Common.hh"
#include "shower_creator.h"
#include "shower_1D.h"
#include "shower_2D.h"
#include "charge_excess_profile.cpp"
#include "units.h"
#include <iostream>
#include <array>
#include <stdio.h>
#include <time.h>
#include <stdexcept>
#include <random>


showers::ShowerCreator::ShowerCreator(
		std::string file_path
		) {
	shower_file = file_path;
	density_profile = environment::IceProfile();
	int reading_complete = 0;
	int N;
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
		int had,
		int i_shower
		) {
	double closest_energy = stored_energies[had][0];
	double energy_diff = std::abs(en - stored_energies[had][0]);
	for (int i=1; i < stored_energies[had].size(); i++) {
		if (std::abs(en - stored_energies[had][i]) < energy_diff) {
			closest_energy = stored_energies[had][i];
			energy_diff = std::abs(en - stored_energies[had][i]);
		}
	}
	int shower_index;
	if (i_shower < 0) {
		std::random_device rand_dev;
		std::default_random_engine generator(rand_dev());
		std::uniform_int_distribution<int> dist(0, ce_profiles[had][closest_energy].size());
		shower_index = dist(generator);
	} else {
		shower_index = i_shower;
	}
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


showers::ShowerCreator2D::ShowerCreator2D(
	std::string file_path
) {
	shower_file = file_path;
	density_profile = environment::IceProfile();
	int reading_complete = 0;
	int N;
	FILE * pFile = fopen(shower_file.c_str(), "rb");
	ce_profiles[0] = std::map<double, std::vector<showers::ChargeExcessProfile2D>>();
	ce_profiles[1] = std::map<double, std::vector<showers::ChargeExcessProfile2D>>();
	
	while (reading_complete == 0) {
		if(fread(&N, sizeof(uint32_t), 1, pFile)==0) {
			reading_complete = 1;
		} else {
    
		showers::ChargeExcessProfile2D profile = read_shower(
				pFile,
				&N
				);
		
		int find_shower = ce_profiles[profile.hadronic].count(profile.energy);
		if (find_shower == 0) {
			ce_profiles[profile.hadronic][profile.energy] = std::vector<showers::ChargeExcessProfile2D>();
			stored_energies[profile.hadronic].push_back(profile.energy);
		}
		ce_profiles[profile.hadronic][profile.energy].push_back(profile);
		}
	}
	
	fclose(pFile);
}

showers::ChargeExcessProfile2D showers::ShowerCreator2D::read_shower(
	FILE *f,
	int *N
	) {
    showers::ChargeExcessProfile2D profile({*N, 10});
    profile.radius.resize(10);  // Resize the radius vector based on the read size

    for (int i = 1; i < profile.radius.size(); i++) {
        profile.radius[i] = profile.radius[i - 1] + 2 * units::cm;
    }


    fread(&profile.hadronic, sizeof(uint32_t), 1, f);
    fread(&profile.energy, sizeof(double), 1, f);
    profile.grammage.resize(*N);

    scalar_t r_max = profile.radius[profile.radius.size() - 1];

    std::vector<double> charge_io;
    charge_io.resize(*N);  // Resize the charge_io vector based on the read size

    fread(&(profile.grammage[0]), sizeof(double), *N, f);
    fread(&(charge_io[0]), sizeof(double), *N, f);

    for (int i = 0; i < profile.grammage.size(); i++) {
        profile.grammage[i] = profile.grammage[i] * units::kilogram / units::square_meter;
        for (int j=0; j < profile.radius.size(); j ++) {
            profile.charge_excess(i, j) = charge_io[i] * (r_max - profile.radius[j]);
        }
    }
    return profile;
}


showers::Shower2D showers::ShowerCreator2D::create_shower(
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
        int i_shower = 4;
        std::cout << "Pick shower " << i_shower << " out of " << ce_profiles[had][closest_energy].size() << " showers. \n";
	return showers::Shower2D(
			pos,
			en,
			zen,
			az,
			ce_profiles[had][closest_energy][i_shower],
			en / closest_energy,
			density_profile
	);
};