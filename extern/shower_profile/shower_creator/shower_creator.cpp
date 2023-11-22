#include "Eisvogel/Common.hh"
#include "Eisvogel/NDArray.hh"
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
#include "H5Cpp.h"


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
	ce_profiles.push_back(
		read_shower(file_path)
	);

}

showers::ChargeExcessProfile2D showers::ShowerCreator2D::read_shower(
	std::string file_path
	) {
	H5::H5File this_file(file_path, H5F_ACC_RDONLY);
    H5::DataSet this_dataset = this_file.openDataSet("ce");
	H5T_class_t type_class = this_dataset.getTypeClass();
	H5::DataSpace dataspace = this_dataset.getSpace();
	int rank = dataspace.getSimpleExtentNdims();
	hsize_t dimensions_out[2];
	int dim = dataspace.getSimpleExtentDims(dimensions_out, NULL);
	hsize_t start[2] = {0, 0};
	dataspace.selectHyperslab(H5S_SELECT_SET, dimensions_out, start);
	H5::DataSpace memspace(rank, dimensions_out);
	memspace.selectHyperslab(H5S_SELECT_SET, dimensions_out, start);
	double data_out[dimensions_out[0]][dimensions_out[1]] = {0};
	this_dataset.read( data_out, H5::PredType::NATIVE_DOUBLE, memspace, dataspace );	
	showers::ChargeExcessProfile2D profile({dimensions_out[0], dimensions_out[1]});
	H5::DataSet grammage_dataset = this_file.openDataSet("depth");
	H5::DataSet radius_dataset = this_file.openDataSet("radius");
	profile.grammage = readDataSet(&grammage_dataset);
	profile.radius = readDataSet(&radius_dataset);
    for(int i=0; i < dimensions_out[0]; i++) {
		for (int j=0; j < dimensions_out[1]; j++) {
			profile.set_charge_excess(i, j, data_out[i][j]);
		}
	}
	return profile;
}

std::vector<double> showers::ShowerCreator2D::readDataSet(
	H5::DataSet *dataset
) {
	H5T_class_t type_class = dataset->getTypeClass();
	H5::DataSpace dataspace = dataset->getSpace();
	int rank = dataspace.getSimpleExtentNdims();
	hsize_t dimensions_out[rank];
	int dim = dataspace.getSimpleExtentDims(dimensions_out, NULL);
	hsize_t start[rank] = {0};
	dataspace.selectHyperslab(H5S_SELECT_SET, dimensions_out, start);
	H5::DataSpace memspace(rank, dimensions_out);
	memspace.selectHyperslab(H5S_SELECT_SET, dimensions_out, start);
	double data_out[dimensions_out[0]];
	dataset->read(data_out, H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
	std::vector<double> vec_out (data_out, data_out + dimensions_out[0]);
	return vec_out;
}

showers::Shower2D showers::ShowerCreator2D::create_shower(
	std::array<float, 3> pos,
		double en,
		double zen,
		double az,
		int had
) {

	return showers::Shower2D(
			pos,
			en,
			zen,
			az,
			ce_profiles[0], //[had][closest_energy][i_shower],
			1, //en / closest_energy,
			density_profile
	);
};
