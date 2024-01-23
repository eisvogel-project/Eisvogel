#include "Eisvogel/DistributedWeightingField.hh"
#include "Eisvogel/Serialization.hh"

DistributedWeightingField::DistributedWeightingField(std::string wf_path, CoordVector start_coords, CoordVector end_coords) : m_wf_path(wf_path) {

  // Assumes that weighting field does not yet exist and directory structure needs to be built up from scratch
  
  // we might have to create the directory structure if this is an empty weighting field
  if(!std::filesystem::exists(m_wf_path)) {
    std::filesystem::create_directory(m_wf_path);
  }
  
  std::string wf_path_E_r = m_wf_path + "/E_r";
  std::string wf_path_E_z = m_wf_path + "/E_z";
  std::string wf_path_E_phi = m_wf_path + "/E_phi";
  
  m_E_r = std::make_shared<storage_t>(wf_path_E_r, 10);
  m_E_z = std::make_shared<storage_t>(wf_path_E_z, 10);
  m_E_phi = std::make_shared<storage_t>(wf_path_E_phi, 10);
  
  m_start_coords = std::make_shared<CoordVector>(start_coords);
  m_end_coords = std::make_shared<CoordVector>(end_coords);
}

DistributedWeightingField::DistributedWeightingField(std::string wf_path) : m_wf_path(wf_path) {
  
  // Assumes that weighting field already exists
  std::string wf_path_E_r = m_wf_path + "/E_r";
  std::string wf_path_E_z = m_wf_path + "/E_z";
  std::string wf_path_E_phi = m_wf_path + "/E_phi";
  
  m_E_r = std::make_shared<storage_t>(wf_path_E_r, 10);
  m_E_z = std::make_shared<storage_t>(wf_path_E_z, 10);
  m_E_phi = std::make_shared<storage_t>(wf_path_E_phi, 10);    
  
  // Read weighting field metadata
  std::fstream ifs;
  std::string meta_path = m_wf_path + "/wf_meta.bin";
  ifs.open(meta_path, std::ios::in | std::ios::binary);
  stor::Serializer iser(ifs);
  m_start_coords = std::make_shared<CoordVector>(iser.deserialize<CoordVector>());
  m_end_coords = std::make_shared<CoordVector>(iser.deserialize<CoordVector>());
}

DistributedWeightingField::~DistributedWeightingField() {
  Flush();
}

void DistributedWeightingField::RegisterChunk(const chunk_t& chunk_E_r, const chunk_t& chunk_E_z, const chunk_t& chunk_E_phi,
					      const IndexVector start_ind) {

  m_E_r -> RegisterChunk(chunk_E_r, start_ind);
  m_E_z -> RegisterChunk(chunk_E_z, start_ind);
  m_E_phi -> RegisterChunk(chunk_E_phi, start_ind);
}

void DistributedWeightingField::Flush() {

  // Flush weighting field metadata
  std::fstream ofs;
  std::string meta_path = m_wf_path + "/wf_meta.bin";
  ofs.open(meta_path, std::ios::out | std::ios::binary);  
  stor::Serializer oser(ofs);
  oser.serialize(*m_start_coords);
  oser.serialize(*m_end_coords);
  
  // Flush the actual data
  m_E_r -> Flush();
  m_E_z -> Flush();
  m_E_phi -> Flush();
}
