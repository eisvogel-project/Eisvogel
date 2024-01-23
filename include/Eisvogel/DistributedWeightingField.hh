#ifndef __DISTRIBUTED_WEIGHTING_FIELD_HH
#define __DISTRIBUTED_WEIGHTING_FIELD_HH

#include <string>
#include <filesystem>
#include <memory>
#include "DistributedNDArray.hh"
#include "Common.hh"
#include "CoordUtils.hh"

#include "Serialization.hh"

class DistributedWeightingField {

public:
  using storage_t = DistributedNDArray<scalar_t, 3>;

public:

  using index_t = std::size_t;
  using shape_t = typename storage_t::shape_t;

  DistributedWeightingField(std::string wf_path, CoordVector start_coords, CoordVector end_coords) : m_wf_path(wf_path) {

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

  DistributedWeightingField(std::string wf_path) : m_wf_path(wf_path) {

    // Assumes that weighting field already exists
    std::string wf_path_E_r = m_wf_path + "/E_r";
    std::string wf_path_E_z = m_wf_path + "/E_z";
    std::string wf_path_E_phi = m_wf_path + "/E_phi";

    m_E_r = std::make_shared<storage_t>(wf_path_E_r, 10);
    m_E_z = std::make_shared<storage_t>(wf_path_E_z, 10);
    m_E_phi = std::make_shared<storage_t>(wf_path_E_phi, 10);    
    
    
  }
  
  const storage_t& E_r() const {return *m_E_r;};
  const storage_t& E_z() const {return *m_E_z;};
  const storage_t& E_phi() const {return *m_E_phi;};

  static inline CoordVector FracIndsToCoord(const CoordVector& frac_inds, const CoordVector& start_coords, const CoordVector& end_coords, 
					    const CoordVector& shape) {
    return start_coords + (end_coords - start_coords) / shape * frac_inds;
  }

  static inline CoordVector CoordToFracInds(const CoordVector& coords, const CoordVector& start_coords, 
					    const CoordVector& end_coords, const CoordVector& shape) {
    return (coords - start_coords) / (end_coords - start_coords) * shape;
  }

  CoordVector getFracInds(const CoordVector& coords) const {
    return CoordToFracInds(coords, *m_start_coords, *m_end_coords, m_E_r -> shape());
  }

  DeltaVector getSamplingIntervals() const {
    return (*m_end_coords - *m_start_coords) / m_E_r -> shape();
  }

private:

  std::string m_wf_path; 

  std::shared_ptr<storage_t> m_E_r;
  std::shared_ptr<storage_t> m_E_z;
  std::shared_ptr<storage_t> m_E_phi;

  std::shared_ptr<CoordVector> m_start_coords;
  std::shared_ptr<CoordVector> m_end_coords;
};

#endif
