#ifndef __DISTRIBUTED_WEIGHTING_FIELD_HH
#define __DISTRIBUTED_WEIGHTING_FIELD_HH

#include <string>
#include <filesystem>
#include <memory>
#include "DistributedNDArray.hh"
#include "Common.hh"
#include "CoordUtils.hh"

class DistributedWeightingField {

public:
  using storage_t = DistributedNDArray<scalar_t, 3>;

public:

  using index_t = std::size_t;
  using shape_t = typename storage_t::shape_t;

  DistributedWeightingField(std::string wf_path, CoordVector start_coords, CoordVector end_coords);
  DistributedWeightingField(std::string wf_path);
  ~DistributedWeightingField();

  void Flush();
  
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
