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
  using chunk_t = DenseNDArray<scalar_t, 3>;

  using index_t = std::size_t;
  using shape_t = typename storage_t::shape_t;

  DistributedWeightingField(std::string wf_path, CoordVector start_coords, CoordVector end_coords);
  DistributedWeightingField(std::string wf_path);
  ~DistributedWeightingField();

  void Flush();

  void RegisterChunk(const chunk_t& chunk_E_r, const chunk_t& chunk_E_z, const chunk_t& chunk_E_phi,
		     const IndexVector chunk_start_inds);

  void MakeIndexPersistent();
  
  storage_t& E_r() const {return *m_E_r;};
  storage_t& E_z() const {return *m_E_z;};
  storage_t& E_phi() const {return *m_E_phi;};

  scalar_t E_r(IndexVector& ind) const {return m_E_r -> operator()(ind);};
  scalar_t E_z(IndexVector& ind) const {return m_E_z -> operator()(ind);};
  scalar_t E_phi(IndexVector& ind) const {return m_E_phi -> operator()(ind);};

  std::size_t shape(std::size_t dim) const {return m_E_r -> shape(dim);};
  std::size_t startInd(std::size_t dim) const {return m_E_r -> startInd(dim);};
  
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

  // TODO: this will be replaced by a single array of vectors (instead of scalars) asap, to keep the different
  // vector components local in memory
  std::shared_ptr<storage_t> m_E_r;
  std::shared_ptr<storage_t> m_E_z;
  std::shared_ptr<storage_t> m_E_phi;

  std::shared_ptr<CoordVector> m_start_coords;
  std::shared_ptr<CoordVector> m_end_coords;
};

#endif
