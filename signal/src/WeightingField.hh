#ifndef __WEIGHTING_FIELD_HH
#define __WEIGHTING_FIELD_HH

#include <type_traits>
#include "NDArray.hh"
#include "Common.hh"
#include "CoordUtils.hh"

#include "Serialization.hh"

class WeightingField {

public:
  using storage_t = DenseNDArray<scalar_t, 3>;

private:
  friend struct stor::Traits<WeightingField>;

public:

  using index_t = std::size_t;
  using shape_t = typename storage_t::shape_t;
 
  WeightingField(storage_t&& E_r, storage_t&& E_z, storage_t&& E_phi,
		 CoordVector start_coords, CoordVector end_coords) :
    m_shape(E_r.shape()), m_E_r(E_r), m_E_z(E_z), m_E_phi(E_phi), m_start_coords(start_coords), m_end_coords(end_coords) {

    if(!(m_E_r.shape() == m_E_z.shape()) && (m_E_z.shape() == m_E_phi.shape())) {
      throw;
    }
  };

  scalar_t E_r(index_t ind_t, index_t ind_r, index_t ind_z) {
    return m_E_r(ind_t, ind_r, ind_z);
  };

  const storage_t& E_r() const {return m_E_r;};
  const storage_t& E_z() const {return m_E_z;};
  const storage_t& E_phi() const {return m_E_phi;};

  scalar_t E_z(index_t ind_t, index_t ind_r, index_t ind_z) {
    return m_E_z(ind_t, ind_r, ind_z);
  };

  scalar_t E_phi(index_t ind_t, index_t ind_r, index_t ind_z) {
    return m_E_phi(ind_t, ind_r, ind_z);
  };

  static inline CoordVector FracIndsToCoord(const CoordVector& frac_inds, const CoordVector& start_coords, const CoordVector& end_coords, 
  					const CoordVector& shape) {
    return start_coords + (end_coords - start_coords) / shape * frac_inds;
  }

  static inline CoordVector CoordToFracInds(const CoordVector& coords, const CoordVector& start_coords, 
					    const CoordVector& end_coords, const CoordVector& shape) {
    return (coords - start_coords) / (end_coords - start_coords) * shape;
  }

  CoordVector getFracInds(const CoordVector& coords) const {

    std::cout << "coords = ";
    for(auto cur: coords) {
      std::cout << cur << "  ";
    }
    std::cout << std::endl;

    std::cout << "m_start_coords = ";
    for(auto cur: m_start_coords) {
      std::cout << cur << "  ";
    }
    std::cout << std::endl;

    std::cout << "m_end_coords = ";
    for(auto cur: m_end_coords) {
      std::cout << cur << "  ";
    }
    std::cout << std::endl;

    std::cout << "m_shape = ";
    for(auto cur: m_shape) {
      std::cout << cur << "  ";
    }
    std::cout << std::endl;

    CoordVector frac_inds = CoordToFracInds(coords, m_start_coords, m_end_coords, m_shape);

    std::cout << "frac_inds = ";
    for(auto cur: frac_inds) {
      std::cout << cur << "  ";
    }
    std::cout << std::endl;

    return frac_inds;
  }

  DeltaVector getSamplingIntervals() const {
    return (m_end_coords - m_start_coords) / m_shape;
  }

private:

  CoordVector m_shape;

  storage_t m_E_r;
  storage_t m_E_z;
  storage_t m_E_phi;

  CoordVector m_start_coords;
  CoordVector m_end_coords;
};

namespace stor {
  
  template <>
  struct Traits<WeightingField> {
    using type = WeightingField;
    using storage_t = typename WeightingField::storage_t;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<storage_t>::serialize(stream, val.m_E_r);
      Traits<storage_t>::serialize(stream, val.m_E_z);
      Traits<storage_t>::serialize(stream, val.m_E_phi);
      Traits<CoordVector>::serialize(stream, val.m_start_coords);
      Traits<CoordVector>::serialize(stream, val.m_end_coords);
    }
    
    static type deserialize(std::iostream& stream) {
      storage_t E_r = Traits<storage_t>::deserialize(stream);
      storage_t E_z = Traits<storage_t>::deserialize(stream);
      storage_t E_phi = Traits<storage_t>::deserialize(stream);
      CoordVector start_coords = Traits<CoordVector>::deserialize(stream);
      CoordVector end_coords = Traits<CoordVector>::deserialize(stream);
      return WeightingField(std::move(E_r), std::move(E_z), std::move(E_phi),
			    start_coords, end_coords);
    }
  };
}

#endif
