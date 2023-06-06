#ifndef __SPARSE_CURRENT_DENSITY_3D_HH
#define __SPARSE_CURRENT_DENSITY_3D_HH

#include <vector>
#include <utility>
#include "CoordUtils.hh"

namespace CU = CoordUtils;

using CurrentElement = std::pair<CoordVector, FieldVector>;

class SparseCurrentDensity3D {

public:
  SparseCurrentDensity3D(const DeltaVector& voxel_size) : m_voxel_size(voxel_size), m_elements({}) { };
  SparseCurrentDensity3D(const DeltaVector& voxel_size, const std::vector<CurrentElement>&& elements) : 
    m_voxel_size(voxel_size), m_elements(elements) { };

  void addCurrentElement(const CoordVector& point, const FieldVector& current_density) {
    addCurrentElement({point, current_density});
  }

  void addCurrentElement(const CurrentElement& element) {
    m_elements.push_back(element);
  }

  scalar_t getVolumeElementTXYZ() const {
    return CU::getT(m_voxel_size) * getVolumeElementXYZ();
  }

  scalar_t getVolumeElementXYZ() const {
    return CU::getX(m_voxel_size) * CU::getY(m_voxel_size) * CU::getZ(m_voxel_size);
  }

  int getNumberofElements() {return m_elements.size();}

  auto begin() const {return m_elements.cbegin();}
  auto end() const {return m_elements.cend();}

private:

  DeltaVector m_voxel_size;
  std::vector<CurrentElement> m_elements;

};

#endif
