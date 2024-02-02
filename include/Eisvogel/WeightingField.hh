#ifndef __WEIGHTING_FIELD_HH
#define __WEIGHTING_FIELD_HH

#include <memory>
#include <string>
#include <filesystem>

#include "Common.hh"
#include "CoordUtils.hh"
#include "Interpolator.hh"
#include "Kernels.hh"

#include "FieldStorage.hh"

template <typename SymmetryT>
class TRZFieldIndexer {

public:

  TRZFieldIndexer(std::filesystem::path index_path);

  void SetShape(IndexVector shape); // temporary, to be removed
  
  CoordVector GetFieldIndexFromCoord(CoordVector pos_txyz);
  IndexVector GetFieldIndex(GridVector inds_trz);

private:
  std::shared_ptr<IndexVector> m_shape;
  std::shared_ptr<CoordVector> m_start_coords;
  std::shared_ptr<CoordVector> m_end_coords;  
};

template <class FieldIndexerT, class FieldStorageT>
class WeightingField {

public:

  WeightingField(std::string wf_path);
  
  template <typename KernelT = KeysCubicInterpolationKernelNew>
  scalar_t Er(CoordVector pos) {    
    return eval<KernelT>(pos, [&](auto& arg){return m_field_storage -> E_r(arg);});
  };

private:
  
  template <typename KernelT, typename GetterT,
	    typename ValueT = std::invoke_result_t<GetterT, IndexVector&>>
  ValueT eval(CoordVector pos, GetterT getter) {

    // Weighting fields are nonzero only for t > 0
    if(CoordUtils::getT(pos) <= 0) {
      return 0.0;
    }

    CoordVector frac_inds = m_field_indexer -> GetFieldIndexFromCoord(pos);

    auto to_interpolate = [&](GridVector& vec) -> ValueT {
      IndexVector inds = m_field_indexer -> GetFieldIndex(vec);
      return getter(inds);
    };

    return InterpolateFuncNew<KernelT>(to_interpolate, frac_inds);        
  }
  
private:

  std::string m_wf_path;
  
  std::shared_ptr<FieldStorageT> m_field_storage;
  std::shared_ptr<FieldIndexerT> m_field_indexer;  
};

#include "WeightingField.hxx"

#endif
