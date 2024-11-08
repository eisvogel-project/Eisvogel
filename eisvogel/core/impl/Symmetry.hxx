namespace SpatialSymmetry {

  template <typename T>
  void Cylindrical<T>::boundary_evaluator(chunk_t& chunk, const RZTSignedIndexVector& chunk_start_ind, const RZTSignedIndexVector& chunk_end_ind,
					  const RZTSignedIndexVector& elem_ind, view_t elem) {

    (void)chunk_end_ind;
    
    const Vector<T, vec_dims> zero_val(0.0);

    // For r < 0, take the value at the corresponding |r| > 0 coordinate
    if(elem_ind.r() < 0) {

      RZTSignedIndexVector ind_wrapped = elem_ind;
      ind_wrapped.r() = std::abs(ind_wrapped.r());
      RZTIndexVector chunk_ind_wrapped = (ind_wrapped - chunk_start_ind).template as_type<std::size_t>();

      assert(chunk.index_within_bounds(chunk_ind_wrapped));  // Make sure we're not accessing out-of-bounds
      elem = chunk[chunk_ind_wrapped];
      return;
    }

    // If this is not an element at the r = 0 boundary, handling it is simpler: there is no field in this direction
    elem = zero_val;
  }
}
