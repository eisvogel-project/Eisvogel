namespace SpatialSymmetry {

  template <typename T>
  void Cylindrical<T>::boundary_evaluator(darr_t& darr, const RZTSignedIndexVector& ind, view_t elem) {
    
    RZTIndexVector darr_shape = darr.GetShape();   
    
    const Vector<T, vec_dims> zero_val(0.0);

    // Put zero field values everywhere outside the positive array boundary
    if(std::cmp_greater_equal(ind.t(), darr_shape.t()) ||
       std::cmp_greater_equal(ind.r(), darr_shape.r()) ||
       std::cmp_greater_equal(ind.z(), darr_shape.z())) {
      elem = zero_val;
      return;
    }

    // There is no field also outside the negative array boundary in the t and z directions
    if((ind.t() < 0) ||
       (ind.z() < 0)) {
      elem = zero_val;
      return;
    }

    // for r < 0, take the value at the corresponding |r| > 0 coordinate
    if((ind.r() < 0)) {
      RZTSignedIndexVector ind_wrapped = ind;
      ind_wrapped.r() = std::abs(ind_wrapped.r());      
      RZTIndexVector ind_to_fetch = ind_wrapped.template as_type<std::size_t>();      
      elem = darr[ind_to_fetch];
      return;
    }    
    throw std::logic_error("Control flow should never reach here.");
  }  
}
