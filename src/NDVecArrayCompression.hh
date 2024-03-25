#pragma once

#include "Eisvogel/MathUtils.hh"

namespace nullsup {

  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  std::size_t calculate_required_buflen(const ArrayT<T, dims, vec_dims>& arr) {  
    // compute worst-case buffer size for a given array
    return arr.GetNumberElements() * vec_dims + MathUtils::IntDivCeil(arr.GetNumberElements(), 2);
  }

  template <std::size_t dims>
  std::size_t calculate_required_buflen(const Vector<std::size_t, dims>& shape, std::size_t vec_dims) {
    std::size_t number_elements = 1;
    for(std::size_t ind = 0; ind < dims; ind++) {
      number_elements *= shape[ind];
    }
    return number_elements * vec_dims + MathUtils::IntDivCeil(number_elements, 2);
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    std::size_t dims, std::size_t vec_dims>
  std::size_t suppress_zero(const ArrayT<float, dims, vec_dims>& arr, std::span<uint32_t>&& buffer) {  
    auto postprocessor = [](const uint32_t& host) -> uint32_t {
      return htonl(host);
    };  
    return suppress_zero(arr, std::move(buffer), postprocessor);
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    std::size_t dims, std::size_t vec_dims>
  std::size_t desuppress_zero(const std::span<uint32_t>&& buffer, ArrayT<float, dims, vec_dims>& arr) {
    auto preprocessor = [](const uint32_t& network) -> uint32_t {
      return ntohl(network);
    };
    return desuppress_zero(std::move(buffer), arr, preprocessor);
  }
  
  // ---------
  
  // returns elements in buffer that need to be considered
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims,
	    typename SerType, class CallableT>
  std::size_t suppress_zero(const ArrayT<T, dims, vec_dims>& arr, std::span<SerType>&& buffer, CallableT&& postprocessor) {
    
    SerType num_zeros = 0;
    auto buffer_it = buffer.begin();
    
    auto zero_suppressor = [&](Vector<std::size_t, dims>& ind) {
      
      if(arr.IsZero(ind)) {      
	// if compression is triggered, assume that this array has may zeroes
	[[likely]];
	
	if(num_zeros == 0) {
	  [[unlikely]];
	  
	  // first zero, put into buffer
	  std::fill_n(std::execution::unseq, buffer_it, vec_dims, postprocessor(0));
	  buffer_it += vec_dims;
	}
	
	// keep counting the streak of zeros
	num_zeros++;
      }
      else {
	if(num_zeros > 0) {	
	  // first non-zero after a streak of zeros: write number count of zeros
	  std::fill_n(std::execution::unseq, buffer_it, 1, postprocessor(num_zeros));
	  buffer_it++;
	  
	  num_zeros = 0;
	}
	
	// copy the array element
	for(float cur_val: arr[ind]) {
	  SerType ser_val = reinterpret_cast<const SerType&>(cur_val);
	  *buffer_it = postprocessor(ser_val);
	  buffer_it++;
	}
      }
    };
    
    loop_over_array_elements(arr, zero_suppressor);
    
    // close any remaining open zeroes
    if(num_zeros > 0) {
      std::fill_n(std::execution::unseq, buffer_it, 1, htonl(num_zeros));
      buffer_it += 1;
    }
    
    return buffer_it - buffer.begin();
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims,
	    typename SerType, class CallableT>
  std::size_t desuppress_zero(const std::span<SerType>&& buffer, ArrayT<T, dims, vec_dims>& arr, CallableT&& preprocessor) {
    
    SerType num_zeros = 0;
    Vector<T, vec_dims> vec_buffer;
    
    auto buffer_it = buffer.begin();
    
    auto zero_desuppressor = [&](Vector<std::size_t, dims>& ind) {
      
      if(num_zeros > 0) {
	[[likely]];
	
	// Fill zero elements into array
	std::fill_n(std::execution::unseq, arr[ind].begin(), vec_dims, 0);
	num_zeros--;
      }
      else {
	
	// Read and fill next element
	for(std::size_t ind = 0; ind < vec_dims; ind++) {	
	  SerType ser_val = preprocessor(*buffer_it);
	  buffer_it++;	
	  std::memcpy(&vec_buffer[ind], &ser_val, sizeof(ser_val));
	}
	
	arr[ind] = vec_buffer;
	
	if(arr.IsZero(ind)) {
	  num_zeros = preprocessor(*buffer_it) - 1; // the first zero has already been written into the array
	  buffer_it++;
	}
      }    
    };
    
    loop_over_array_elements(arr, zero_desuppressor);
    
    return buffer_it - buffer.begin();
  }
}
