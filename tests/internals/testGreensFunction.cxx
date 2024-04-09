#include "GreensFunction.hh"
#include "Symmetry.hh"
#include "Eisvogel/IteratorUtils.hh"

using T = float;
constexpr std::size_t dims = 3;
constexpr std::size_t vec_dims = 2;
using darr_t = DistributedNDVecArray<NDVecArray, T, dims, vec_dims>;
using view_t = typename DistributedNDVecArray<NDVecArray, T, dims, vec_dims>::view_t;

template <std::size_t dims, std::size_t vec_dims, class CallableT>
void fill_array(NDVecArray<T, dims, vec_dims>& to_fill, const Vector<std::size_t, dims>& start_ind, CallableT&& fill_function) {

  auto worker = [&](const Vector<std::size_t, dims>& ind) {
    to_fill[ind] = fill_function(start_ind + ind);
  };  
  IteratorUtils::index_loop_over_array_elements(to_fill, worker);
}

// evaluates a linear function
template <std::size_t dims, std::size_t vec_dims>
Vector<scalar_t, vec_dims> linear(const Vector<scalar_t, dims>& pos, const Vector<scalar_t, dims>& coeffs) {

  scalar_t lin_val = VectorUtils::inner(pos, coeffs);  
  Vector<scalar_t, vec_dims> retvec;
  for(std::size_t i = 0; i < vec_dims; i++) {
    retvec[i] = (i + 1) * lin_val;
  }
  
  return retvec;
}

template <std::size_t dims, std::size_t vec_dims, class CallableT>
void register_chunks(DistributedNDVecArray<NDVecArray, T, dims, vec_dims>& darr, const Vector<std::size_t, dims>& start_ind, const Vector<std::size_t, dims>& end_ind,
		     const Vector<std::size_t, dims>& chunk_shape, CallableT&& filler) {

  auto chunk_registerer = [&](const Vector<std::size_t, dims>& chunk_start, const Vector<std::size_t, dims>& chunk_end) {

    // prepare filled array
    NDVecArray<T, dims, vec_dims> array_buffer(chunk_shape);    
    fill_array(array_buffer, chunk_start, filler);

    // register as chunk
    std::cout << "Registering chunk: " << chunk_start << " -> " << chunk_end << std::endl;
    darr.RegisterChunk(array_buffer, chunk_start);
  };
  
  IteratorUtils::index_loop_over_chunks(start_ind, end_ind, chunk_shape, chunk_registerer);  
}

int main(int argc, char* argv[]) {
  
  std::filesystem::path workdir = "./darr_test";
  if(std::filesystem::exists(workdir)) {
    std::filesystem::remove_all(workdir);
  }

  // prepare empty distributed array
  Vector<std::size_t, dims> init_cache_el_shape(10);
  Vector<std::size_t, dims> streamer_chunk_shape(stor::INFTY);
  streamer_chunk_shape[0] = 1;
  
  darr_t darr(workdir, 5, init_cache_el_shape, streamer_chunk_shape);

  // fill distributed array with values
  Vector<std::size_t, dims> chunk_size(100);
  Vector<std::size_t, dims> start_ind(0);
  Vector<std::size_t, dims> end_ind(400);
  
  Vector<scalar_t, dims> coeffs(1.0f);
  auto filler = [&](const Vector<std::size_t, dims>& ind){
    return linear<dims, vec_dims>(ind.template as_type<scalar_t>(), coeffs);
  };
  register_chunks(darr, start_ind, end_ind, chunk_size, filler);

  // transpose axes to go from (t, z, r) to (r, z, t)
  darr.SwapAxes<0, 2>();
  
  // reshape the chunks to make sure the overlap is respected for the interpolation
  
  std::filesystem::path workdir_tmp = "./darr_test_tmp";
  Vector<std::size_t, dims> requested_chunk_size(100);  
  darr.RebuildChunks(requested_chunk_size, workdir_tmp, 2, SpatialSymmetry::Cylindrical<T, vec_dims>::boundary_evaluator);
  
  auto val = darr[{250u, 250u, 250u}];
  for(auto cur : val) {
    std::cout << cur << std::endl;
  }  
  
  RZTCoordVector start_pos{0.0f, -400.0f, 0.0f};
  RZTCoordVector end_pos{100.0f, 400.0f, 100.0f};
  RZTCoordVector sample_interval{1.0f, 1.0f, 1.0f};

  std::cout << start_pos.r() << std::endl;
  
  std::cout << "done" << std::endl;  
}
