#include <chrono>
#include "GreensFunction.hh"
#include "Symmetry.hh"
#include "IteratorUtils.hh"
#include "Current.hh"
#include "Interpolation.hh"

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

// evaluates a toy Green's function in (t, r, z) coordinates
template <std::size_t dims, std::size_t vec_dims>
Vector<scalar_t, vec_dims> toy_green(const Vector<scalar_t, dims>& pos, const Vector<scalar_t, dims>& coeffs) {

  scalar_t lin_val = VectorUtils::inner(pos, coeffs);  
  Vector<scalar_t, vec_dims> retvec;
  for(std::size_t i = 0; i < vec_dims; i++) {
    retvec[i] = (i + 1) * lin_val * pos[0] * pos[0] / 100.0;
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

template <class CallableT>
scalar_t integrate_simple(scalar_t t_sig, const LineCurrentSegment& seg, CallableT&& green) {

  auto get_F_cart = [](const RZFieldVector& field, const XYZCoordVector& pos) -> XYZFieldVector {
    scalar_t x = pos.x();
    scalar_t y = pos.y();
    scalar_t r_xy = std::sqrt(pos.x() * pos.x() + pos.y() * pos.y());

    scalar_t cos_phi = 0.0, sin_phi = 0.0;
    if(r_xy > 0) {
      cos_phi = x / r_xy;
      sin_phi = y / r_xy;
    }
    
    return XYZFieldVector{
      field.r() * cos_phi,
      field.r() * sin_phi,
      field.z()
    };
  };
  
  scalar_t stepsize = 0.01f;
  
  XYZCoordVector seg_vel = (seg.end_pos - seg.start_pos) / (seg.end_time - seg.start_time);
  XYZCoordVector source_xyz = seg_vel * seg.charge;

  scalar_t sigval = 0.0;
  for(scalar_t t_p = seg.start_time; t_p < seg.end_time; t_p += stepsize) {
    XYZCoordVector cur_pos = seg.start_pos + seg_vel * (t_p - seg.start_time);
    XYZFieldVector green_val = get_F_cart(green(t_sig - t_p, cur_pos), cur_pos);

    for(std::size_t i = 0; i < 3; i++) {
      sigval += -green_val[i] * source_xyz[i];
    }
  }
  sigval *= stepsize;
  
  return sigval;
}

int main(int argc, char* argv[]) {

  (void)argc;
  (void)argv;
  
  std::filesystem::path workdir = "./darr_test";
  if(std::filesystem::exists(workdir)) {
    std::filesystem::remove_all(workdir);
  }

  Vector<scalar_t, dims> coeffs(1.0f);
  
  // prepare Green's function for testing purposes
  {
    Vector<std::size_t, dims> init_cache_el_shape(10);
    Vector<std::size_t, dims> streamer_chunk_shape(stor::INFTY);
    streamer_chunk_shape[0] = 1;
    
    darr_t darr(workdir, 5, init_cache_el_shape, streamer_chunk_shape);
    
    // fill distributed array with values
    Vector<std::size_t, dims> chunk_size(100);
    Vector<std::size_t, dims> start_ind(0);
    Vector<std::size_t, dims> end_ind(400);
    
    auto filler = [&](const Vector<std::size_t, dims>& ind){
      return toy_green<dims, vec_dims>(ind.template as_type<scalar_t>(), coeffs);
    };
    register_chunks(darr, start_ind, end_ind, chunk_size, filler);
    
    // transpose axes to go from (t, z, r) to (r, z, t)
    darr.SwapAxes<0, 2>();
    
    // reshape the chunks to make sure the overlap is respected for the interpolation
    
    std::filesystem::path workdir_tmp = "./darr_test_tmp";
    Vector<std::size_t, dims> requested_chunk_size(100);  
    darr.RebuildChunks(requested_chunk_size, workdir_tmp, 2, SpatialSymmetry::Cylindrical<T>::boundary_evaluator);
    
    auto val = darr[{250u, 250u, 250u}];
    for(auto cur : val) {
      std::cout << cur << std::endl;
    }  
    
    std::cout << darr.GetShape() << std::endl;
    
    RZTCoordVector start_pos{0.0f, 0.0f, 0.0f};
    RZTCoordVector end_pos{400.0f, 400.0f, 400.0f};
    RZTCoordVector sample_interval{1.0f, 1.0f, 1.0f};
    
    CylindricalGreensFunction(start_pos, end_pos, sample_interval, std::move(darr));
  }

  // read Green's function for evaluation purposes
  {
    CylindricalGreensFunction gf(workdir, 5);

    std::size_t num_samples = 150;
    
    std::vector<scalar_t> accum(num_samples, 0.0);
    // RZCoordVector cur_pt{123.0f, 123.0f};
    // RZCoordVector cur_pt{28.2843f, 20.0f};
    RZCoordVector cur_pt{28.2843f, 20.0f};
    // scalar_t t_start = 123.0f;
    scalar_t t_start = 0.001f;
    // scalar_t t_start = 16.801;
    scalar_t t_samp = 0.1;

    // RZFieldVector source{1.0f, 1.0f};
    RZFieldVector source{7.07107f, 0.0f};

    gf.accumulate_inner_product<Interpolation::Kernel::Keys>(cur_pt, t_start, t_samp, num_samples,
							     source, accum.begin());

    // test closure
    std::cout << "Testing closure ... ";
    for(std::size_t sample_ind = 0; sample_ind < num_samples; sample_ind++) {
      RZTCoordVector pos{t_start + t_samp * sample_ind, cur_pt.r(), cur_pt.z()};
      Vector<scalar_t, vec_dims> field = toy_green<dims, vec_dims>(pos, coeffs);
      
      scalar_t inner_product = 0.0;
      for(std::size_t i = 0; i < vec_dims; i++) {
	inner_product += field[i] * source[i];
      }

      if(std::fabs((inner_product - accum[sample_ind]) / inner_product) > 1e-3) {
	std::cout << "---------" << std::endl;
	std::cout << "Problem!" << std::endl;
	std::cout << "pos = " << pos << std::endl;
	std::cout << "field = " << field << std::endl;
	std::cout << "interpolated = " << accum[sample_ind] << " vs expected = " << inner_product << std::endl;
	std::cout << "---------" << std::endl;
      }
    }
    std::cout << "OK!" << std::endl;
  }
  
  // test some integration against a line current
  {
    CylindricalGreensFunction gf(workdir, 5);

    scalar_t t_sig_start = -100.0f;
    scalar_t t_sig_samp = 1.0f;
    std::size_t num_samples = 240;
    std::vector<scalar_t> signal_buffer(num_samples);

    // build current segment
    // scalar_t track_start_time = 120.0f;
    scalar_t track_start_time = -50.0f;
    // scalar_t track_end_time = 350.0f;
    scalar_t track_end_time = 50.0f;
    scalar_t track_charge = 1.0f;
    // LineCurrentSegment track(XYZCoordVector{10.0f, 10.0f, 20.0f},    // track start position
    // 			     XYZCoordVector{20.0f, 20.0f, 40.0f},   // track end position
    // 			     track_start_time, track_end_time, track_charge);

    LineCurrentSegment track(XYZCoordVector{-100.0f, 0.0f, 40.0f},    // track start position
			     XYZCoordVector{100.0f, 0.0f, 40.0f},   // track end position
			     track_start_time, track_end_time, track_charge);
    
    gf.apply_accumulate<Interpolation::Kernel::Keys>(track, t_sig_start, t_sig_samp, num_samples, signal_buffer);
    std::fill(signal_buffer.begin(), signal_buffer.end(), 0.0);
    
    auto start = std::chrono::high_resolution_clock::now();
    
    // integrate the current against the Green's function
    gf.apply_accumulate<Interpolation::Kernel::Keys>(track, t_sig_start, t_sig_samp, num_samples, signal_buffer);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);  

    auto green_for_eval = [&](scalar_t t, const XYZCoordVector& pos) {

      if(t < 0) {
	return Vector2D<scalar_t>{0.0f, 0.0f};
      }
      
      Vector<scalar_t, 3> pos_trz(0);
      pos_trz[0] = t;
      pos_trz[1] = std::sqrt(pos.x() * pos.x() + pos.y() * pos.y());
      pos_trz[2] = pos.z();
      
      return toy_green<3, 2>(pos_trz, coeffs);
    };
    
    std::cout << "BBBB got final signal BBBB" << std::endl;    
    for(std::size_t i = 0; i < signal_buffer.size(); i++) {

      scalar_t exp_sig = integrate_simple(t_sig_start + i * t_sig_samp, track, green_for_eval);
      
      std::cout << "t = " << t_sig_start + i * t_sig_samp << " , signal = " << signal_buffer[i] << ", expected = " << exp_sig << " deviation = " << (signal_buffer[i]-exp_sig)/exp_sig << std::endl;
    }
    std::cout << std::endl;
    std::cout << "took " << duration << std::endl;
  }
  
  std::cout << "done" << std::endl;  
}
