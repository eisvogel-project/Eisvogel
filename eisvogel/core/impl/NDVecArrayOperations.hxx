#include "IteratorUtils.hh"

namespace {

  template <std::size_t dims>
  Vector<std::size_t, dims> get_first_downsampled_ind_in_range(
      const Vector<std::size_t, dims>& start_ind,
      const Vector<std::size_t, dims>& downsampling)
  {
    return (start_ind + downsampling - 1) / downsampling;
  }

  template <std::size_t dims>
  Vector<std::size_t, dims> get_last_downsampled_ind_in_range(
      const Vector<std::size_t, dims>& end_ind,
      const Vector<std::size_t, dims>& downsampling)
  {
    return (end_ind - 1) / downsampling + 1; // upper range ends are exclusive
  }

  template <std::size_t dims>
  Vector<std::size_t, dims> get_shape(
      const Vector<std::size_t, dims>& downsampled_start_ind,
      const Vector<std::size_t, dims>& downsampled_end_ind)
  {
    return VectorUtils::max(downsampled_end_ind, downsampled_start_ind) -
           downsampled_start_ind;
  }
} // namespace

namespace Downsampling {

  template <std::size_t dims>
  Vector<std::size_t, dims> get_downsampled_shape(
      const Vector<std::size_t, dims>& start_ind,
      const Vector<std::size_t, dims>& shape,
      const Vector<std::size_t, dims>& downsampling)
  {
    return get_shape(
        get_first_downsampled_ind_in_range(start_ind, downsampling),
        get_last_downsampled_ind_in_range(start_ind + shape, downsampling));
  }

  template <typename T, std::size_t dims, std::size_t vec_dims>
  void downsample(const NDVecArray<T, dims, vec_dims>& to_downsample,
                  const Vector<std::size_t, dims> start_ind,
                  const Vector<std::size_t, dims> downsampling,
                  NDVecArray<T, dims, vec_dims>& downsampled,
                  Vector<std::size_t, dims>& downsampled_start_ind)
  {

    using ind_t = typename NDVecArray<T, dims, vec_dims>::ind_t;
    using shape_t = typename NDVecArray<T, dims, vec_dims>::shape_t;

    ind_t end_ind = start_ind + to_downsample.GetShape();

    // The start and end indices on the downsampled grid that are covered by the
    // input array
    downsampled_start_ind =
        get_first_downsampled_ind_in_range(start_ind, downsampling);
    ind_t downsampled_end_ind =
        get_last_downsampled_ind_in_range(end_ind, downsampling);

    // Prepare output buffer with the correct shape (keeping in mind that the
    // range might be empty along some directions)
    shape_t downsampled_shape =
        get_shape(downsampled_start_ind, downsampled_end_ind);
    downsampled.resize(downsampled_shape);

    // Iterate over downsampled indices in the range of this array
    auto downsampler = [&](const ind_t& downsampled_ind) {
      ind_t dest_ind = downsampled_ind - downsampled_start_ind;
      ind_t src_ind = downsampled_ind * downsampling - start_ind;

      assert(downsampled.has_index(dest_ind));
      assert(to_downsample.has_index(src_ind));

      downsampled[dest_ind] = to_downsample[src_ind];
    };
    IteratorUtils::index_loop_over_elements(downsampled_start_ind,
                                            downsampled_end_ind, downsampler);
  }
} // namespace Downsampling
