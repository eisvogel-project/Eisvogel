#include <stdexcept>
#include <cassert>
#include <ios>

#include "Eisvogel/IteratorUtils.hh"
#include "NDVecArraySerialization.hh"

namespace stor {

  // Metadata structures for serialization
  // Metadata for array
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  struct NDVecArrayStreamerMetadata {
    
    using shape_t = typename NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::shape_t;
    using stride_t = typename NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::stride_t;   
    
    NDVecArrayStreamerMetadata(const shape_t ser_chunk_shape, const shape_t array_shape, const AccessMode access_mode) :
      ser_chunk_shape(ser_chunk_shape), array_shape(array_shape), access_mode(access_mode) { };
    
    shape_t ser_chunk_shape;
    shape_t array_shape;
    AccessMode access_mode;
  };
  
  // Metadata for each serialization chunk
  struct NDVecArrayStreamerChunkMetadata {
    
    NDVecArrayStreamerChunkMetadata(const StreamerMode ser_mode, const std::size_t chunk_size) :
      ser_mode(ser_mode), chunk_size(chunk_size) { };
    
    StreamerMode ser_mode;
    std::size_t chunk_size;
  };    
  
  // Prescriptions for serializing metadata
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  struct Traits<NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>> {    
    using type = NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>;
    using shape_t = typename type::shape_t;
    using stride_t = typename type::stride_t;
    
    static void serialize(std::iostream& stream, const type& val) {
      Traits<shape_t>::serialize(stream, val.ser_chunk_shape);
      Traits<shape_t>::serialize(stream, val.array_shape);
      Traits<std::size_t>::serialize(stream, static_cast<std::size_t>(val.access_mode));
    }

    static type deserialize(std::iostream& stream) {
      shape_t ser_chunk_shape(Traits<shape_t>::deserialize(stream));
      shape_t array_shape(Traits<shape_t>::deserialize(stream));
      AccessMode access_mode{Traits<std::size_t>::deserialize(stream)};

      return type(ser_chunk_shape, array_shape, access_mode);
    }
  };

  template <>
  struct Traits<NDVecArrayStreamerChunkMetadata> {
    using type = NDVecArrayStreamerChunkMetadata;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<std::size_t>::serialize(stream, static_cast<std::size_t>(val.ser_mode));
      Traits<std::size_t>::serialize(stream, val.chunk_size);
    }

    static type deserialize(std::iostream& stream) {
      StreamerMode ser_mode{Traits<std::size_t>::deserialize(stream)};
      std::size_t chunk_size = Traits<std::size_t>::deserialize(stream);

      return type(ser_mode, chunk_size);
    }
  };

  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::write_buffer(std::size_t num_elems, std::fstream& stream) {
    std::span<ser_type> to_write(m_ser_buffer -> begin(), num_elems);
    stream.write((char*)(&to_write[0]), to_write.size_bytes());
  }

  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::read_into_buffer(std::size_t num_elems, std::fstream& stream) {
    std::span<ser_type> to_read(m_ser_buffer -> begin(), num_elems);
    stream.read((char*)(&to_read[0]), to_read.size_bytes());
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::NDVecArrayStreamer(std::size_t initial_buffer_size) {
    m_ser_buffer = std::make_shared<buffer_t>(initial_buffer_size);
  }

  // tests if a serialization chunk is well-formed, i.e. its extent along at least one direction is explicitly specified
  template <std::size_t dims>
  bool is_serialization_chunk_well_formed(const Vector<std::size_t, dims>& ser_chunk_shape) {
    std::size_t num_dimensions_specified = 0;
    for(std::size_t ind = 0; ind < dims; ind++) {
      if(ser_chunk_shape[ind] != stor::INFTY) {
	num_dimensions_specified++;
      }
    }
    return num_dimensions_specified >= 1;
  }

  // tests if during the serialization of an array with `array_shape` all serialization chunks will have the same shape (as specified by the requested `ser_chunk_shape`,
  // which will be different from the actual serialization chunk size in case stor::INFTY is used in its specification)
  template <std::size_t dims>
  bool permits_equally_sized_serialization_chunks(const Vector<std::size_t, dims>& array_shape, const Vector<std::size_t, dims>& ser_chunk_shape) {
    for(std::size_t ind = 0; ind < dims; ind++) {
      if((ser_chunk_shape[ind] < array_shape[ind]) && (array_shape[ind] % ser_chunk_shape[ind] != 0)) {
	return false;
      }
    }
    return true;
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::serialize(std::fstream& stream, const type& val, const shape_t& ser_chunk_shape, const StreamerMode& mode) {

    // check if the specified serialization chunk size is well-formed
    if(!is_serialization_chunk_well_formed(ser_chunk_shape)) {
      throw std::runtime_error("Error: requested serialization chunk is not well-formed!");
    }
    
    // check if the passed array can be described by a number of equally-sized serialization chunks ...
    // ... which determines whether or not appending slices on-disk is allowed later
    AccessMode access_mode = permits_equally_sized_serialization_chunks(val.m_shape, ser_chunk_shape) ? AccessMode::modification_allowed : AccessMode::modification_not_allowed;
    
    // build array-wide metadata
    // don't directly use the strides the array comes with: these may refer to a view!
    shape_t array_shape = val.m_shape;
    NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims> meta(ser_chunk_shape, array_shape, access_mode);

    // serialize array-wide metadata
    Traits<NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>>::serialize(stream, meta);

    // serialize array chunks
    serialize_all_chunks(stream, val, ser_chunk_shape, mode);
  }

  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::serialize_all_chunks(std::fstream& stream, const type& val, const shape_t& ser_chunk_shape, const StreamerMode& mode) {
    
    switch(mode) {
      
    case StreamerMode::dense:
      serialize_all_chunks_dense(stream, val, ser_chunk_shape);
      break;

    case StreamerMode::null_suppressed:
      serialize_all_chunks_null_suppressed(stream, val, ser_chunk_shape);
      break;

    case StreamerMode::automatic:
      // TODO: implement this so it decides on a per-chunk basis whether to go for dense or null-suppressed storage
      throw std::runtime_error("Error: streamer mode 'automatic' not implemented yet.");
      break;      
    }    
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::deserialize(std::fstream& stream, type& val) {

    // deserialize array-wide metadata
    NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims> meta = Traits<NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>>::deserialize(stream);
    
    // prepare new array of the correct shape
    val.resize(meta.array_shape);

    // loop over chunks
    auto chunk_deserializer = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>& chunk_end) -> void {

      // get chunk metadata
      NDVecArrayStreamerChunkMetadata chunk_meta = Traits<NDVecArrayStreamerChunkMetadata>::deserialize(stream);

      // make sure buffer is large enough and read data
      m_ser_buffer -> resize(chunk_meta.chunk_size);
      read_into_buffer(chunk_meta.chunk_size, stream);
      
      // deserialize data and fill into array
      type val_view = val.View(chunk_begin, chunk_end);
      
      std::size_t elems_read;
      switch(chunk_meta.ser_mode) {
	
      case StreamerMode::dense:
	  elems_read = dense::from_buffer(val_view, std::span<ser_type>(*m_ser_buffer));
	  break;
	
      case StreamerMode::null_suppressed:
	  elems_read = nullsup::desuppress_null(std::span<ser_type>(*m_ser_buffer), val_view);
	  break;	
      }
      assert(elems_read == chunk_meta.chunk_size);
    };

    IteratorUtils::index_loop_over_array_chunks(val, meta.ser_chunk_shape, chunk_deserializer);    
  }

  // tests if a chunk with `ser_chunk_shape` is a `slice` of the array with shape `array_shape`
  // A chunk is a `slice` if (as the name suggests) fully "slices through" the array along a particular direction
  template <std::size_t dims>
  bool chunk_is_slice(const Vector<std::size_t, dims>& array_shape, const Vector<std::size_t, dims>& ser_chunk_shape) {
    std::size_t num_dimensions_not_sliced = 0;
    for(std::size_t ind = 0; ind < dims; ind++) {
      if(ser_chunk_shape[ind] < array_shape[ind]) {
	num_dimensions_not_sliced++;
      }
    }
    // this chunk is a slice if at most one direction is not fully sliced through
    return num_dimensions_not_sliced <= 1;
  }

  // determines the axis along which the serialization chunk loop proceeds, assuming that the serialization chunk is a slice and that there is a unique
  // serialization axis (which is then defined as the unique axis that the serialization chunk does not fully slice through)
  template <std::size_t dims>
  std::size_t determine_unique_serialization_axis(const Vector<std::size_t, dims>& array_shape, const Vector<std::size_t, dims>& ser_chunk_shape) {
    for(std::size_t ind = 0; ind < dims; ind++) {
      if(ser_chunk_shape[ind] < array_shape[ind]) {
	return ind;
      }
    }
    throw std::logic_error("This function assumes a serialization chunk strictly smaller than the full array!");
  }
  
  // tests if the serialization axis is compatible with the `axis` along which the new slice should be inserted
  template <std::size_t dims>
  bool ser_axis_matches_insertion_axis(const Vector<std::size_t, dims>& array_shape, const Vector<std::size_t, dims>& ser_chunk_shape, std::size_t axis) {

    // if the serialization chunk contains the full array, any insertion axis is allowed ...
    if(array_shape == ser_chunk_shape) {
      return true;
    }
    
    // ... otherwise, there is a unique serialization axis, which must match the insertion axis
    if(determine_unique_serialization_axis(array_shape, ser_chunk_shape) == axis) {
      return true;
    }

    return false;
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::append_slice(std::fstream& stream, const type& slice, std::size_t axis, const StreamerMode& mode) {

    // keep track of current (i.e. beginning) location on stream
    std::streampos init_pos = stream.tellg();
    
    // deserialize array-wide metadata
    NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims> meta = Traits<NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>>::deserialize(stream);

    if(meta.access_mode == AccessMode::modification_not_allowed) {
      throw std::runtime_error("Error: trying to modify an array that cannot be modified!");
    }

    // need to make sure that the new slice has the right dimensions for appending to the array along the prescribed axis
    if(!ArrayT<T, dims, vec_dims>::template ShapeAllowsConcatenation(meta.array_shape, slice.m_shape, axis)) {
      throw std::runtime_error("Error: dimensions not compatible for concatenation!");
    }

    // --------
    // TODO: Finish implementing checks
    // must satisfy that it(part1) + it(part2) === it(part1 + part2) for the given chunk size and dimensions of part1, part2

    // Better conditions:
    // 1) must have no gaps and no automatically-growing chunks along concatenation dimension: i.e. must not have partial chunks along concat dimension
    //       (ser_chunk_shape[concat_dim] <= array_shape[concat_dim]) && (array_shape[concat_dim] % ser_chunk_shape[concat_dim] == 0)
    // 2) concat dimension must not be after the outermost changing dimension
    // --------
    
    // // check if the serialization chunk is a `slice`, i.e. the chunk fully "slices through" the array along a particular dimension
    // if(!chunk_is_slice(meta.array_shape, meta.ser_chunk_shape)) {
    //   throw std::runtime_error("Error: serialization chunk is not a slice, cannot append!");
    // }
    
    // // need to make sure that the serialization axis is the same as the axis for appending the new slice
    // if(!ser_axis_matches_insertion_axis(meta.array_shape, meta.ser_chunk_shape, axis)) {
    //   throw std::runtime_error("Error: serialization axis does not match requested insertion axis!");
    // }

    // --------
        
    // update metadata
    // check if the passed `slice` fits exactly into one or more serialization chunks; no further modifications are allowed if it doesn't
    meta.access_mode = permits_equally_sized_serialization_chunks(slice.m_shape, meta.ser_chunk_shape) ? AccessMode::modification_allowed : AccessMode::modification_not_allowed;

    // update array size after `slice` is appended
    meta.array_shape[axis] += slice.m_shape[axis];

    // update metadata on disk
    stream.seekp(init_pos);
    Traits<NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>>::serialize(stream, meta);

    // append new slice to the end of the file
    stream.seekp(0, std::ios_base::end);
    serialize_all_chunks(stream, slice, meta.ser_chunk_shape, mode);
  }  

  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::mark_as_final(std::fstream& stream) {

    // keep track of current (i.e. beginning) location on stream
    std::streampos init_pos = stream.tellg();

    // deserialize array-wide metadata
    NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims> meta = Traits<NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>>::deserialize(stream);

    // mark as final if not done already
    if(meta.access_mode != AccessMode::modification_not_allowed) {
      meta.access_mode = AccessMode::modification_not_allowed;
      stream.seekp(init_pos);
      Traits<NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>>::serialize(stream, meta);
    }
  }
  
  // --------
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::serialize_all_chunks_dense(std::fstream& stream, const type& val, const shape_t& ser_chunk_shape) {
        
    auto chunk_serializer = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>& chunk_end) -> void {

      // make sure serialization buffer is large enough
      std::size_t ser_buflen = NDVecArray<T, dims, vec_dims>::ComputeVolume(chunk_end - chunk_begin);
      m_ser_buffer -> resize(ser_buflen);
      
      // fill serialization buffer
      std::size_t elems_written = dense::to_buffer(val.View(chunk_begin, chunk_end), std::span<ser_type>(*m_ser_buffer));

      // prepare and serialize chunk metadata
      NDVecArrayStreamerChunkMetadata chunk_meta(StreamerMode::dense, elems_written);
      Traits<NDVecArrayStreamerChunkMetadata>::serialize(stream, chunk_meta);
      
      // serialize chunk data from buffer
      write_buffer(elems_written, stream);
    };
      
    IteratorUtils::index_loop_over_array_chunks(val, ser_chunk_shape, chunk_serializer);
  }

  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::serialize_all_chunks_null_suppressed(std::fstream& stream, const type& val, const shape_t& ser_chunk_shape) {

    auto chunk_serializer = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>& chunk_end) -> void {

      // make sure serialization buffer is large enough
      std::size_t ser_buflen = nullsup::calculate_required_buflen(chunk_end - chunk_begin, vec_dims);
      m_ser_buffer -> resize(ser_buflen);

      // fill serialization buffer
      std::size_t elems_written = nullsup::suppress_null(val.View(chunk_begin, chunk_end), std::span<ser_type>(*m_ser_buffer));
      
      // prepare and serialize chunk metadata
      NDVecArrayStreamerChunkMetadata chunk_meta(StreamerMode::null_suppressed, elems_written);
      Traits<NDVecArrayStreamerChunkMetadata>::serialize(stream, chunk_meta);

      // serialize chunk data from buffer
      write_buffer(elems_written, stream);      
    };
    
    IteratorUtils::index_loop_over_array_chunks(val, ser_chunk_shape, chunk_serializer);    
  }  
}

