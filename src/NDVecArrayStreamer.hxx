#include <stdexcept>
#include <cassert>

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
    
    NDVecArrayStreamerMetadata(const shape_t chunk_size, const shape_t array_shape, const AccessMode access_mode) :
      chunk_size(chunk_size), array_shape(array_shape), access_mode(access_mode) { };
    
    shape_t chunk_size;
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
      Traits<shape_t>::serialize(stream, val.chunk_size);
      Traits<shape_t>::serialize(stream, val.array_shape);
      Traits<std::size_t>::serialize(stream, static_cast<std::size_t>(val.access_mode));
    }

    static type deserialize(std::iostream& stream) {
      shape_t chunk_size(Traits<shape_t>::deserialize(stream));
      shape_t array_shape(Traits<shape_t>::deserialize(stream));
      AccessMode access_mode{Traits<std::size_t>::deserialize(stream)};

      return type(chunk_size, array_shape, access_mode);
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
  bool is_serialization_chunk_well_formed(const Vector<std::size_t, dims>& chunk_size) {
    std::size_t num_dimensions_specified = 0;
    for(std::size_t ind = 0; ind < dims; ind++) {
      if(chunk_size[ind] != stor::INFTY) {
	num_dimensions_specified++;
      }
    }
    return num_dimensions_specified >= 1;
  }

  // tests if during the serialization of an array with `array_shape` all serialization chunks will have the same shape (as specified by the requested `chunk_size`,
  // which will be different from the actual serialization chunk size in case stor::INFTY is used in its specification)
  template <std::size_t dims>
  bool permits_equally_sized_serialization_chunks(const Vector<std::size_t, dims>& array_shape, const Vector<std::size_t, dims>& chunk_size) {
    for(std::size_t ind = 0; ind < dims; ind++) {
      if((chunk_size[ind] < array_shape[ind]) && (array_shape[ind] % chunk_size[ind] != 0)) {
	return false;
      }
    }
    return true;
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::serialize(std::fstream& stream, const type& val, const shape_t& chunk_size, const StreamerMode& mode) {

    // check if the specified serialization chunk size is well-formed
    if(!is_serialization_chunk_well_formed(chunk_size)) {
      throw std::runtime_error("Error: requested serialization chunk is not well-formed!");
    }
    
    // check if the passed array can be described by a number of equally-sized serialization chunks ...
    // ... which determines whether or not appending slices on-disk is allowed later
    AccessMode access_mode = permits_equally_sized_serialization_chunks(val.m_shape, chunk_size) ? AccessMode::modification_allowed : AccessMode::modification_not_allowed;
    
    // build array-wide metadata
    // don't directly use the strides the array comes with: these may refer to a view!
    shape_t array_shape = val.m_shape;
    NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims> meta(chunk_size, array_shape, access_mode);

    // serialize array-wide metadata
    Traits<NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>>::serialize(stream, meta);
    
    switch(mode) {
      
    case StreamerMode::dense:
      serialize_all_chunks_dense(stream, val, chunk_size);
      break;

    case StreamerMode::null_suppressed:
      serialize_all_chunks_null_suppressed(stream, val, chunk_size);
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
      m_ser_buffer -> reserve(chunk_meta.chunk_size);
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

    loop_over_array_chunks(val, meta.chunk_size, chunk_deserializer);    
  }

  // tests if a chunk with `chunk_size` is a `slice` of the array with shape `array_shape`
  // A chunk is a `slice` if (as the name suggests) fully "slices through" the array along a particular direction
  template <std::size_t dims>
  bool chunk_is_slice(const Vector<std::size_t, dims>& array_shape, const Vector<std::size_t, dims>& chunk_size) {
    std::size_t num_dimensions_not_sliced = 0;
    for(std::size_t ind = 0; ind < dims; ind++) {
      if(chunk_size[ind] < array_shape[ind]) {
	num_dimensions_not_sliced++;
      }
    }
    // this chunk is a slice if at most one direction is not fully sliced through
    return num_dimensions_not_sliced <= 1;
  }
  
  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::append_slice(std::fstream& stream, const type& chunk, const StreamerMode& mode) {

    // keep track of current (i.e. beginning) location on stream
    std::streampos init_pos = stream.tellg();
    
    // deserialize array-wide metadata
    NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims> meta = Traits<NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>>::deserialize(stream);

    if(meta.access_mode == AccessMode::modification_not_allowed) {
      throw std::runtime_error("Error: trying to modify an array that cannot be modified!");
    }

    // check if the serialization chunk is also a `slice`, i.e. the chunk fully "slices through" the array along a particular dimension
    if(!chunk_is_slice(meta.array_shape, chunk.m_shape)) {
      throw std::runtime_error("Error: passed array chunk is not a slice, cannot append!");
    }

    // update metadata
    // check if the passed `chunk` fits exactly into one or more serialization chunks; no further modifications are allowed if it doesn't
    meta.access_mode = permits_equally_sized_serialization_chunks(chunk.m_shape, meta.chunk_size) ? AccessMode::modification_allowed : AccessMode::modification_not_allowed;

    // update array size after `chunk` is appended
    meta.array_shape += chunk.m_shape;
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
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::serialize_all_chunks_dense(std::fstream& stream, const type& val, const shape_t& chunk_size) {
        
    auto chunk_serializer = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>& chunk_end) -> void {

      // make sure serialization buffer is large enough
      std::size_t ser_buflen = NDVecArray<T, dims, vec_dims>::ComputeVolume(chunk_end - chunk_begin);
      m_ser_buffer -> reserve(ser_buflen); 
      
      // fill serialization buffer
      std::size_t elems_written = dense::to_buffer(val.View(chunk_begin, chunk_end), std::span<ser_type>(*m_ser_buffer));

      // prepare and serialize chunk metadata
      NDVecArrayStreamerChunkMetadata chunk_meta(StreamerMode::dense, elems_written);
      Traits<NDVecArrayStreamerChunkMetadata>::serialize(stream, chunk_meta);
      
      // serialize chunk data from buffer
      write_buffer(elems_written, stream);
    };
      
    loop_over_array_chunks(val, chunk_size, chunk_serializer);
  }

  template <template<typename, std::size_t, std::size_t> class ArrayT,
	    typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<ArrayT, T, dims, vec_dims>::serialize_all_chunks_null_suppressed(std::fstream& stream, const type& val, const shape_t& chunk_size) {

    auto chunk_serializer = [&](const Vector<std::size_t, dims>& chunk_begin, const Vector<std::size_t, dims>& chunk_end) -> void {

      // make sure serialization buffer is large enough
      std::size_t ser_buflen = nullsup::calculate_required_buflen(chunk_end - chunk_begin, vec_dims);
      m_ser_buffer -> reserve(ser_buflen);

      // fill serialization buffer
      std::size_t elems_written = nullsup::suppress_null(val.View(chunk_begin, chunk_end), std::span<ser_type>(*m_ser_buffer));

      // prepare and serialize chunk metadata
      NDVecArrayStreamerChunkMetadata chunk_meta(StreamerMode::null_suppressed, elems_written);
      Traits<NDVecArrayStreamerChunkMetadata>::serialize(stream, chunk_meta);

      // serialize chunk data from buffer
      write_buffer(elems_written, stream);      
    };
    
    loop_over_array_chunks(val, chunk_size, chunk_serializer);    
  }  
}

