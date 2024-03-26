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
    
    NDVecArrayStreamerMetadata(const shape_t chunk_size, const shape_t array_shape, const stride_t strides, const std::size_t offset, const AccessMode access_mode) :
      chunk_size(chunk_size), array_shape(array_shape), strides(strides), offset(offset), access_mode(access_mode) { };
    
    shape_t chunk_size;
    shape_t array_shape;
    stride_t strides;
    std::size_t offset;
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
      Traits<stride_t>::serialize(stream, val.strides);
      Traits<std::size_t>::serialize(stream, val.offset);
      Traits<std::size_t>::serialize(stream, static_cast<std::size_t>(val.access_mode));
    }

    static type deserialize(std::iostream& stream) {
      shape_t chunk_size(Traits<shape_t>::deserialize(stream));
      shape_t array_shape(Traits<shape_t>::deserialize(stream));
      stride_t strides(Traits<stride_t>::deserialize(stream));
      std::size_t offset = Traits<std::size_t>::deserialize(stream);
      AccessMode access_mode{Traits<std::size_t>::deserialize(stream)};

      return type(chunk_size, array_shape, strides, offset, access_mode);
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
    stride_t array_strides = type::ComputeStrides(array_shape);
    std::size_t offset = 0;
    
    NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims> meta(chunk_size, array_shape, array_strides, offset, access_mode);

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

    // deserialize array-wide metadata
    NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims> meta = Traits<NDVecArrayStreamerMetadata<ArrayT, T, dims, vec_dims>>::deserialize(stream);

    // check if the serialization chunk is also a `slice`, i.e. the chunk fully "slices through" the array along a particular dimension

    // check if the passed `chunk` fits exactly into one or more serialization chunks
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

