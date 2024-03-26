#include <stdexcept>

#include "Eisvogel/IteratorUtils.hh"
#include "NDVecArraySerialization.hh"

namespace stor {

  // Metadata structures for serialization
  // Metadata for array
  template <typename T, std::size_t dims, std::size_t vec_dims>
  struct NDVecArrayStreamerMetadata {
    
    using shape_t = typename NDVecArrayStreamer<T, dims, vec_dims>::shape_t;
    using stride_t = typename NDVecArrayStreamer<T, dims, vec_dims>::stride_t;   
    
    NDVecArrayStreamerMetadata(const shape_t chunk_size, const shape_t array_shape, const stride_t strides, const std::size_t offset) :
      chunk_size(chunk_size), array_shape(array_shape), strides(strides), offset(offset) { };
    
    shape_t chunk_size;
    shape_t array_shape;
    stride_t strides;
    std::size_t offset;
  };
  
  // Metadata for each serialization chunk
  struct NDVecArrayStreamerChunkMetadata {
    
    NDVecArrayStreamerChunkMetadata(const StreamerMode ser_mode, const std::size_t chunk_size) :
      ser_mode(ser_mode), chunk_size(chunk_size) { };
    
    StreamerMode ser_mode;
    std::size_t chunk_size;
  };    
  
  // Prescriptions for serializing metadata
  template <typename T, std::size_t dims, std::size_t vec_dims>
  struct Traits<NDVecArrayStreamerMetadata<T, dims, vec_dims>> {    
    using type = NDVecArrayStreamerMetadata<T, dims, vec_dims>;
    using shape_t = typename type::shape_t;
    using stride_t = typename type::stride_t;
    
    static void serialize(std::iostream& stream, const type& val) {
      Traits<shape_t>::serialize(stream, val.chunk_size);
      Traits<shape_t>::serialize(stream, val.array_shape);
      Traits<stride_t>::serialize(stream, val.strides);
      Traits<size_t>::serialize(stream, val.offset);
    }

    static type deserialize(std::iostream& stream) {
      shape_t chunk_size(Traits<shape_t>::deserialize(stream));
      shape_t array_shape(Traits<shape_t>::deserialize(stream));
      stride_t strides(Traits<stride_t>::deserialize(stream));
      std::size_t offset = Traits<size_t>::deserialize(stream);

      return type(chunk_size, array_shape, strides, offset);
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

  template <typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<T, dims, vec_dims>::write_buffer(std::size_t num_elems, std::fstream& stream) {
    std::span<ser_type> to_write(m_ser_buffer -> begin(), num_elems);
    stream.write((char*)(&to_write[0]), to_write.size_bytes());
  }

  template <typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<T, dims, vec_dims>::read_into_buffer(std::size_t num_elems, std::fstream& stream) {
    std::span<ser_type> to_read(m_ser_buffer -> begin(), num_elems);
    stream.read((char*)(&to_read[0]), to_read.size_bytes());
  }
  
  template <typename T, std::size_t dims, std::size_t vec_dims>
  NDVecArrayStreamer<T, dims, vec_dims>::NDVecArrayStreamer(std::size_t initial_buffer_size) {
    m_ser_buffer = std::make_shared<buffer_t>(initial_buffer_size);
  }

  template <typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<T, dims, vec_dims>::serialize(std::fstream& stream, const type& val, const shape_t& chunk_size, const StreamerMode& mode) {

    // build array-wide metadata
    // don't directly use the strides the array comes with: these may refer to a view!
    shape_t array_shape = val.m_shape;
    stride_t array_strides = type::ComputeStrides(array_shape);
    std::size_t offset = 0;
    
    NDVecArrayStreamerMetadata<T, dims, vec_dims> meta(chunk_size, array_shape, array_strides, offset);

    // serialize array-wide metadata
    Traits<NDVecArrayStreamerMetadata<T, dims, vec_dims>>::serialize(stream, meta);
    
    switch(mode) {
      
    case StreamerMode::dense:
      serialize_all_chunks_dense(stream, val, chunk_size);
      break;

    case StreamerMode::zero_suppressed:
      serialize_all_chunks_zero_suppressed(stream, val, chunk_size);
      break;

    case StreamerMode::automatic:
      // TODO: implement this so it decides on a per-chunk basis whether to go for dense or zero-suppressed storage
      throw std::runtime_error("Error: streamer mode 'automatic' not implemented yet.");
      break;      
    }    
  }

  template <typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<T, dims, vec_dims>::deserialize(std::fstream& stream, type& val) {

    // deserialize array-wide metadata
    NDVecArrayStreamerMetadata<T, dims, vec_dims> meta = Traits<NDVecArrayStreamerMetadata<T, dims, vec_dims>>::deserialize(stream);

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
      std::size_t elems_read = dense::from_buffer(val_view, std::span<ser_type>(*m_ser_buffer));
    };

    loop_over_array_chunks(val, meta.chunk_size, chunk_deserializer);    
  }

  template <typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<T, dims, vec_dims>::append_slice(std::fstream& stream, const type& chunk, const StreamerMode& mode) {

  }  

  // --------
  
  template <typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<T, dims, vec_dims>::serialize_all_chunks_dense(std::fstream& stream, const type& val, const shape_t& chunk_size) {
        
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
  
  template <typename T, std::size_t dims, std::size_t vec_dims>
  NDVecArrayStreamer<T, dims, vec_dims>::type NDVecArrayStreamer<T, dims, vec_dims>::deserialize_chunk_dense(std::fstream& stream) {
    
    shape_t shape = Traits<shape_t>::deserialize(stream);
    stride_t strides = Traits<stride_t>::deserialize(stream);
    std::size_t offset = Traits<std::size_t>::deserialize(stream);
    data_t data = Traits<data_t>::deserialize(stream);
    return type(shape, strides, offset, std::move(data));
  }

  // ----------------------------------------
  
  template <typename T, std::size_t dims, std::size_t vec_dims>
  void NDVecArrayStreamer<T, dims, vec_dims>::serialize_all_chunks_zero_suppressed(std::fstream& stream, const type& val, const shape_t& chunk_size) {

    Traits<shape_t>::serialize(stream, val.m_shape);
    Traits<stride_t>::serialize(stream, val.m_strides);
    Traits<std::size_t>::serialize(stream, val.m_offset);
    Traits<data_t>::serialize(stream, *val.m_data);
    
    // loop_over_elements();
    
  }

  template <typename T, std::size_t dims, std::size_t vec_dims>
  NDVecArrayStreamer<T, dims, vec_dims>::type NDVecArrayStreamer<T, dims, vec_dims>::deserialize_chunk_zero_suppressed(std::fstream& stream) {
    
    shape_t shape = Traits<shape_t>::deserialize(stream);
    stride_t strides = Traits<stride_t>::deserialize(stream);
    std::size_t offset = Traits<std::size_t>::deserialize(stream);
    data_t data = Traits<data_t>::deserialize(stream);
    return type(shape, strides, offset, std::move(data));
  }
}

