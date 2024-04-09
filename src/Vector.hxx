// Serialization
namespace stor {

  template <typename T, std::size_t vec_dims>
  struct Traits<Vector<T, vec_dims>> {
    using type = Vector<T, vec_dims>;
    using data_t = typename type::data_t;

    static void serialize(std::iostream& stream, const type& val) {
      Traits<data_t>::serialize(stream, val.m_data);
    }

    static type deserialize(std::iostream& stream) {
      data_t data = Traits<data_t>::deserialize(stream);
      return Vector<T, vec_dims>(std::move(data));
    }
  };    
}

// The wrapped vectors (de)serialize just like their base classes
namespace stor {

  template <typename T>
  struct Traits<RZTVector<T>> {
    using type = RZTVector<T>;
    
    static void serialize(std::iostream& stream, const type& val) {
      Traits<Vector3D<T>>::serialize(stream, val);
    }
    static type deserialize(std::iostream& stream) {
      return Traits<Vector3D<T>>::deserialize(stream);
    }
  };  
}

namespace VectorUtils {

  template <typename T, std::size_t vec_dims>
  T inner(const Vector<T, vec_dims>& vec_a, const Vector<T, vec_dims>& vec_b) {
    T result = 0;
    result = std::transform_reduce(std::execution::unseq, vec_a.begin(), vec_a.end(), vec_b.begin(), 0.0, std::plus<>(), std::multiplies<>());
    return result;
  }
  
  // `min` and `max` for vectors holding integers (can be signed or unsigned)
  template <std::integral T1, std::integral T2, std::size_t vec_dims>
  Vector<T1, vec_dims> min(const Vector<T1, vec_dims>& vec_a, const Vector<T2, vec_dims>& vec_b) {
    
    Vector<T1, vec_dims> result;
    
    auto take_min = [](const T1& a, const T2& b) -> T1 {
      return std::cmp_less(a, b) ? a : static_cast<T1>(b);
    };
    
    std::transform(std::execution::unseq, vec_a.begin(), vec_a.end(), vec_b.begin(), result.begin(), take_min);
    return result;
  }
  
  template <std::integral T1, std::integral T2, std::size_t vec_dims>
  Vector<T1, vec_dims> max(const Vector<T1, vec_dims>& vec_a, const Vector<T2, vec_dims>& vec_b) {
    
    Vector<T1, vec_dims> result;
    
    auto take_max = [](const T1& a, const T2& b) -> T1 {
      return std::cmp_greater(a, b) ? a : static_cast<T1>(b);
    };
    
    std::transform(std::execution::unseq, vec_a.begin(), vec_a.end(), vec_b.begin(), result.begin(), take_max);
    return result;
  }   
}

