namespace stor {

  template <typename T>
  void DefaultSerializer::serialize_impl(std::fstream& stream, const T& value)
  {
    Traits<T>::serialize(stream, value);
  }

  template <typename T>
  T DefaultSerializer::deserialize_impl(std::fstream& stream)
  {
    return Traits<T>::deserialize(stream);
  }

} // namespace stor
