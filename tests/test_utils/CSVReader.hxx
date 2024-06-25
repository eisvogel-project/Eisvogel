#include <vector>
#include <string>
#include <sstream>

template <typename T>
class CSVRow {

public:
  void parse_row(std::istream& stream);
  const T& operator[](std::size_t ind) const;

private:
  std::vector<T> m_row_data;
};

template <typename T>
std::istream& operator>>(std::istream& stream, CSVRow<T>& row)
{
  row.parse_row(stream);
  return stream;
}

// ---------------

namespace {

  template <typename T>
  struct StringTraits;

  template <>
  struct StringTraits<float> {
    using type = float;

    static type from_string(const std::string& in) { return std::stof(in); }
  };
} // namespace

template <typename T>
void CSVRow<T>::parse_row(std::istream& stream)
{

  m_row_data.clear();

  std::string cur_line;
  std::getline(stream, cur_line);

  if (cur_line.empty()) {
    return;
  }

  std::stringstream linestream(cur_line);
  std::string cur_field;
  while (linestream.good()) {
    std::getline(linestream, cur_field, ',');
    m_row_data.push_back(StringTraits<T>::from_string(cur_field));
  }
}

template <typename T>
const T& CSVRow<T>::operator[](std::size_t ind) const
{
  return m_row_data[ind];
}

template <typename T>
CSVReader<T>::CSVReader(std::filesystem::path path)
    : m_stream(path)
{
}

template <typename T>
void CSVReader<T>::read_column(std::size_t col_ind, std::vector<T>& dest)
{

  dest.clear();

  // rewind stream
  m_stream.clear();
  m_stream.seekg(0);

  CSVRow<float> row;
  while (m_stream >> row) { dest.push_back(row[col_ind]); }
}
