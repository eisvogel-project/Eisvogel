#pragma once

#include <fstream>
#include <filesystem>

template <typename T>
class CSVReader {

public:  
  CSVReader(std::filesystem::path path);
  void read_column(std::size_t col_ind, std::vector<T>& dest);

private:
  std::ifstream m_stream;  
};

#include "CSVReader.hxx"
