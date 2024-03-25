#pragma once

#include <filesystem>
#include <fstream>

#include "NDVecArrayStreamer.hh"

class ChunkIndex {

public:

  
  
};

template <std::size_t dims, std::size_t vec_dims>
class DistributedNDVecArray {

public:

  DistributedNDVecArray(std::filesystem::path dirpath, std::size_t cache_size);

private:

  const std::filesystem::path m_dirpath;
  const std::size_t m_cache_size;
  
};

#include "DistributedNDVecArray.hxx"
