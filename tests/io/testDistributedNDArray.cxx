#include "Eisvogel/DistributedNDArray.hh"
#include "Eisvogel/NDArray.hh"

#include <fstream>

int main(int argc, char* argv[]) {

  DenseNDArray<float, 2> chunk1({2, 2}, 0.0);
  chunk1(0, 0) = 1.0;
  chunk1(0, 1) = 2.0;
  chunk1(1, 1) = 4.0;
  chunk1.print();

  IndexVector start_ind1({0, 0});
  start_ind1.print();

  DenseNDArray<float, 2> chunk2({2, 2}, 0.0);
  chunk2(0, 0) = -1.0;
  chunk2(0, 1) = -2.0;
  chunk2(1, 1) = -4.0;
  chunk2.print();

  IndexVector start_ind2({2, 0});
  start_ind2.print();
  
  DistributedNDArray<float, 2> arr("./distarr/", 10);
  arr.RegisterChunk(chunk1, start_ind1);
  arr.RegisterChunk(chunk2, start_ind2);

  ChunkMetadata chunk_meta1("bla", start_ind1, {1, 1});

  std::fstream ofs;
  ofs.open("meta_ser.bin", std::ios::out | std::ios::binary);  
  stor::Serializer oser(ofs);

  std::string test_string = "test string";

  oser.serialize(chunk_meta1);
  ofs.close();


  std::fstream ifs; 
  ifs.open("meta_ser.bin", std::ios::in | std::ios::binary);
  stor::Serializer iser(ifs);
  ChunkMetadata read_meta = iser.deserialize<ChunkMetadata>();

  std::cout << read_meta.filename << std::endl;
  read_meta.start_ind.print();
  read_meta.stop_ind.print();
  
}
