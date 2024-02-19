#include "DistributedNDArray.hh"
#include "DenseNDArray.hh"

#include <fstream>

int main(int argc, char* argv[]) {

  DenseNDArray<float, 2> chunk1({2, 2}, 0.0);
  chunk1(0, 0) = 1.0;
  chunk1(0, 1) = 2.0;
  chunk1(1, 1) = 4.0;
  chunk1.print();

  DenseNDArray<float, 2> chunk1_1({2, 2}, 0.0);
  chunk1_1.print();

  DenseNDArray<bool, 2> comp = chunk1_1 < chunk1;
  comp.print();
  
  IndexVector start_ind1({0, 0});
  start_ind1.print();

  DenseNDArray<float, 2> chunk2({2, 2}, 0.0);
  chunk2(0, 0) = -1.0;
  chunk2(0, 1) = -2.0;
  chunk2(1, 1) = -4.0;
  chunk2.print();

  IndexVector start_ind2({2, 0});
  start_ind2.print();

  std::shared_ptr<stor::DefaultSerializer> ser = std::make_shared<stor::DefaultSerializer>();

  auto to_keep = [](float value) -> bool {
    return true;
  };
  
  DistributedScalarNDArray<float, 2> darr_save("./distarr/", 10, *ser);
  
  darr_save.RegisterChunk(chunk1, start_ind1);
  darr_save.RegisterChunk(chunk2, start_ind2);
  darr_save.MakeIndexPersistent();
  
  DistributedScalarNDArray<float, 2> darr_load("./distarr/", 10, *ser);

  // IndexVector requested_chunk_size = {10, 10};
  // darr_load.RebuildChunks(requested_chunk_size);
  // darr_load.RebuildChunks(requested_chunk_size);

  std::cout << "HHH BEFORE MERGING HHH" << std::endl;
  darr_load.printChunks();

  darr_load.MergeChunks(0, 10);

  std::cout << "HHH AFTER MERGING HHH" << std::endl;
  darr_load.printChunks();
  
  // IndexVector acc_ind1 = {1,1};
  // std::cout << darr_load(acc_ind1) << std::endl;
  // std::cout << darr_load(acc_ind1) << std::endl;

  // IndexVector acc_ind2 = {2,1};
  // std::cout << darr_load(acc_ind2) << std::endl;
  // std::cout << darr_load(acc_ind2) << std::endl;

  // std::cout << "shape:" << std::endl;  
  // for(auto cur : darr_load.shape()) {
  //   std::cout << cur << std::endl;
  // }
}
