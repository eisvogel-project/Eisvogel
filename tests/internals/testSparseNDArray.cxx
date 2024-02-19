#include <iostream>
#include <fstream>
#include "Serialization.hh"
#include "DenseNDArray.hh"
#include "SparseNDArray.hh"

int main(int argc, char* argv[]) {

  DenseNDArray<float, 2> darr({3, 3}, 0.0);

  darr(0, 0) = 1.0;
  darr(1, 0) = 0.4;
  darr(1, 1) = 1.4;
  
  auto to_keep = [](float value) -> bool {
    return abs(value) > 1e-6;
  };
  
  SparseNDArray<float, 2> sparr = SparseNDArray<float, 2>::From(darr, to_keep, 0.0);

  std::string ser_path = "sparse_ser_test.bin";
  
  std::fstream ofs;
  ofs.open(ser_path, std::ios::out | std::ios::binary);  
  stor::DefaultSerializer oser;
  oser.serialize(ofs, sparr);
  ofs.close();

  std::fstream ifs;
  ifs.open(ser_path, std::ios::in | std::ios::binary);
  stor::DefaultSerializer iser; 
  SparseNDArray<float, 2> sparr_res = iser.deserialize<SparseNDArray<float, 2>>(ifs);
  ifs.close();
  
  IndexVector acc_ind1 = {1,1};
  std::cout << sparr_res(acc_ind1) << std::endl;

  DenseNDArray<float, 2> darr_res = DenseNDArray<float, 2>::From(sparr_res);
  std::cout << darr_res(acc_ind1) << std::endl;
  
  std::cout << "done" << std::endl;
}
