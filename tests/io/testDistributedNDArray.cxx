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
  
  DistributedNDArray<float, 2> darr("./distarr/", 10);
  darr.RegisterChunk(chunk1, start_ind1);
  darr.RegisterChunk(chunk2, start_ind2);
  darr.Flush();
  
}
