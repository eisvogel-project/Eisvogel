#include <string>
#include "DistributedNDArray.hh"

int main(int argc, char* argv[]) {

  if(argc < 2) {
    throw std::runtime_error("Error: need to pass path to output directory!");
  }

  std::string wf_path = argv[1];

  std::shared_ptr<stor::DefaultSerializer> ser = std::make_shared<stor::DefaultSerializer>();
  DistributedScalarNDArray<scalar_t, 3> darr(wf_path, 10, *ser);

  darr.printChunks();
}
