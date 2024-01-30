#include "Eisvogel/IteratorUtils.hh"

bool isInIndexRange(const IndexVector& inds, const IndexVector& start_inds, const IndexVector& stop_inds) {

  if((inds.size() != start_inds.size()) || (inds.size() != stop_inds.size())) {
    throw;
  }
  
  for(std::size_t dim = 0; dim < inds.size(); dim++) {
    if((inds(dim) < start_inds(dim)) || (inds(dim) >= stop_inds(dim))) {
      return false;
    }
  }

  return true;
}
