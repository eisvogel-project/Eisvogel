#pragma once

#include <vector>

template <class IndexT, class PayloadT, std::size_t dims, std::size_t max_elems = 5>
class RTree {

public:
  
  RTree();

  void AddElement(const PayloadT& elem, const IndexT& start_ind, const IndexT& end_ind);

  // Rebuild the tree and rebalance the nodes, if needed
  void Rebuild();
  
  const PayloadT& Search(const IndexT& ind) const;
  PayloadT& Search(const IndexT& ind);
  std::vector<std::reference_wrapper<const PayloadT&>> Search(const IndexT& start_ind, const IndexT& end_ind);

private:

  struct Node {

  };
  
private:

  // Contiguous storage for all data members
  std::vector<PayloadT> m_data;
};

#include "RTree.hxx"
