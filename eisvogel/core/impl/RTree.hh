#pragma once

#include <vector>

template <class IndexT, class PayloadT, std::size_t dims, std::size_t max_elems = 5>
class RTree {

public:

  // Constructs empty tree and reserves `init_slot_size` slots to hold elements
  RTree(std::size_t init_slot_size);

  void AddElement(const PayloadT& elem, const IndexT& start_ind, const IndexT& end_ind);

  // Rebuild the tree and rebalance the nodes, if needed
  void Rebuild();
  
  const PayloadT& Search(const IndexT& ind) const;
  PayloadT& Search(const IndexT& ind);
  std::vector<std::reference_wrapper<const PayloadT&>> Search(const IndexT& start_ind, const IndexT& end_ind);

private:

  // A general bounding box with start and end coordinates
  struct BoundingBox {
    IndexT start_ind;
    IndexT end_ind;
  };

  // Internal tree node
  template <class BoundedT>
  struct Node : BoundedT {

    // List of pointers to items that are children of this node
    std::vector<BoundedT*> items;
  };

  // Leaf node of the tree
  template <class BoundedT>
  struct Leaf : BoundedT {

    // Tree leaf node with pointer to the data
    PayloadT* data;
  };

  using RTreeNode = Node<BoundingBox>;
  using RTreeLeaf = Leaf<BoundingBox>;

private:

  // Methods to request new slots for nodes, leaves, and the actual data
  RTreeNode* get_new_node();
  RTreeLeaf* get_new_leaf();
  PayloadT* get_new_payload();
  
private:

  // Contiguous storage for all internal tree nodes and tree leaves
  std::vector<RTreeNode> m_nodes;
  std::vector<RTreeLeaf> m_leaves;
  
  // Contiguous storage for all data members
  std::vector<PayloadT> m_data;
};

#include "RTree.hxx"
