template <class IndexT, class PayloadT, std::size_t dims, std::size_t max_elems>
RTree<IndexT, PayloadT, dims, max_elems>::RTree(std::size_t init_slot_size) : m_nodes(init_slot_size), m_leaves(init_slot_size), m_data(init_slot_size) { }

template <class IndexT, class PayloadT, std::size_t dims, std::size_t max_elems>
void RTree<IndexT, PayloadT, dims, max_elems>::AddElement(const PayloadT& elem, const IndexT& start_ind, const IndexT& end_ind) {

}

template <class IndexT, class PayloadT, std::size_t dims, std::size_t max_elems>
void RTree<IndexT, PayloadT, dims, max_elems>::Rebuild() {

}

// template <class IndexT, class PayloadT, std::size_t dims, std::size_t max_elems>
// const PayloadT& RTree<IndexT, PayloadT, dims, max_elems>::Search(const IndexT& ind) {

// }

