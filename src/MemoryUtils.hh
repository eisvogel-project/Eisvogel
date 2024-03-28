#pragma once

// allocator that does not perform any initialization
// see https://stackoverflow.com/questions/21028299/is-this-behavior-of-vectorresizesize-type-n-under-c11-and-boost-container
template <class T>
class no_init_alloc : public std::allocator<T> {
public:
  using std::allocator<T>::allocator;
  template <class U, class... Args> void construct(U*, Args&&...) {}
};
