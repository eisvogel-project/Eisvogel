#pragma once

// allocator that does not perform any initialization
// see https://stackoverflow.com/questions/21028299/is-this-behavior-of-vectorresizesize-type-n-under-c11-and-boost-container
template <typename T, typename A = std::allocator<T>>
class no_init_alloc : public A {
  
  typedef std::allocator_traits<A> a_t;
  
public:
  
  template <typename U> struct rebind {
    using other = no_init_alloc<U, typename a_t::template rebind_alloc<U>>;
  };

  using A::A;

  template <typename U>
  void construct(U* ptr)
    noexcept(std::is_nothrow_default_constructible<U>::value) {
    ::new(static_cast<void*>(ptr)) U;
  }
  
  template <typename U, typename...Args>
  void construct(U* ptr, Args&&... args) {
    a_t::construct(static_cast<A&>(*this),
                   ptr, std::forward<Args>(args)...);
  }
};
