#pragma once

#include <stdlib.h>

template<typename T, size_t S>
class SignedArray {
private:
  std::array<T, 2 * S + 1> elements;

public:
  SignedArray() 
    : elements{} 
  {
  }

  void set(const SignedArray& other)
  {
    std::copy(other.elements.begin(), other.elements.end(), elements.begin());
  }

  // This exists so that we can index clusters by their size, and have clusters of negative size
  T& operator[](int index) 
  {
    return elements[S + index];
  }

  size_t size() 
  {
    return S;
  }
};