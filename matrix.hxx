#ifndef _MATRIX_HXX_
#define _MATRIX_HXX_

#include <vector>


template <typename T>
class UTMatrix
{
public:
  typedef std::vector<T> container_type;

  UTMatrix (size_t n)
    : N (n), data (elems_from_n (), T ())
  { }


  typename container_type::reference 
  operator () (size_t x, size_t y)
  {
    return data[index_for_xy (x, y)];
  }

  typename container_type::const_reference 
  operator () (size_t x, size_t y) const
  {
    return data[index_for_xy (x, y)];
  }

  size_t 
  getN () const
  {
    return N;
  }
  
protected:
  size_t 
  elems_from_n () const
  {
    return (N * (N + 1ul)) / 2;
  }

  size_t 
  index_for_xy (size_t x, size_t y) const
  {
    if (x < 1 || x > N || y < 1 || y > N)
      abort ();
    
    if (y > x)
      std::swap (x, y);
    return ((y - 1) * (N + (N - (y-1-1)))) / 2
      + (x - (y - 1)) - 1;
  }


  size_t const N;
  container_type data;
};


#endif // _MATRIX_HXX_
