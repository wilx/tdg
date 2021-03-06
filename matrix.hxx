/*
Copyright (c) 2003-2007, Václav Haisman

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
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
