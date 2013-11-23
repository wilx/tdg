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
#ifndef MPI_BUFFER_H
#define MPI_BUFFER_H

#include <vector>
#include <cassert>
#include <stdexcept>
#include <cstdlib>
#include <mpi.h>


template <typename T>
void * type_to_mpi_type (T const &)
{ 
  abort ();
}

#define TYPE_TO_MPI_TYPE(t, mpit) \
  inline MPI_Datatype type_to_mpi_type (t const &) \
  { return mpit; }

TYPE_TO_MPI_TYPE (int, MPI_INT);
TYPE_TO_MPI_TYPE (unsigned, MPI_UNSIGNED);
TYPE_TO_MPI_TYPE (char, MPI_CHAR);
TYPE_TO_MPI_TYPE (signed char, MPI_CHAR);
TYPE_TO_MPI_TYPE (unsigned char, MPI_UNSIGNED_CHAR);
TYPE_TO_MPI_TYPE (short, MPI_SHORT);
TYPE_TO_MPI_TYPE (unsigned short, MPI_UNSIGNED_SHORT);
TYPE_TO_MPI_TYPE (long, MPI_LONG);
TYPE_TO_MPI_TYPE (unsigned long, MPI_UNSIGNED_LONG);
TYPE_TO_MPI_TYPE (float, MPI_FLOAT);
TYPE_TO_MPI_TYPE (double, MPI_DOUBLE);

#undef TYPE_TO_MPI_TYPE


class mpi_obuffer
  : protected std::vector<char>
{
public:
  mpi_obuffer ()
    : pos (0)
  {
    reserve (128);
  }

  mpi_obuffer (mpi_obuffer const & buf)
    : std::vector<char> (buf), pos (buf.pos)
  { }

  template <typename T>
  mpi_obuffer & 
  operator << (T x) 
  {
    resize (size () + sizeof (T));
    int ret = MPI_Pack (&x, 1, type_to_mpi_type (T ()), &*begin (), size (), 
                        &pos, MPI_COMM_WORLD);
    if (ret != MPI_SUCCESS)
      throw std::runtime_error ("MPI_Pack() error");

    return *this;
  }

  char const * 
  data () const
  {
    if (size ())
      return &*begin ();
    else
      return 0;
  }

  char * 
  data ()
  {
    if (size ())
      return &*begin ();
    else
      return 0;
  }

  mpi_obuffer & 
  clear ()
  {
    std::vector<char>::clear ();
    pos = 0;
    return *this;
  }

  ssize_t 
  get_pos () const
  {
    return pos;
  }

  using std::vector<char>::size;

protected:
  ssize_t pos;
};


class mpi_ibuffer 
{
public:
  mpi_ibuffer (void * b, ssize_t s, ssize_t p = 0)
    : buf (b), sz (s), pos (p)
  { }

  mpi_ibuffer (mpi_ibuffer const & b)
    : buf (b.buf), sz (b.sz), pos (b.pos)
  { }

  mpi_ibuffer ()
    : buf (0), sz (0), pos (0)
  { }

  template <typename T>
  mpi_ibuffer &
  operator >> (T & dest)
  {
    int ret = MPI_Unpack (buf, sz, &pos, &dest, 1, type_to_mpi_type (T ()),
                          MPI_COMM_WORLD);
    if (ret != MPI_SUCCESS)
      throw std::runtime_error ("MPI_Unpack error");
    return *this;
  }

  mpi_ibuffer & 
  reset ()
  {
    pos = 0;
    return *this;
  }

  ssize_t 
  size () const
  {
    return sz;
  }

  ssize_t 
  get_pos () const
  {
    return pos;
  }

  void *
  get_buffer () const
  {
    return buf;
  }

  mpi_ibuffer &
  attach (void * b, ssize_t s, ssize_t p = 0)
  {
    buf = b;
    sz = s;
    pos = p;
    return *this;
  }

protected:
  void * buf;
  ssize_t sz;
  ssize_t pos;
};


#endif // MPI_BUFFER_H
