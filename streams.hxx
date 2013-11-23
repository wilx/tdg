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
#ifndef STREAMS_HXX
#define STREAMS_HXX

#include <iostream>


struct null_ostream 
{
  template <class T> 
  null_ostream & 
  operator << (T const &) 
  { 
    return *this; 
  }

  template <class T> 
  null_ostream & 
  operator << (T const *) 
  { 
    return *this; 
  }

  null_ostream & 
  operator << (std::ostream & (*) (std::ostream &)) 
  { 
    return *this; 
  }
};


class prefixed_ostream
{
public:
  prefixed_ostream (std::string const & p, std::ostream & os)
    : pref (p), ostr (os)
  { }

  prefixed_ostream (std::ostream & os)
    : pref (""), ostr (os)
  { }

  template <typename T>
  std::ostream &
  operator << (T const & x)
  {
    ostr << pref << x;
    return ostr;
  }

  template <typename T>
  std::ostream &
  operator << (T const * x)
  {
    ostr << pref << x;
    return ostr;
  }

  std::ostream & 
  operator << (std::ostream & (*f) (std::ostream &)) 
  { 
    ostr << f << pref;
    return ostr;
  }
  
  prefixed_ostream &
  set_prefix (std::string const & p)
  {
    pref = p;
    return *this;
  }

protected:
  std::string pref;
  std::ostream & ostr;
};


#endif // STREAMS_HXX
