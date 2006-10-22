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
