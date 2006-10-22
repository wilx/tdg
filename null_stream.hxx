#ifndef STREAMS_HXX
#define STREAMS_HXX

#include <iostream>


class null_ostream
  : public std::ostream
{
public:
  template <typename T>
  virtual std::null_ostream & 
  operator << (T const &)
  {
    return *this;
  }
};


class prefixed_ostream
//: public std::ostream
{
public:
  prefixed_ostream (std::string const & p, std::ostream & os)
    : pref (p), ostr (os)
  { }

  template <typename T>
  virtual std::ostream &
  operator << (T const & x)
  {
    os << pref << x;
    return os;
  }
  
protected:
  std::string const pref;
  std::ostream & ostr;
};


#endif STREAMS_HXX
