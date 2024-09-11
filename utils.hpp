#include <sstream>
#include <string> 

inline int intParse(std::string &str, double defvalue = -1){
    int val;
    std::istringstream iss(str); 
    if(!(iss >> val)) return defvalue;
    return val;
}
inline int doubleParse(std::string &str, double defvalue = -1){
    double val;
    std::istringstream iss(str);
    if(!(iss >> val)) return defvalue;
    return val;
}

template<typename var>
inline std::string stringParse(const var val){
    std::stringstream ss;
    ss << val;
    return ss.str();
}

static unsigned long xor128() 
{
  static unsigned long x=123456789, y=362436069, z=521288629, w=88675123;
  unsigned long t;
  t=(x^(x<<11)); x=y; y=z; z=w;
  return (w=(w^(w>>19))^(t^(t>>8)));
}
