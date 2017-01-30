#ifndef CUSH_SIGN_H_
#define CUSH_SIGN_H_

#include <config.h>

namespace cush
{
template <typename type> 
INLINE COMMON int sign(type value) 
{
  return (type(0) < value) - (value < type(0));
}
}

#endif
