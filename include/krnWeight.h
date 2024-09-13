#ifndef _KRNWEIGHT_H
#define _KRNWEIGHT_H

#include<cmath>

using std::tgamma;
using std::pow;

/* Weight function for L2(T) with Koornwinder basis, where
   T is the standard right triangle */

template<typename T>
inline double 
krnWeight  ( T x, T y, T a, T b, Tc  )
{
  T w = tgamma(a + b + c + 1.5) / ( tgamma(a + 0.5) * tgamma(b + 0.5) * tgamma(c + 0.5) );
  return pow(x, a - 0.5) * pow(y, b - 0.5) * pow(1 - x - y, c - 0.5) * w / 2.0;
}

#endif
