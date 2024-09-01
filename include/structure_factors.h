#ifndef _COMMON_H
#define _COMMON_H

#include<cmath>
#include"exceptions.h"

using std::tgamma;
using std::pow;
using std::abs;
using std::sqrt;

/* 
   Generate the normalization factors for Koornwinder 
   polynomials with parameter (a,b,c). Note, the 
   basis parameters appear as (a-1/2,b-1/2,c-1/2) in the 
   definintion of the polynomials, so the legendre analog
   weight of (0,0,0) corresponds to a=b=c=1/2. 

  Inputs:
    N - polynomial normalizations generated up to order N-1
    a,b,c - Koornwinder parameters
    H - Floating type T pointer to N^2 storage (heap allocated)

  Outputs:
    H - Populated as NxN upper triangular matrix 
      - each column has normalizations for $P_i^j$,
      - where the column is j, and i <= j
      - That is, each column j corresponds to the space of 
      - homogeneous polynomials of degree j
      - NOTE: stored in column major order
*/
   
template<typename T>
inline void structure_factors_tri(const unsigned int N,
                              const T a, const T b, const T c, 
                              T* H)
{
  if (!H || N == 0)
  {
    exitErr("scale array must not be empty, and N cannot be 0!");
  }

  T kap = abs(a+b+c);
  T wabc = tgamma(kap+1.5) / 
           ( tgamma(a+0.5) * tgamma(b+0.5) * tgamma(c+0.5) );

  for (unsigned int nn = 0; nn < N; ++nn)
  {
    #pragma omp simd
    for (unsigned int kk = 0; kk <= nn; ++kk)
    {
      H[kk + N*nn] = sqrt( 
                             wabc/( (2*nn+kap+0.5) * (2*kk+b+c) ) *
                             tgamma(nn+kk+b+c+1) * tgamma(nn-kk+a+0.5) * 
                             tgamma(kk+b+0.5) * tgamma(kk+c+0.5) / 
                              ( tgamma(nn-kk+1) * tgamma(kk+1) * 
                                tgamma(nn+kk+kap+0.5) * tgamma(kk+b+c)
                              )
                           );
    }
  }
}


/*
template<typename T>
inline void structure_factors(const unsigned int N,
                              const T a, const T b, 
                              T* H)
{
  #pragma omp simd
  for (unsigned int i = 0; i < N; ++i)
  {
    H[i] = pow(2, a+b+1) * tgamma(i+a+1) * tgamma(i+b+1) /
           ( (2*i+a+b+1) * tgamma(i+a+b+1) * tgamma(i+1) );
  }
  H[0] = pow(2, a+b+1) * tgamma(a+1) * tgamma(b+1) /
         tgamma(a+b+2);
}
*/

#endif
