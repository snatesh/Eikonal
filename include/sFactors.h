#ifndef _SFACTORS_H
#define _SFACTORS_H

#include<cmath>
#include<iostream>
#include<omp.h>

using std::tgamma;
using std::pow;
using std::abs;
using std::sqrt;


/*
   The normalization factors for Jacobi
   polynomials on the tetrahedron with parameter (a,b,c,d). 
   Note, the basis parameters appear as (a-1/2,b-1/2,c-1/2,d-1/2) 
   in the definintion of the polynomials, so the legendre analog
   weight of (0,0,0,0) corresponds to a=b=c=d=1/2. 

*/  

template<typename T>
inline T pochhammer(const T x, const unsigned int n)
{
  T prod = 1;
  # pragma omp simd
  for (unsigned int k = 0; k < n; ++k)
  {
    prod *= (x-k);
  }
  return prod;
}

template<typename T>
inline T sFactors ( const unsigned int n, 
                    const unsigned int k,
                    const unsigned int j,
                    const T a, const T b,
                    const T c, const T d  )
{
  T fac = 
    ((1. + 2.*(-0.5 + c + d + j))*(1. + 2.*(b + c + d + j + k))*
    (1. + 2.*(0.5 + a + b + c + d + k + n))*pochhammer(0.5 + a, -1.*k + n)*
    pochhammer(0.5 + b, -1.*j + k)*pochhammer(0.5 + c, j)*
    pochhammer(0.5 + d, j)*pochhammer(1. + c + d, j + k)*
    pochhammer(1.5 + b + c + d, k + n))/((1. + 2.*(-0.5 + c + d + 2.*j))*
    (1. + 2.*(b + c + d + 2.*j + 2.*(-1.*j + k)))*
    (1. + 2.*(0.5 + a + b + c + d + 2.*k + 2.*(-1.*k + n)))*
    tgamma(j+1)*tgamma(-1.*j + k+1)*
    tgamma(-1.*k + n+1)*pochhammer(1. + c + d, j)*pochhammer(1.5 + b + c + d, j + k)*
    pochhammer(2. + a + b + c + d, k + n));
  return fac;
}

/*
   Generate the normalization factors for Koornwinder 
   polynomials with parameter (a,b,c). Note, the 
   basis parameters appear as (a-1/2,b-1/2,c-1/2) in the 
   definintion of the polynomials, so the legendre analog
   weight of (0,0,0) corresponds to a=b=c=1/2. 

  Inputs:
    N - polynomial normalizations generated up to order N-1
    a,b,c - Koornwinder parameters
    H - Floating type T pointer to N^2 storage

  Outputs:
    H - Populated as NxN upper triangular matrix 
      - each column has normalizations for $P_i^j$,
      - where the column is j, and i <= j
      - That is, each column j corresponds to normalizations for the space of 
      - homogeneous polynomials of degree j
      - NOTE: stored in column major order 
*/   
template<typename T>
inline void sFactors (  const unsigned int N,
                        const T a, const T b, const T c, 
                        T* H  )
{
  if (!H || N == 0)
  {
    std::cerr << "scale array must not be empty";
    std::cerr << " and N cannot be 0! Exiting ...\n";
    exit(1);
  }

  T kap = (a+b+c);
  T wabc = tgamma(kap+1.5) / 
           ( tgamma(a+0.5) * tgamma(b+0.5) * tgamma(c+0.5) );

  for (unsigned int nn = 0; nn < N; ++nn)
  {
    #pragma omp simd
    for (unsigned int kk = 0; kk <= nn; ++kk)
    {
      H[kk + N*nn] = 
        sqrt( 
              wabc/( (2*nn+kap+0.5) * (2*kk+b+c) ) *
              tgamma(nn+kk+b+c+1) * tgamma(nn-kk+a+0.5) * 
              tgamma(kk+b+0.5) * tgamma(kk+c+0.5) / 
              ( 
                tgamma(nn-kk+1) * tgamma(kk+1) * 
                tgamma(nn+kk+kap+0.5) * tgamma(kk+b+c)
              )
            );
    }
  }
}


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


template void sFactors<double>  ( const unsigned int N,
                                  const double a, const double b, 
                                  const double c, double* H  );

#endif
