#ifndef _JPOLY_H
#define _JPOLY_H
#include<cstddef>
#include<iostream>

/*
  Generate the Jacobi (a,b) Vandermonde matrix
  given an array of x of Nx points and number of polynomials N.
  The first N (order=0,..,N-1) polynomials are evaluated
  on the grid x. Each column j (zero-based) corresponds to 
  the j^th order Jacobi polynomial evaluated on x. The output is an
  Nx (length of x) by Np (number of poly) Vandermonde matrix.

  Inputs:
    x   - Floating point type T pointer to Nx storage.
        - these are Nx 1D points 
    Nx  - number of points in x (length(x))
    Np  - number of polynomials 
    a,b - Jacoby poly parameters
    V   - Floating point type T pointer to Nx*Np storage 
 
  Outputs:
    V - populated as Nx x Np matrix stored in column major order
      - NOTE: these are not normalizez. One must call structure_factors()
      -       to obtain the normalization coefficients. 
*/
template<typename T>
inline void jPoly(T* x, unsigned int Nx, unsigned int N, const T a, const T b, T* V)
{
  T apb = a + b; 
  T aa  = a * a; 
  T bb  = b * b;
  N = N-1;
  #pragma omp simd
  for (unsigned int i = 0; i < Nx; ++i) 
  { 
    V[i] = 1.0; 
  }
  
  if (N > 0)
  {
    T* v = &(V[Nx]);
    #pragma omp simd
    for (unsigned int i = 0; i < Nx; ++i)
    {
      v[i] = 0.5 * ( 2.0 * (a + 1.0) + (apb + 2.0) * ( x[i] - 1) );
    }
  }
  
  T k2, k2apb, q1, q2, q3, q4;

  for (unsigned int k = 2; k <= N; ++k)
  {
    k2 = 2.0 * k;
    k2apb = k2 + apb;
    q1 =  k2 * (k + apb) * (k2apb - 2.0); 
    q2 = (k2apb - 1.0) * (aa - bb);
    q3 = (k2apb - 2.0) * (k2apb - 1.0) * k2apb;
    q4 =  2.0 * (k + a - 1.0) * (k + b - 1.0) * k2apb;
    T* vkp1 = &V[Nx*k];
    T* vk   = &V[Nx*(k-1)];
    T* vkm1 = &V[Nx*(k-2)];
    #pragma omp simd
    for (unsigned int i = 0; i < Nx; ++i)
    {
      vkp1[i] =  ( (q2 + q3 * x[i]) * vk[i] - q4 * vkm1[i] ) / q1; 
    }
  }
}

template<typename T>
inline void jPoly_tri(T* X, T* Y, T* H, unsigned int n, T a, T b, T c)
{
  // total num of polys up to degree n in d dimensions is nchoosek(n+d,n)
  unsigned int N = static_cast<unsigned int>(0.5 * (n + 1) * (n + 2));

}

#endif
