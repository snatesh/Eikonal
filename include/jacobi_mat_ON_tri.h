#ifndef _JACOBI_MAT_ON_TRI_H
#define _JAOCBI_MAT_ON_TRI_H
#include"structure_factors.h"

template<typename T>
inline void jacobi_mat_ON_tri(unsigned int n, T a, T b, T c, T* Jn1, T* Jn2)
{
  unsigned int N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
  T kap = abs(a + b + c);
  T *H, *A, *B, *C, *D, *E, *F, *G;
  H = (T*) calloc((n+3)*(n+3), sizeof(T));
  A = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  B = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  C = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  D = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  E = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  F = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  G = (T*) calloc((n+1)*(n+1), sizeof(T)); 

  structure_factors_tri<T>(n+3, a, b, c, H);

  for (unsigned int nn = 0; nn <= n; ++nn)
  {
    for (unsigned int kk = 0; kk <= nn; ++kk)
    {
      A[kk + nn*(n+1)] = 
        (H[kk + (nn+1)*(n+3)] / H[kk + nn*(n+3)]) * (nn - kk + 1) * (nn + kk + kap + 0.5) /
        ( (2 * nn + kap + 0.5) * (2 * nn + kap + 1.5) );
       
      B[kk + nn*(n+1)] = 
        0.5 - ( pow(b + c + 2 * kk, 2) - pow(a-0.5, 2) ) /
        ( 2 * (2 * nn + kap - 0.5) * (2 * nn + kap + 1.5) );
      
      if (nn > 0)
      {
        for (unsigned int kk1 = 1; kk1 <= nn; ++kk1)
        {
          C[kk1-1 + nn*(n+1)] = 
            ( (H[kk1-1 + (nn+1)*(n+3)] / H[kk1 + nn*(n+3)]) * (nn - kk1 + 1) * 
              (nn - kk1 + 2) * (kk1 + b - 0.5) * (kk1 + c -0.5) 
            ) /
            ( (2 * nn + kap + 0.5) * (2 * nn + kap + 1.5) * 
              (2 * kk1 + b + c) * (2 * kk1 + b + c - 1.0) 
            );          
        }
      }
      
      E[kk + nn*(n+1)] = 
        -( 1 + ( pow(b - 0.5, 2) - pow(c-1/2,2)) /
          ( (2 * kk + b + c + 1.0) * (2 * kk + b + c -1.0) ) 
         ) * A[kk + nn*(n+1)] / 2.0;

      F[kk + nn*(n+1)] = 
        (1 + ( pow(b - 0.5, 2) - pow(c - 0.5, 2) ) /
          ( (2 * kk + b + c +1.0) * (2 * kk + b + c -1.0) ) 
        ) * (1.0 - B[kk + nn*(n+1)]) / 2.0;  
     
      if (a == 0.5 && b == 0.5 && c == 0.5)
      {
        E[nn*(n+1)] = -A[nn*(n+1)] / 2.0;
        F[nn*(n+1)] = (1.0 - B[nn*(n+1)]) / 2.0;
      } 
      
      D[kk + nn*(n+1)] = 
        (H[kk+1 + (nn+1)*(nn+3)] / H[kk + nn*(n+3)]) * (nn + kk + kap + 0.5) *
        (nn + kk + kap + 1.5) * (kk + 1) * (kk + b + c) /
        ( (2 * nn + kap + 0.5) * (2 * nn + kap + 1.5) * (2 * kk + b + c) * 
          (2 * kk + b + c + 1.0) );

      G[kk + nn*(n+1)] = 
        (-2.0 * H[kk+1 + nn*(n+3)] / (H[kk + nn*(n+3)])) *
        (nn - kk + a - 0.5) * (nn + kk + kap + 0.5) * (kk + 1) * (kk + b + c) /
        ( (2 * nn + kap - 0.5) * (2 * nn + kap + 1.5) * (2 * kk + b + c) * (2 * kk + b + c + 1.0));
    }
  }
// TODO: Figure out what's wrong with Jn2 
  unsigned int inds[2] = {1,1};
  unsigned int inds1[2];
  for (unsigned int nn = 0; nn <= (n-1); ++nn)
  {
    T* A1 = (T*) calloc((nn+1)*(nn+2), sizeof(T));
    T* A2 = (T*) calloc((nn+1)*(nn+2), sizeof(T));
    T* B1 = (T*) calloc((nn+1)*(nn+1), sizeof(T));
    T* B2 = (T*) calloc((nn+1)*(nn+1), sizeof(T));
    for (unsigned int ii = 0; ii <= nn; ++ii) 
    { 
      A1[ii + ii*(nn+1)]      = A[ii + nn*(n+1)]; 
      B1[ii + ii*(nn+1)]      = B[ii + nn*(n+1)];
      A2[ii + ii*(nn+1)]      = E[ii + nn*(n+1)];
      A2[ii + (ii+1)*(nn+1)]  = D[ii + nn*(n+1)];
      B2[ii + ii*(nn+1)]      = F[ii + nn*(n+1)];
      if (ii < nn)
      {
        A2[ii+1 + ii*(nn+1)]  = C[ii + nn*(n+1)];
        B2[ii + (ii+1)*(nn+1)]= G[ii + nn*(n+1)];
        B2[ii+1 + ii*(nn+1)]  = G[ii + nn*(n+1)];
      }
    }
    if (nn < (n-1))
    {
      unsigned int i = 0;
      for (unsigned int col = inds[0]; col <= inds[1]; ++col)
      {
        for (unsigned int row = inds[0]; row <= inds[1]; ++row)
        {
          Jn1[(row-1) + N*(col-1)] = B1[i]; 
          Jn2[(row-1) + N*(col-1)] = B2[i];
          i += 1; 
        }
      }
      inds1[0] = inds[0] + nn + 1;
      inds1[1] = inds[1] + nn + 2;
      i = 0;
      for (unsigned int col = inds1[0]; col <= inds1[1]; ++col)
      {
        for (unsigned int row = inds[0]; row <= inds[1]; ++row)
        {
          Jn1[(row-1) + N*(col-1)] = A1[i];
          Jn1[(col-1) + N*(row-1)] = A1[i];
          Jn2[(row-1) + N*(col-1)] = A2[i];
          Jn2[(col-1) + N*(row-1)] = A2[i];
          i += 1;
        }
      }
      inds[0] += nn + 1;
      inds[1] += nn + 2;
    }
    else if (nn == (n-1))
    { 
      unsigned int i = 0;
      for (unsigned int col = inds[0]; col <= inds[1]; ++col)
      {
        for (unsigned int row = inds[0]; row <= inds[1]; ++row)
        {
          Jn1[(row-1) + N*(col-1)] = B1[i];
          Jn2[(row-1) + N*(col-1)] = B2[i];
          i += 1;
        }
      }
    }
    free(A1);
    free(B1);
    free(A2);
    free(B2);
  } 

  free(H);
  free(A); 
  free(B); 
  free(C); 
  free(D); 
  free(E); 
  free(F); 
  free(G); 
}
#endif
