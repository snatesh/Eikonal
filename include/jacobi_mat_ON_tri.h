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
              (2 * kk1 + b + c) * (2 * kk1 + b + c - 1) 
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
        (H[kk+1 + nn*(nn+3)] / H[kk + nn*(n+3)]) * (nn + kk + kap + 0.5) *
        (nn + kk + kap + 1.5) * (kk + 1) * (kk + b + c) /
        ( (2 * nn + kap + 0.5) * (2 * nn + kap + 1.5) * (2 * kk + b + c) * 
          (2 * kk + b + c + 1.0) );

      G[kk + nn*(n+1)] = 
        (-2.0 * H[kk+1 + nn*(n+3)] / (H[kk + nn*(n+3)])) *
        (nn - kk + a - 0.5) * (nn + kk + kap + 0.5) * (kk + 1) * (kk + b + c) /
        ( (2 * nn + kap - 0.5) * (2 * nn + kap + 1.5) * (2 * kk + b + c) * (2 * kk + b + c + 1.0));
    }
  }
 
  unsigned int inds[2] = {1,1};
  unsigned int inds1[2];
  T *A1, *A2, *B1, *B2;
  for (unsigned int nn = 1; nn <= (n-1); ++nn)
  {
    A1 = (T*) calloc(nn*(nn+1), sizeof(T));
    A2 = (T*) calloc(nn*(nn+1), sizeof(T));
    B1 = (T*) calloc(nn*nn, sizeof(T));
    B2 = (T*) calloc(nn*nn, sizeof(T));

    for (unsigned int ii = 0; ii <= nn; ++ii) 
    { 
      A1[ii + ii*nn] = A[ii + nn*(n+1)]; 
      B1[ii + ii*nn] = B[ii + nn*(n+1)];
      A2[ii + ii*nn] = E[ii + nn*(n+1)];
      A2[ii + (ii+1)*nn] = D[ii + nn*(n+1)];
      A2[ii+1 + ii*nn] = C[ii + nn*(n+1)];
      B2[ii + ii*nn] = F[ii + nn*(n+1)];
      B2[ii + (ii+1)*nn] = G[ii + nn*(n+1)];
      B2[ii+1 + ii*nn] = G[ii + nn*(n+1)];
    }

    for (unsigned int row = inds[0]; row <= inds[1]; ++ row)
    {
      for (unsigned int col = inds[0]; col <= inds[1]; ++ col)
      {
        Jn1[row + N*col] = B1[row + nn*col]; 
        Jn2[row + N*col] = B2[row + nn*col]; 
      }
    }
    inds1[0] = inds[0] + nn;
    inds1[1] = inds[1] + nn + 1;
    for (unsigned int row = inds[0]; row <= inds[1]; ++row)
    {
      for (unsigned int col = inds1[0]; col <= inds1[1]; ++col)
      {
        Jn1[row + N*col] = A1[row + nn*col];
        Jn2[row + N*col] = A2[row + nn*col];
      }
    }
    
    for (unsigned int row = inds1[0]; row <= inds1[1]; ++row)
    {
      for (unsigned int col = inds[0]; col <= inds[1]; ++col)
      {
        Jn1[row + N*col] = A1[col + (nn+1)*row];
        Jn2[row + N*col] = A2[col + (nn+1)*row];
      }
    }
    inds[0] += nn;
    inds[1] += nn + 1;
    
    free(A1);
    free(B1);
    free(A2);
    free(B1);
  } 

  for (unsigned int row = inds[0]; row <= inds[1]; ++row)
  {
    for (unsigned int col = inds1[0]; col <= inds1[1]; ++col)
    {
      Jn1[row + N*col] = B1[row + n*col];
      Jn2[row + N*col] = B2[row + n*col];
    }
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
