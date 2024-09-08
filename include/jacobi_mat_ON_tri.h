#ifndef _JACOBI_MAT_ON_TRI_H
#define _JAOCBI_MAT_ON_TRI_H
#include<structure_factors.h>
#include<iomanip>

/* TODO: use jacobi_mat_ON_tri below as constructor later */
template <typename T> 
struct Jn 
{
  unsigned int n; 
  unsigned int N; 
  T *H;
  T *A, *B, *C, *D, *E, *F, *G;
  Jn(unsigned int _n): n(_n)
  {
    N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
  }
}; 



template<typename T>
inline void jacobi_mat_ON_tri(unsigned int n, T a, T b, T c, T* Jn1, T* Jn2)
{
  unsigned int N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
  T kap = abs(a + b + c);

  // normalization constants
  T *H = (T*) calloc((n+3)*(n+3), sizeof(T));
  structure_factors_tri<T>(n+3, a, b, c, H);
  /* temporary upper triangular matrices for
     constructing coefficient matrices in recurrence */
  T *A, *B, *C, *D, *E, *F, *G;
  A = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  B = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  C = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  D = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  E = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  F = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  G = (T*) calloc((n+1)*(n+1), sizeof(T)); 
  /* block array of coefficient matrices in 3-term recurrence*/
  T** A1 = (T**) malloc((n+1) * sizeof(T**));
  T** A2 = (T**) malloc((n+1) * sizeof(T**));
  T** B1 = (T**) malloc((n+1) * sizeof(T**));
  T** B2 = (T**) malloc((n+1) * sizeof(T**));

  // populated tmp up tri mats
  for (unsigned int nn = 0; nn <= n; ++nn)
  {
    for (unsigned int kk = 0; kk <= nn; ++kk)
    {
      A[kk + nn*(n+1)] = 
        (H[kk + (nn+1)*(n+3)] / H[kk + nn*(n+3)]) * (nn - kk + 1.0) * (nn + kk + kap + 0.5) /
        ( (2.0 * nn + kap + 0.5) * (2.0 * nn + kap + 1.5) );
       
      B[kk + nn*(n+1)] = 
        0.5 - ( pow(b + c + 2.0 * kk, 2.0) - pow(a - 0.5, 2.0) ) /
        ( 2.0 * (2.0 * nn + kap - 0.5) * (2.0 * nn + kap + 1.5) );
      
      if (nn > 0)
      {
        for (unsigned int kk1 = 1; kk1 <= nn; ++kk1)
        {
          C[kk1-1 + nn*(n+1)] = 
            ( (H[kk1-1 + (nn+1)*(n+3)] / H[kk1 + nn*(n+3)]) * (nn - kk1 + 1.0) * 
              (nn - kk1 + 2.0) * (kk1 + b - 0.5) * (kk1 + c - 0.5) ) /
            ( (2.0 * nn + kap + 0.5) * (2.0 * nn + kap + 1.5) * 
              (2.0 * kk1 + b + c) * (2.0 * kk1 + b + c - 1.0) );          
        }
      }
      
      E[kk + nn*(n+1)] = 
        -( 1.0 + ( pow(b - 0.5, 2.0) - pow(c - 0.5, 2.0) ) /
          ( (2.0 * kk + b + c + 1.0) * (2.0 * kk + b + c - 1.0) ) ) * A[kk + nn*(n+1)] / 2.0;

      F[kk + nn*(n+1)] = 
        ( 1.0 + ( pow(b - 0.5, 2.0) - pow(c - 0.5, 2.0) ) /
          ( (2.0 * kk + b + c + 1.0) * (2.0 * kk + b + c - 1.0) ) ) * (1.0 - B[kk + nn*(n+1)]) / 2.0;  
     
      if (a == 0.5 && b == 0.5 && c == 0.5)
      {
        E[nn*(n+1)] = -A[nn*(n+1)] / 2.0;
        F[nn*(n+1)] = (1.0 - B[nn*(n+1)]) / 2.0;
      } 
      
      D[kk + nn*(n+1)] = 
        (H[kk+1 + (nn+1)*(n+3)] / H[kk + nn*(n+3)]) * (nn + kk + kap + 0.5) *
        (nn + kk + kap + 1.5) * (kk + 1.0) * (kk + b + c) /
        ( (2.0 * nn + kap + 0.5) * (2.0 * nn + kap + 1.5) * (2.0 * kk + b + c) * 
          (2.0 * kk + b + c + 1.0) );

      G[kk + nn*(n+1)] = 
        (-2.0 * H[kk+1 + nn*(n+3)] / (H[kk + nn*(n+3)]) ) *
        (nn - kk + a - 0.5) * (nn + kk + kap + 0.5) * (kk + 1) * (kk + b + c) /
        ( (2.0 * nn + kap - 0.5) * (2.0 * nn + kap + 1.5) * (2.0 * kk + b + c) * (2.0 * kk + b + c + 1.0));
    }
  }
  // prepoulate blocks for easy assignment into Jn_i
  // efficiency/redundancy doesn't matter 
  // - amortize out precompute, go for readablity instead
  for (unsigned int nn = 0; nn <= (n-1); ++nn)
  {
    A1[nn]    = (T*) calloc((nn+1)*(nn+2), sizeof(T));
    A2[nn]    = (T*) calloc((nn+1)*(nn+2), sizeof(T));
    B1[nn]    = (T*) calloc((nn+1)*(nn+1), sizeof(T));
    B2[nn]    = (T*) calloc((nn+1)*(nn+1), sizeof(T));
    for (unsigned int ii = 0; ii <= nn; ++ii) 
    { 
      A1[nn][ii + ii*(nn+1)]        = A[ii + nn*(n+1)];
      B1[nn][ii + ii*(nn+1)]        = B[ii + nn*(n+1)];
      A2[nn][ii + ii*(nn+1)]        = E[ii + nn*(n+1)];
      A2[nn][ii + (ii+1)*(nn+1)]    = D[ii + nn*(n+1)];
      B2[nn][ii + ii*(nn+1)]        = F[ii + nn*(n+1)];
      if (ii < nn)
      {
        A2[nn][ii+1 + ii*(nn+1)]      = C[ii + nn*(n+1)];
        B2[nn][ii + (ii+1)*(nn+1)]    = G[ii + nn*(n+1)];
        B2[nn][ii+1 + ii*(nn+1)]      = G[ii + nn*(n+1)];
      }
    }
  }
  // copy into Jn_i
  unsigned int inds[2] = {1,1};
  unsigned int inds1[2];
  unsigned int i;
  for (unsigned int nn = 1; nn <= (n-1); ++nn)
  {
    i = 0;
    for (unsigned int col = inds[0]; col <= inds[1]; ++col)
    {
      for (unsigned int row = inds[0]; row <= inds[1]; ++row)
      {
        Jn1[(row-1) + N*(col-1)] = B1[nn-1][i]; 
        Jn2[(row-1) + N*(col-1)] = B2[nn-1][i];
        i += 1; 
      }
    }
    inds1[0] = inds[0] + nn;
    inds1[1] = inds[1] + nn + 1;
    i = 0;
    for (unsigned int col = inds1[0]; col <= inds1[1]; ++col)
    {
      for (unsigned int row = inds[0]; row <= inds[1]; ++row)
      {
        Jn1[(row-1) + N*(col-1)] = A1[nn-1][i];
        Jn1[(col-1) + N*(row-1)] = A1[nn-1][i];
        Jn2[(row-1) + N*(col-1)] = A2[nn-1][i];
        Jn2[(col-1) + N*(row-1)] = A2[nn-1][i];
        
        i += 1;
      }
    }
    inds[0] += nn;
    inds[1] += nn + 1;
  }
  // get last block
  i = 0;
  for (unsigned int col = inds[0]; col <= inds[1]; ++col)
  {
    for (unsigned int row = inds[0]; row <= inds[1]; ++row)
    {
      Jn1[(row-1) + N*(col-1)] = B1[n-1][i]; 
      Jn2[(row-1) + N*(col-1)] = B2[n-1][i];
      i += 1;
    }
  }
  
  for (unsigned int nn = 0; nn < n; ++nn)
  {
    free(A1[nn]); free(A2[nn]);
    free(B1[nn]); free(B2[nn]);
  }
  free(A1); free(B1); 
  free(A2); free(B2); 
  free(H); free(A); 
  free(B); free(C); 
  free(D); free(E); 
  free(F); free(G); 

}
   
#endif
