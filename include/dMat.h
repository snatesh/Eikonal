#ifndef _DMAT_H
#define _DMAT_H
#include<iostream>

template<typename T>
inline void dMat(T a, T b, T c, const T* H, const T* H1, unsigned int n, unsigned int mode, T* D)
{
  a = a - 0.5; b = b - 0.5; c = c - 0.5;
  unsigned int N = static_cast<unsigned int>(0.5 * (n + 1) * (n + 2)); 
  // Dx
  if (mode == 0)
  {
    D[N] = (a + b + c + 3.0) * H1[0] / H[n+2];
    D[2*N] = (1.0 + b) * (b + c + 3.0) * H1[0] / ( (b + c + 3.0) * H[1 + (n+2)] );
    unsigned int row, col; row = 4; col = 2;
    for (unsigned int blk = 2; blk <= n; ++blk)
    {
      T cPnm1k[blk+1];
      T cPnm1km1[blk+2];
      for (unsigned int kk = 0; kk <= blk; ++kk)
      {
        cPnm1k[kk] = 
          H1[kk + (n+2)*(blk-1)]  * (blk + kk + a + b + c + 2.0) * (kk + b + c + 1.0) / 
          ( H[kk + (n+2)*blk] * (2.0 * kk + b + c + 1.0) );
        cPnm1km1[0] = 0.0;
        cPnm1km1[kk+1] = 
          H1[kk + (n+2)*(blk-1)]  * (kk + 1.0 + b) * (blk + kk + b + c + 2.0) / 
          ( H[kk+1 + (n+2)*blk] * (2.0 * kk + b + c + 3.0) );
      }
      cPnm1k[blk] = 0.0;
      for (unsigned int j = 0; j <= blk; ++j)
      {
        // write it as transpose
        D[col+j-1 + N*(row+j-1)] = cPnm1k[j];
        D[col+j-2 + N*(row+j-1)] = cPnm1km1[j];
      } 
      row = row+blk+1;
      col = col+blk;
    }
  }
  // Dy
  else if (mode == 1)
  {
    D[2*N] = (b + c + 2.0) * H1[0] / H[1 + (n+2)];
    unsigned int row, col; row = 4; col = 2;
    for (unsigned int blk = 2; blk <= n; ++blk)
    {
      T cPnm1km1[blk+2];
      for (unsigned int kk = 0; kk <= blk; ++kk)
      {
        cPnm1km1[0] = 0.0;
        cPnm1km1[kk+1] = 
          H1[kk + (n+2)*(blk-1)]  * (kk + b + c + 2.0) / H[kk+1 + (n+2)*blk];
      }
      for (unsigned int j = 0; j <= blk; ++j)
      {
        // write it as transpose
        D[col+j-2 + N*(row+j-1)] = cPnm1km1[j];
      } 
      row = row+blk+1;
      col = col+blk;
    }
  }
  else
  {
    std::cerr << "invalid mode (0,1)" << std::endl;
    exit(1);
  }
}

template void dMat<double>  ( double a, double b, double c, 
                              const double* H, const double* H1, 
                              unsigned int n, unsigned int mode, 
                              double* D );

#endif
