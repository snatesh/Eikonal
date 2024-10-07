#ifndef _LMAT_H
#define _LMAT_H
#include<iostream>

template<typename T>
inline void lMat(T a, T b, T c, const T* H, const T* H1, unsigned int n, unsigned int mode, T* L)
{
  a = a - 0.5; b = b - 0.5; c = c - 0.5;
  unsigned int N = static_cast<unsigned int>(0.5 * (n + 1) * (n + 2)); 
  if (mode == 0)
  {
    unsigned int row = 1, col = 1;
    for (unsigned int blk = 0; blk <= n; ++blk)
    {
      T cPnk[blk+1];
      T cPnp1k[blk+1];
      for (unsigned int kk = 0; kk <= blk; ++kk)
      {
        cPnk[kk] = 
          H1[kk + (n+2)*blk] * (blk - kk + a) /
          ( H[kk + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) );
        cPnp1k[kk] = 
          H1[kk + (n+2)*(blk+1)] * (blk - kk + 1) /
          ( H[kk + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) );
      }
      for (unsigned int j = 0; j <= blk; ++j)
      {
        if (row + j <= N)
        {
          if (col + j <= N)
          {
            L[col+j-1 + N*(row+j-1)] = cPnk[j]; 
          }
          if (col + j + blk + 1 <= N)
          {
            L[col+j+blk + N*(row+j-1)] = cPnp1k[j];
          }
        }
      }
      row = row + blk + 1;
      col = col + blk + 1;
    }
  }
  else if (mode == 1)
  {
    L[0] = b*H1[0] / ( (a + b + c + 2.0) * H[0] );
    L[1] = -b*H1[n+2] / ( (b + c + 1.0) * (a + b + c + 2.0) * H[0] );
    L[2] = H1[n+3] / ( (b + c + 1.0) * H[0] );
    unsigned int row = 2, col = 2;
    for (unsigned int blk = 1; blk <= n; ++blk)
    {
      T cPnk[blk+1];
      T cPnkp1[blk+1];
      T cPnp1k[blk+1];
      T cPnp1kp1[blk+1];
      for (unsigned int kk = 0; kk <= blk; ++kk)
      {
        cPnk[kk] = 
          H1[kk + (n+2)*blk] * (kk + b) * (blk + kk + b + c + 1.0) /
          ( H[kk + (n+2)*blk] * (2.0 * kk + b + c + 1.0) * (2.0 * blk + a + b + c + 2.0) );
        cPnkp1[kk] = 
          -H1[kk+1 + (n+2)*blk] * (kk + 1) * (blk - kk + a) /
          ( H[kk + (n+2)*blk] * (2.0 * kk + b + c + 1.0) * (2.0 * blk + a + b + c + 2.0) );
        cPnp1k[kk] = 
          -H1[kk + (n+2)*(blk+1)] * (kk + b) * (blk - kk + 1.0) / 
          ( H[kk + (n+2)*blk] * (2.0 * kk + b + c + 1.0) * (2.0 * blk + a + b + c + 2.0) );
        cPnp1kp1[kk] = 
          H1[kk+1 + (n+2)*(blk+1)] * (kk + 1.0) * (blk + kk + a + b + c + 2.0 ) / 
          ( H[kk + (n+2)*blk] * (2.0 * kk + b + c + 1.0) * (2.0 * blk + a + b + c + 2.0) ); 

      }
      cPnkp1[blk] = 0;
      
      for (unsigned int j = 0; j <= blk; ++j)
      {
        if (row + j <= N)
        {
          if (col + j <= N)
          {
            L[col+j-1 + N*(row+j-1)] = cPnk[j];
          }
          if (col + j + 1 <= N)
          {
            L[col+j + N*(row+j-1)] = cPnkp1[j];
          }
          if (col + j + blk + 1 <= N)
          {
            L[col+j+blk + N*(row+j-1)] = cPnp1k[j];
          }
          if (col + j + blk + 2 <= N)
          {
            L[col+j+blk+1 + N*(row+j-1)] = cPnp1kp1[j];
          }
        }
      }
      row = row + blk + 1;
      col = col + blk + 1; 
    }
  }
  else
  {
    std::cerr << "invalid mode (0,1)" << std::endl;
    exit(1);
  }
  

}


template void lMat<double>  ( double a, double b, double c, 
                              const double* H, const double* H1, 
                              unsigned int n, unsigned int mode, 
                              double* D);

#endif
