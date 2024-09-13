#ifndef _KMAT_H 
#define _KMAT_H
#include<iostream>

template<typename T>
inline void kMat  ( T a, T b, T c, const T* H,
                    const T* H1, unsigned int n, 
                    unsigned int mode, T* K )
{
  a = a - 0.5; b = b - 0.5; c = c - 0.5;
  unsigned int N = static_cast<unsigned int>(0.5 * (n + 1) * (n + 2)); 
  // a + 1
  if (mode == 0)
  {
    K[0] = 1.0; 
    unsigned int row, col; row = 2; col = 1;
    for (unsigned int blk = 1; blk <= n; ++blk)
    {
      T cPnk[blk+1];
      T cPnm1k[blk+1];
      for (unsigned int kk = 0; kk <= blk; ++kk)
      {
        cPnk[kk] = 
          H1[kk + (n+2)*blk] * (blk + kk + a + b + c + 2.0) / 
          ( H[kk + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) );
        cPnm1k[kk] = 
          H1[kk + (n+2)*(blk-1)]  * (blk + kk + b + c + 1.0) / 
          ( H[kk + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) );
      }
      cPnm1k[blk] = 0;
      for (unsigned int j = 0; j <= blk; ++j)
      {
        // write it as transpose
        K[col+j-1 + N*(row+j-1)] = cPnm1k[j];
        K[col+j+blk-1 + N*(row+j-1)] = cPnk[j];
      } 
      row = row+blk+1;
      col = col+blk;
    }
   
  }
  // b + 1
  else if (mode == 1)
  {
    // writing as transpose
    K[0] = 1.0;

    K[N] = 
      -(1.0 + a) * (b + c + 1.0) * H1[0] / 
      ( (2.0 + a + b + c + 2.0) * (b + c + 1.0) * H[n+2] );
 
    K[1 + N] =  
      (1.0 + a + b + c + 2.0) * (b + c + 1.0) * H1[n+2] /
      ( (2.0 + a + b + c + 2.0) * (b + c + 1.0) * H[n+2] );
 
    K[2*N] = 
      (1.0 + c) * (2.0 + b + c + 1.0) * H1[0] /
      ( (2.0 + a + b + c + 2.0) * (2.0 + b + c + 1.0) * H[1 + (n+2)]);
 
    K[1 + 2*N] = 
      -(1.0 + c) * H1[n+2] / 
      ( (2.0 + a + b + c + 2.0) * (2.0 + b + c + 1.0) * H[1 + (n+2)] );
 
    K[2 + 2*N] = 
      (2.0 + a + b + c + 2.0) * (1.0 + b + c + 1.0) * H1[1 + (n+2)] / 
      ( (2.0 + a + b + c + 2.0) * (2.0 + b + c + 1.0) * H[1 + (n+2)] );
    
    unsigned int row, col; row = 4; col = 2;
    for (unsigned int blk = 2; blk <= n; ++blk)
    {
      T cPnk[blk+1];
      T cPnm1k[blk+1];
      T cPnkm1[blk+2];
      T cPnm1km1[blk+2];
      for (unsigned int kk = 0; kk <= blk; ++kk)
      {
        cPnk[kk] = 
          H1[kk + (n+2)*blk] * (blk + kk + a + b + c + 2.0) * (kk + b + c + 1.0) / 
          ( H[kk + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) * (2.0 * kk + b + c + 1.0) );
        cPnm1k[kk] = 
          -H1[kk + (n+2)*(blk-1)]  * (blk - kk + a) * (kk + b + c + 1.0) / 
          ( H[kk + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) * (2.0 * kk + b + c + 1.0) );
        cPnkm1[0] = 0.0;
        cPnkm1[kk+1] = 
          -H1[kk + (n+2)*blk] * (kk + 1.0 + c) * (blk - kk) / 
          ( H[kk+1 + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) * (2.0 * kk + b + c + 3.0) );
        cPnm1km1[0] = 0.0;
        cPnm1km1[kk+1] = 
          H1[kk + (n+2)*(blk-1)]  * (kk + 1.0 + c) * (blk + kk + b + c + 2.0) / 
          ( H[kk+1 + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) * (2.0 * kk + b + c + 3.0) );
      }
      cPnm1k[blk] = 0.0;
      for (unsigned int j = 0; j <= blk; ++j)
      {
        // write it as transpose
        K[col+j-1 + N*(row+j-1)] = cPnm1k[j];
        K[col+j+blk-1 + N*(row+j-1)] = cPnk[j];
        K[col+j-2 + N*(row+j-1)] = cPnm1km1[j];
        K[col+j+blk-2 + N*(row+j-1)] = cPnkm1[j];
      } 
      row = row+blk+1;
      col = col+blk;
    }
  }
  // c + 1
  else if (mode == 2)
  {
    // writing as transpose
    K[0] = 1.0;

    K[N] = 
      -(1.0 + a) * (b + c + 1.0) * H1[0] / 
      ( (2.0 + a + b + c + 2.0) * (b + c + 1.0) * H[n+2] );
 
    K[1 + N] =  
      (1.0 + a + b + c + 2.0) * (b + c + 1.0) * H1[n+2] /
      ( (2.0 + a + b + c + 2.0) * (b + c + 1.0) * H[n+2] );
 
    K[2*N] = 
      -(1.0 + b) * (2.0 + b + c + 1.0) * H1[0] /
      ( (2.0 + a + b + c + 2.0) * (2.0 + b + c + 1.0) * H[1 + (n+2)]);
 
    K[1 + 2*N] = 
      (1.0 + b) * H1[n+2] / 
      ( (2.0 + a + b + c + 2.0) * (2.0 + b + c + 1.0) * H[1 + (n+2)] );
 
    K[2 + 2*N] = 
      (2.0 + a + b + c + 2.0) * (1.0 + b + c + 1.0) * H1[1 + (n+2)] / 
      ( (2.0 + a + b + c + 2.0) * (2.0 + b + c + 1.0) * H[1 + (n+2)] );
    
    unsigned int row, col; row = 4; col = 2;
    for (unsigned int blk = 2; blk <= n; ++blk)
    {
      T cPnk[blk+1];
      T cPnm1k[blk+1];
      T cPnkm1[blk+2];
      T cPnm1km1[blk+2];
      for (unsigned int kk = 0; kk <= blk; ++kk)
      {
        cPnk[kk] = 
          H1[kk + (n+2)*blk] * (blk + kk + a + b + c + 2.0) * (kk + b + c + 1.0) / 
          ( H[kk + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) * (2.0 * kk + b + c + 1.0) );
        cPnm1k[kk] = 
          -H1[kk + (n+2)*(blk-1)]  * (blk - kk + a) * (kk + b + c + 1.0) / 
          ( H[kk + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) * (2.0 * kk + b + c + 1.0) );
        cPnkm1[0] = 0.0;
        cPnkm1[kk+1] = 
          H1[kk + (n+2)*blk] * (kk + 1.0 + b) * (blk - kk) / 
          ( H[kk+1 + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) * (2.0 * kk + b + c + 3.0) );
        cPnm1km1[0] = 0.0;
        cPnm1km1[kk+1] = 
          -H1[kk + (n+2)*(blk-1)]  * (kk + 1.0 + b) * (blk + kk + b + c + 2.0) / 
          ( H[kk+1 + (n+2)*blk] * (2.0 * blk + a + b + c + 2.0) * (2.0 * kk + b + c + 3.0) );
      }
      cPnm1k[blk] = 0.0;
      for (unsigned int j = 0; j <= blk; ++j)
      {
        // write it as transpose
        K[col+j-1 + N*(row+j-1)] = cPnm1k[j];
        K[col+j+blk-1 + N*(row+j-1)] = cPnk[j];
        K[col+j-2 + N*(row+j-1)] = cPnm1km1[j];
        K[col+j+blk-2 + N*(row+j-1)] = cPnkm1[j];
      } 
      row = row+blk+1;
      col = col+blk;
    }

  }
  else
  {
    std::cerr << "invalid mode (0,1,2)" << std::endl; 
    exit(1);
  }
}

template void kMat<double>( double a, double b, double c, 
                            const double* H, const double* H1, 
                            unsigned int n, unsigned int mode, 
                            double* K );

#endif
