#ifndef _PROMOTION_MAT_TRI_H 
#define _PROMOTION_MAT_TRI_H
#include<iostream>

template<typename T>
inline void promotion_mat_tri(T a, T b, T c, const T* H, const T* H1, unsigned int n, unsigned int mode, T* K)
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

  }
  // c + 1
  else if (mode == 2)
  {

  }
  else
  {
    std::cerr << "invalid mode (0,1,2)" << std::endl; 
  }

}

//function K = promotion_mat_tri(a,b,c,H,H1,mode)
//n = size(H,1)-2;
//a = a-1/2; b = b-1/2; c = c-1/2; 
//N = (n+1)*(n+2)/2;
//K = zeros(N); 
//% a+1
//if mode == 0
//  K(1,1) = 1; 
//  row = 2; col = 1;
//  for blk = 1:n  
//    kk = 0:blk;
//    cPnk = H1(kk+1,blk+1)'.*(blk+kk+a+b+c+2)./(H(kk+1,blk+1)'.*(2*blk+a+b+c+2));
//    cPnm1k = H1(kk+1,blk)'.*(blk+kk+b+c+1)./(H(kk+1,blk+1)'.*(2*blk+a+b+c+2));
//    cPnm1k(blk+1) = 0; 
//    for j = 0:blk
//      K(row+j,col+j) = cPnm1k(j+1);
//      K(row+j,col+j+blk) = cPnk(j+1);
//    end
//    row = row+blk+1;
//    col = col+blk;
//  end
//% b+1
//elseif mode == 1
//  K(1,1) = 1; 
//  K(2,1) = -(1+a)*(b+c+1)*H1(1,1)/((2+a+b+c+2)*(b+c+1)*H(1,2)); 
//  K(2,2) = (1+a+b+c+2)*(b+c+1)*H1(1,2)/((2+a+b+c+2)*(b+c+1)*H(1,2));
//  K(3,1) = (1+c)*(2+b+c+1)*H1(1,1)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
//  K(3,2) = -(1+c)*H1(1,2)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
//  K(3,3) = (2+a+b+c+2)*(1+b+c+1)*H1(2,2)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
//  row = 4; col = 2;
//  for blk = 2:n  
//    kk = 0:blk;
//    cPnk = H1(kk+1,blk+1)'.*(blk+kk+a+b+c+2).*(kk+b+c+1) ./ ...
//      (H(kk+1,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+b+c+1));
//    
//    cPnm1k = -H1(kk+1,blk)'.*(blk-kk+a).*(kk+b+c+1) ./ ... 
//      (H(kk+1,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+b+c+1));
//    
//    cPnkm1 = [0,-H1(kk+1,blk+1)'.*(kk+1+c).*(blk-kk-1+1) ./ ...
//      (H(kk+2,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+2+b+c+1))];
//    
//    cPnm1km1 = [0,H1(kk+1,blk)'.*(kk+1+c).*(blk+kk+1+b+c+1) ./ ... 
//      (H(kk+2,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+2+b+c+1))];
//    
//    cPnm1k(blk+1) = 0; 
//    for j = 0:blk
//      K(row+j,col+j) = cPnm1k(j+1);
//      K(row+j,col+j+blk) = cPnk(j+1);
//      K(row+j,col+j-1) = cPnm1km1(j+1);
//      K(row+j,col+j+blk-1) = cPnkm1(j+1);
//    end
//    row = row+blk+1;
//    col = col+blk;
//  end 
//% c+1
//elseif mode == 2
//  K(1,1) = 1; 
//  K(2,1) = -(1+a)*(b+c+1)*H1(1,1)/((2+a+b+c+2)*(b+c+1)*H(1,2)); 
//  K(2,2) = (1+a+b+c+2)*(b+c+1)*H1(1,2)/((2+a+b+c+2)*(b+c+1)*H(1,2));
//  K(3,1) = -(1+b)*(2+b+c+1)*H1(1,1)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
//  K(3,2) = (1+b)*H1(1,2)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
//  K(3,3) = (2+a+b+c+2)*(1+b+c+1)*H1(2,2)/((2+a+b+c+2)*(2+b+c+1)*H(2,2));
//  row = 4; col = 2;
//  for blk = 2:n  
//    kk = 0:blk;
//    cPnk = H1(kk+1,blk+1)'.*(blk+kk+a+b+c+2).*(kk+b+c+1) ./ ...
//      (H(kk+1,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+b+c+1));
//    
//    cPnm1k = -H1(kk+1,blk)'.*(blk-kk+a).*(kk+b+c+1) ./ ... 
//      (H(kk+1,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+b+c+1));
//    
//    cPnkm1 = [0,H1(kk+1,blk+1)'.*(kk+1+b).*(blk-kk-1+1) ./ ...
//      (H(kk+2,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+2+b+c+1))];
//    
//    cPnm1km1 = [0,-H1(kk+1,blk)'.*(kk+1+b).*(blk+kk+1+b+c+1) ./ ... 
//      (H(kk+2,blk+1)'.*(2*blk+a+b+c+2).*(2*kk+2+b+c+1))];
//    
//    cPnm1k(blk+1) = 0; 
//    for j = 0:blk
//      K(row+j,col+j) = cPnm1k(j+1);
//      K(row+j,col+j+blk) = cPnk(j+1);
//      K(row+j,col+j-1) = cPnm1km1(j+1);
//      K(row+j,col+j+blk-1) = cPnkm1(j+1);
//    end
//    row = row+blk+1;
//    col = col+blk;
//  end 
//else
//  error('invalid mode (0,1,2)');
//end
//% we write as transpose 
//% so P_{a+1,b,c}^T K = P_{a,b,c}
//K = K';
//end  

#endif
