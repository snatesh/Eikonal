#ifndef _JEVD_H
#define _JEVD_H
#include<cmath>
#include<omp.h>
#include<iostream>
#include<iomanip>
#include<cblas.h>

/*
  Joint diagonalization (possibly
  approximate) of REAL matrices.
  
  This function minimizes a joint diagonality criterion
  through n matrices of size m by m.
  
  Input :
  - J is matrix of size m x mn which is a concatenation
    n matrices of size m x m (block-column) such that
    J = [J1 J2 .. Jn]
  - threshold is the approximate diagonality measure tolerance
    (default to sqrt(mach_eps).
  Output :
  - V is an m by m orthogonal matrix - the approximate orthogonalizer
    of all J1,..,Jn
  - E is the concatenation of d approximately diagonal m by m matrices
    such that E = [E1 .. En], and J1 = V E1 V' .. Jd = V En V'
  
  If J1 .. Jn have a common eigenstructure, i.e. a common set
  of orthonormal eigenvectors, then E1 .. En are exactly diagonal
  matrices with eigenvalues of J1 .. Jn on the diagonals.
  
  The value of threshold is problem dependent, and should be 
  experimented with for your case.
  
 The algorithm is explained in:

@article{SC-siam,
   author = "Jean-Fran\c{c}ois Cardoso and Antoine Souloumiac",
   journal = "{SIAM} J. Mat. Anal. Appl.",
   title = "Jacobi angles for simultaneous diagonalization",
   pages = "161--164",
   volume = "17",
   number = "1",
   month = jan,
   year = {1995}}

*/

template<typename T>
struct jointDiag
{
  // n mxm matrices
  unsigned int m, n, nm, nthreads;
  T thresh;
  T *J, *V, *Vorth;
  bool hasV; 
  
  void checkOrth()
  {
    if (hasV)
    {
      if (!Vorth)
      {
        Vorth = (T*) calloc(m*m, sizeof(T));
        cblas_dgemm ( CblasColMajor, 
                      CblasNoTrans, 
                      CblasTrans, 
                      m, m, m, 
                      1.0, V, m, 
                      V, m, 0.0, Vorth, m );
      }
    }
    
  }
   
  void printVVt()
  {
    if (hasV)
    {
      for (unsigned int i = 0; i < m; ++i)
      {
        for (unsigned int j = 0; j < m; ++j)
        {
          std::cout << std::setw(10);
          std::cout << Vorth[i + m*j] << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
  }

  void printEigs()
  {
    for (unsigned int i = 0; i < m; ++i)
    {
      std::cout << "(" << J[i + i*m] << ", ";
      std::cout << J[i + m*(i + m)] << ")" << std::endl;
    }
    
  } 

  jointDiag ( unsigned int _m,
              unsigned int _n,
              T _thresh, T* _J, 
              unsigned int _nthreads, 
              bool _hasV = false  )
    : m(_m), n(_n), nm(_m*_n), 
      nthreads(_nthreads),
      thresh(_thresh), J(_J), hasV(_hasV)
  {

    bool go = 1;
    if (hasV) 
    { 
      V = (T*) calloc(m*m, sizeof(T)); 
      for (unsigned int i = 0; i < m; ++i)
      {
        V[i + i*m] = 1.0;
      }
    }
   
    unsigned int iter = 0; 
    while (go)
    {
      go = 0; iter += 1;

      T ton, toff, theta, c, s;

      #pragma omp parallel for collapse(2) num_threads(nthreads)
      for (unsigned int p = 1; p <= m-1; ++p)
      {
        for (unsigned int q = p+1; q <=m; ++q)
        {
          /* TODO: need to define custom packing for J to 
                   effectively use simd optimziations. 
                   keep dumb packing for now */
          T App[n], Aqq[n], Apq[n], Aqp[n];
          T g1[n], g2[n];
          #pragma omp simd
          for (unsigned int nn = 0; nn < n; ++nn)
          {
            App[nn] = J[p-1 + m*(p-1+nn*m)];
            Aqq[nn] = J[q-1 + m*(q-1+nn*m)];
            Apq[nn] = J[p-1 + m*(q-1+nn*m)];
            Aqp[nn] = J[q-1 + m*(p-1+nn*m)]; 
          }        
          #pragma omp simd
          for (unsigned int nn = 0; nn < n; ++nn)
          {
            g1[nn] = App[nn] - Aqq[nn];
            g2[nn] = Apq[nn] + Aqp[nn];
          }
          T G00 = 0, G01 = 0, G11 = 0;
          #pragma omp simd reduction(+:G00,G01,G11)
          for (unsigned int nn = 0; nn < n; ++nn) 
          { 
            G00 += g1[nn] * g1[nn];
            G11 += g2[nn] * g2[nn];
            G01 += g1[nn] * g2[nn];
          }

          ton = G00 - G11; toff = G01 * 2.0; 
          theta = 0.5 * std::atan2( toff, ton + std::sqrt(ton * ton + toff * toff) );
          c = std::cos(theta); s = std::sin(theta);
          go = ( go || (std::abs(s) > thresh) );
          //std::cout << "JEVD abs(s) = " << std::abs(s) << std::endl;
          if (std::abs(s) > thresh)
          {
            for (unsigned int nn = 0; nn < n; ++nn)
            {
              T Mp[m], Mq[m];
              #pragma omp simd
              for (unsigned i = 0; i < m; ++i)
              {
                Mp[i] = J[i + m*(p-1+nn*m)];
                Mq[i] = J[i + m*(q-1+nn*m)];
              }
              #pragma omp simd
              for (unsigned int i = 0; i < m; ++i)
              {
                J[i + m*(p-1+nn*m)] = c*Mp[i] + s*Mq[i];
                J[i + m*(q-1+nn*m)] = c*Mq[i] - s*Mp[i];
              } 
            }
            T rowp[nm], rowq[nm];
            #pragma omp simd
            for (unsigned int i = 0; i < nm; ++i)
            {
              rowp[i] = J[(p-1) + m*i];
              rowq[i] = J[(q-1) + m*i];
            }
            #pragma omp simd
            for (unsigned int i = 0; i < nm; ++i)
            {
              J[(p-1) + m*i] = c * rowp[i] + s * rowq[i];
              J[(q-1) + m*i] = c * rowq[i] - s * rowp[i];
            }
           
            if (hasV)
            { 
              T* vp = &V[(p-1)*m];
              T* vq = &V[(q-1)*m];
              T tmp;
              #pragma omp simd
              for (unsigned int i = 0; i < m; ++i)
              {
                tmp = vp[i];
                vp[i] = c * vp[i] + s * vq[i];
                vq[i] = c * vq[i] - s * tmp;
              }
            }
          }
        }
      }
    }
  }

  ~jointDiag() 
  {
    if (hasV)
    { 
      if (V) { free(V); V = 0; } 
      if (Vorth) { free(Vorth); Vorth = 0;} 
    } 
  }

};

template struct jointDiag<double>;

#endif
