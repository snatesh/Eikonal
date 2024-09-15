#ifndef _JEVD_H
#define _JEVD_H
#include<cmath>
#include<omp.h>


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
  
  If J1 .. Jd have a common eigenstructure, i.e. a common set
  of orthonormal eigenvectors, then E1 .. Ed are diagonal
  matrices with eigenvalues of J1 .. Jd on the diagonals.
  
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


  // BELOW ONLY WORKS FOR n=2
*/

template<typename T>
struct jointDiag
{
  // n mxm matrices
  unsigned int m, n, nm, nthreads;
  T thresh;
  T *J, *E;
  

  jointDiag ( unsigned int _m,
              unsigned int _n,
              T _thresh, T* _J, 
              unsigned int _nthreads)
    : m(_m), n(_n), nm(_m*_n), 
      nthreads(_nthreads),
      thresh(_thresh), J(_J)
  {
    E = (T*) calloc(m*nm, sizeof(T));

    bool go = 1;
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for
      for (unsigned int i = 0; i < m*nm; ++i)
      {
        E[i] = J[i];
      }
      while (go)
      {
        T ton, toff, theta, c, s;
        #pragma omp for collapse(2) schedule(dynamic)
        for (unsigned int p = 1; p <= m-1; ++p)
        {
          for (unsigned int q = p+1; q <=m; ++q)
          {
            /* TODO: need to define custom packing for J to 
                     effectively use simd optimziations. 
                     keep dumb packing for now */
            T App[n], Aqq[n], Apq[n], Aqp[n];
            T g1[n], g2[n], g[2*n];
            #pragma omp simd
            for (unsigned int nn = 0; nn < n; ++nn)
            {
              App[nn] = J[p-1 + nm*(p-1+nn)];
              Aqq[nn] = J[q-1 + nm*(q-1+nn)];
              Apq[nn] = J[p-1 + nm*(q-1+nn)];
              Aqp[nn] = J[q-1 + nm*(p-1+nn)]; 
            }        
            #pragma omp simd
            for (unsigned int nn = 0; nn < n; ++nn)
            {
              g1[nn] = App[nn] - Aqq[nn];
              g2[nn] = Apq[nn] - Aqp[nn];
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
            if (std::abs(s) > thresh)
            {
              for (unsigned int nn = 0; nn < n; ++nn)
              {
                T Mp[m], Mq[m];
                #pragma omp simd
                for (unsigned i = 0; i < m; ++i)
                {
                  Mp[i] = J[i + nm*(p-1+nn)];
                  Mq[i] = J[i + nm*(q-1+nn)];
                }
                #pragma omp simd
                for (unsigned int i = 0; i < m; ++i)
                {
                  E[i + nm*(p-1+nn)] = c*Mp[i] + s*Mq[i];
                  E[i + nm*(q-1+nn)] = c*Mq[i] - s*Mp[i];
                } 
              }
              T rowp[nm], rowq[nm];
              #pragma omp simd
              for (unsigned int i = 0; i < nm; ++i)
              {
                rowp[i] = E[(p-1) + nm*i];
                rowq[i] = E[(q-1) + nm*i];
              }
              #pragma omp simd
              for (unsigned int i = 0; i < nm; ++i)
              {
                E[(p-1) + nm*i] = c * rowp[i] + s * rowq[i];
                E[(q-1) + nm*i] = c * rowq[i] + s * rowp[i];
              }
        
            }
          
          }
        }
      }
    }
  }

  ~jointDiag() { free(E); E = 0; }

};

template struct jointDiag<double>;

#endif
