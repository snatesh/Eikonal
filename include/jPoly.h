#ifndef _JPOLY_H
#define _JPOLY_H
#include<cstddef>
#include<iostream>
#include<cmath>
#include<omp.h>
#include<sFactors.h>



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
template<typename T> inline void jpoly  ( T* x, unsigned int Nx, unsigned int Np, 
                                          const T a, const T b, T* V  )
{
  T apb = a + b; 
  T aa  = a * a; 
  T bb  = b * b;
  Np = Np-1;
  #pragma omp simd
  for (unsigned int i = 0; i < Nx; ++i) 
  { 
    V[i] = 1.0; 
  }
  
  if (Np > 0)
  {
    T* v = &(V[Nx]);
    #pragma omp simd
    for (unsigned int i = 0; i < Nx; ++i)
    {
      v[i] = 0.5 * ( 2.0 * (a + 1.0) + (apb + 2.0) * ( x[i] - 1) );
    }
  }
  
  T k2, k2apb, q1, q2, q3, q4;

  for (unsigned int k = 2; k <= Np; ++k)
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

/*
  Generate the normalized Koornwinder (a,b,c) Vandermonde matrix
  given an arrays X,Y of Nx points and number fo polynomial degrees n.
  The first N = dim(P_n) = [(n + d) choose n] polynomials are evaluated
  at the points (X,Y) up to total degree n-1. 
  The Nx by Np (number of poly (sum_{k=0}^n dim(P_k))) Vandermonde 
  matrix is stored in data member V.
*/
template<typename T>
class jPoly
{
  public:
    unsigned int Nx, Np, n; 
    T a, b, c;
    const T *X, *Y;
    T *H, *V;
    unsigned int nthreads;

    void init  ()
    {
      this->Np = static_cast<unsigned int>(0.5 * (n + 1) * (n + 2));
      this->V = (double*) calloc(Nx*Np, sizeof(double));
      this->H = (double*) calloc((n+1)*(n+1), sizeof(double));
      sFactors(n+1, a, b, c, this->H);
      this->ydx  = (T*) calloc(Nx, sizeof(T));
      this->x2m1 = (T*) calloc(Nx, sizeof(T));
      this->mxp1 = (T*) calloc(Nx, sizeof(T));
      this->Pk = (T*) calloc(Nx*(n+1), sizeof(T));
      this->Pnmk = (T*) calloc(Nx*(n+1)*(n+1), sizeof(T));
    }
    
    void computeV(const T* _X, const T* _Y)
    {
      this->X = _X; this->Y = _Y;
      #pragma omp parallel num_threads(nthreads)
      {
        #pragma omp for
        for (unsigned int i = 0; i < Nx; ++i)
        {
          ydx[i]  = 2.0 * Y[i] / (1 - X[i]) - 1;
          x2m1[i] = 2.0 * X[i] - 1;
          mxp1[i] = 1.0 - X[i];
        }

        #pragma omp for
        for (unsigned int kk = 0; kk <= n; ++kk) 
        {
          T* pnmk = &Pnmk[Nx*(n+1)*kk];
          jpoly<T>(x2m1, Nx, n + 1, 2.0 * kk + b + c, a - 0.5, pnmk);
        }
      }  
      jpoly<T>(this->ydx, Nx, n + 1, c - 0.5, b - 0.5, this->Pk);
      
      unsigned int ind = 1;
      for (unsigned int nn = 0; nn <= n; ++nn)
      {
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (unsigned int kk = 0; kk <= nn; ++kk)
        {
          T* v = &V[Nx*(ind+kk-1)];
          T* pnmk = &Pnmk[Nx*(n+1)*kk];
          #pragma omp simd
          for (unsigned int i = 0; i < Nx; ++i)
          {
            v[i] = 1.0 / H[kk + (n+1)*nn] * pnmk[i + Nx*(nn-kk)] * 
                   std::pow(mxp1[i],kk) * Pk[i + Nx*kk];
          }
        }
        ind = ind + nn + 1;
      }
    }

    jPoly ( unsigned int _Nx, unsigned int _n,
            T _a, T _b, T _c, unsigned int _nthreads  )
      : Nx(_Nx), n(_n), 
        a(_a), b(_b), c(_c), 
        nthreads(_nthreads) { init(); }



    jPoly ( const T* _X, const T* _Y, 
            unsigned int _Nx, 
            unsigned int _n,
            T _a, T _b, T _c, 
            unsigned int _nthreads  )
      : X(_X), Y(_Y), Nx(_Nx), 
        a(_a), b(_b), c(_c), 
        nthreads(_nthreads) 
    {
      init();
      computeV(this->X, this->Y);
    }

    ~jPoly  ()
    {
      free(H); free(V);
      free(ydx); free(x2m1); free(mxp1);
      free(Pk); free(Pnmk); 
    }

  private: 
    T *ydx, *x2m1, *mxp1;
    T *Pk, *Pnmk;

};

template void jpoly<double> ( double* x, unsigned int Nx, unsigned int Np, 
                              const double a, const double b, double* V );
template class jPoly<double>;


#endif
