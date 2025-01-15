#ifndef _JPOLY_H
#define _JPOLY_H
#include<cstddef>
#include<iostream>
#include<cmath>
#include<omp.h>
#include<sFactors.hh>
#include<cblas.h>
#include<type_traits>
#include<vector>

/* 
   Evaluate a Jacobi (a,b) polynomial at a point x in [-1,1].
   This uses the definition of the polynomials in terms of 
   the hypergoemetric 2F1 series. 
  
   This is not a stable way of evaluating the polynomials.
   We loose nearly all digits (relative to the order of
   accuracy of other methods which call this function)
   for polynomial degree n > 40, largely due to recurrent
   products in the falling factorial and hypergeometric
   functions.

   Really speaking, we shouldn't ever need to go higher
   than n = 20, in terms of problems which demand high order, 
   and also practicality of computational cost. 

*/
template<typename T> 
inline T jpoly ( const T a, const T b, 
                 const unsigned int n,
                 const T x )
{
  T poch = pochhammer<T>(a+1, n) / tgamma(n+1);
  T twoF1 = hypergeometric<T>(-1.0*n, n+a+b+1, a+1, (1.0-x)/2.0);
  return poch * twoF1;   
} 



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
template<typename T> 
inline void jpoly  ( const T* x, unsigned int Nx, unsigned int Np, 
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
    T a, b, c, d;
    const T *X = 0, *Y = 0, *Z = 0;
    T *H = 0, *V = 0;
    const unsigned int dim;
    unsigned int nthreads;
    bool weighted = false;
    void makeWeighted() { this->weighted = true; }


    jPoly ( unsigned int _Nx, unsigned int _n,
            T _a, T _b, T _c, T _d, 
            unsigned int _nthreads  )
      : Nx(_Nx), n(_n), 
        a(_a), b(_b), c(_c), d(_d), dim(3),
        nthreads(_nthreads) { init(); }

    jPoly ( unsigned int _Nx, unsigned int _n,
            T _a, T _b, T _c, unsigned int _nthreads  )
      : Nx(_Nx), n(_n), 
        a(_a), b(_b), c(_c), dim(2),
        nthreads(_nthreads) { init(); }

    jPoly ( unsigned int _Nx, unsigned int _n,
            T _a, T _b, unsigned int _nthreads  )
      : Nx(_Nx), n(_n),
        a(_a), b(_b), dim(1),
        nthreads(_nthreads) { init(); }

    jPoly ( const T* _X, const T* _Y, 
            unsigned int _Nx, 
            unsigned int _n,
            T _a, T _b, T _c, 
            unsigned int _nthreads  )
      : X(_X), Y(_Y), Nx(_Nx), n(_n), 
        a(_a), b(_b), c(_c), dim(2), 
        nthreads(_nthreads) 
    {
      init();
      computeV(this->X, this->Y);
    }

    jPoly ( const T* _X, 
            unsigned int _Nx,
            unsigned int _n,
            T _a, T _b,
            unsigned int _nthreads  )
      : X(_X), Nx(_Nx), n(_n), 
        a(_a), b(_b), dim(1),
        nthreads(_nthreads)
    {
      init();
      computeV(this->X);
    }

    void init  ()
    {
      if (this->dim == 2)
      {
        this->Np = static_cast<unsigned int>(0.5 * (n + 1) * (n + 2));
        this->V = (T*) calloc(Nx*Np, sizeof(T));
        this->H = (T*) calloc((n+1)*(n+1), sizeof(T));
        sFactors(n+1, a, b, c, this->H);
        this->ydx  = (T*) calloc(Nx, sizeof(T));
        this->x2m1 = (T*) calloc(Nx, sizeof(T));
        this->mxp1 = (T*) calloc(Nx, sizeof(T));
        this->Pk = (T*) calloc(Nx*(n+1), sizeof(T));
        this->Pnmk = (T*) calloc(Nx*(n+1)*(n+1), sizeof(T));
      }
      else if (this->dim == 3)
      {
        this->Np = dimPI3(n);
        this->V = (T*) calloc(Nx*Np, sizeof(T));
      }
      else if (this->dim == 1)
      {
        this->Np = n;
        this->V = (T*) calloc(Nx*Np, sizeof(T));
      }
    }

 
    void computeV(const T* _X, const T* _Y = 0, const T* _Z = 0)
    {
      if (this->dim == 3)
      {
        if (_Z == 0) { std::cout << "Z values empty! Exiting ..\n"; exit(1); }
        this->X =_X; this->Y = _Y; this->Z = _Z;
        #pragma omp parallel num_threads(nthreads)
        { 
          unsigned int blockind = 0;
          unsigned int blockcol = 0;
          for (unsigned int nn = 0; nn <= n; ++nn)
          {
            blockcol = 0;
            for (unsigned int kk = 0; kk <= nn; ++kk)
            {
              for (unsigned int jj = 0; jj <= kk; ++jj)
              {
                #pragma omp for 
                for (unsigned int i = 0; i < Nx; ++i)
                {
                  V[i + Nx*(blockind+blockcol)] = 
                    ( 1.0 / sFactors(nn,kk,jj,a,b,c,d) ) * 
                    jpoly<T>(2.0*kk+b+c+d+0.5, a-0.5, nn-kk, 2*X[i]-1) *
                    jpoly<T>(2.0*jj+c+d, b-0.5, kk-jj, 2*Y[i]/(1-X[i]) - 1) *
                    jpoly<T>(d-0.5, c-0.5, jj, 2*Z[i]/(1-X[i]-Y[i])-1) *
                    std::pow(1-X[i], kk-jj) * 
                    std::pow(1-X[i]-Y[i], jj);
                }
                blockcol += 1;
              }

            }
            blockind = blockind + rn3(nn);
          }
        } 

      }
      

      else if (this->dim == 2)
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

        // get bndry pts
        std::vector<int> x1inds, x0inds, y0inds, y1mxinds;
        for (unsigned int i = 0; i < Nx; ++i)
        {
          if (std::abs(X[i]-1) < 1e-15) {x1inds.push_back(i);}
          if (std::abs(X[i]) < 1e-15) {x0inds.push_back(i);}
          if (std::abs(Y[i]) < 1e-15) {y0inds.push_back(i);}
          if (std::abs(Y[i]-(1-X[i])) < 1e-15) {y1mxinds.push_back(i);}
        }

        int nx1 = x1inds.size();
        int nx0 = x0inds.size();
        int ny0 = y0inds.size();
        int ny1mx = y1mxinds.size();
        //std::cout << nx1 << " " << nx0 << " " << ny0 << " " << ny1mx << std::endl;
        if (not this->weighted)
        { 
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
              // substitute limit relation for X=1
              if (nx1 > 0)
              {
                for (unsigned int i = 0; i < nx1; ++i)
                {
                  if (kk != 0)
                  {
                    v[x1inds[i]] =  
                      ( 1.0 / H[kk + (n+1)*nn] * tgamma(nn+kk+b+c+1) / 
                      ( tgamma(nn-kk+1) * tgamma(nn+kk+b+c-(nn-kk)+1) ) ) *
                      ( pochhammer<T>(kk+c+b, kk) / (pow(2.0,kk) * tgamma(kk+1)) ) *
                      pow(2.0*Y[x1inds[i]], kk);
                  }
                  else
                  {
                    v[x1inds[i]] = 1;
                  }
                }
              }
              //if (nx0 > 0)
              //{
              //  #pragma omp simd
              //  for (unsigned int i = 0; i < nx0; ++i)
              //  {
              //    v[x0inds[i]] = 
              //      1.0 / H[kk + (n+1)*nn] *
              //      ( std::pow(-1.0, nn-kk) * tgamma(nn-kk+a-0.5+1) / 
              //      ( tgamma(nn-kk+1) * tgamma(nn-kk+a-0.5-(nn-kk)+1) ) ) *
              //      Pk[x0inds[i] + Nx*kk];
              //  }
              //}
              //if (ny0 > 0)
              //{
              //  #pragma omp simd
              //  for (unsigned int i = 0; i < ny0; ++i)
              //  {
              //    v[y0inds[i]] = 
              //      1.0 / H[kk + (n+1)*nn] *
              //      pnmk[y0inds[i] + Nx*(nn-kk)] * std::pow(mxp1[y0inds[i]], kk) *
              //      std::pow(-1.0, kk) * tgamma(kk+b-0.5+1) /
              //      ( tgamma(kk+1) * tgamma(kk+b-0.5-(kk)+1));
              //  }
              //}
              //if (ny1mx > 0)
              //{
              //  #pragma omp simd
              //  for (unsigned int i = 0; i < ny1mx; ++i)
              //  {
              //    v[y1mxinds[i]] = 
              //      1.0 / H[kk + (n+1)*nn] *
              //      pnmk[y1mxinds[i] + Nx*(nn-kk)] * std::pow(mxp1[y1mxinds[i]], kk) * 
              //      tgamma(kk+c-0.5+1) / 
              //      ( tgamma(kk+1) * tgamma(kk+c-0.5-(kk)+1));
              //  }
              //}
            }
            ind = ind + nn + 1;
          }
        }
        else if (this->weighted)
        {
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
                       std::pow(mxp1[i],kk) * Pk[i + Nx*kk] *
                       std::pow(X[i], a-0.5) * 
                       std::pow(Y[i], b-0.5) *
                       std::pow(1-X[i]-Y[i], c-0.5);
              }

            }
            ind = ind + nn + 1;
          }
        }
      }
      else if (this->dim == 1)
      {
        this->X = _X;
        jpoly<T>(X, Nx, Np, a, b, V);
      }
    }

    void computeCoeffs(T* F, T* X, T* Y, T* W, T* cF, T* _V = nullptr)
    {
      if (this->dim == 2)
      {
        this->computeV(X,Y);
        #pragma omp simd
        for (unsigned int i = 0; i < Nx; ++i)
        {
          F[i] *= W[i];
        }
        if (std::is_same_v<T, double>)
        {
          cblas_dgemv ( CblasColMajor, CblasTrans,
                        Nx, Np, 1.0, (double*) this->V, Nx, (double*) F, 1, 0.0, (double*) cF, 1);
          if (_V) { cblas_dcopy(Nx*Np, (double*) this->V, 1, (double*) _V, 1); }
        }
        else if (std::is_same_v<T, float>)
        {
          cblas_sgemv ( CblasColMajor, CblasTrans,
                        Nx, Np, 1.0, (float*) this->V, Nx, (float*) F, 1, 0.0, (float*) cF, 1);
          if (_V) { cblas_scopy(Nx*Np, (float*) this->V, 1, (float*) _V, 1); }
        }

      }
      else
      {
        std::cerr << "ERROR: Coefficient expansion is supported only for dim=2\n"; 
        exit(1);
      }
    }
    
    // compute only the first (m+1)*(m+2)/2 coeffs with m <= n
    void computeCoeffsM(T* F, T* X, T* Y, T* W, T* cF, unsigned int m, bool hasweights)
    {
      if (m > n)
      {
        std::cerr << "m <= n reaquired for truncated coeff computation\n";
        exit(1);
      } 
      
      if (this->dim == 2)
      {
        unsigned int Mp = static_cast<unsigned int>(0.5 * (m + 1) * (m + 2));
        this->computeV(X,Y);
        if (not hasweights)
        {
          #pragma omp simd
          for (unsigned int i = 0; i < Nx; ++i)
          {
            F[i] *= W[i];
          }
        }
        if (std::is_same_v<T, double>)
        {
          cblas_dgemv ( CblasColMajor, CblasTrans,
                        Nx, Mp, 1.0, (double*) this->V, Nx, (double*) F, 1, 0.0, (double*) cF, 1);
        }
        else if (std::is_same_v<T, float>)
        {
          cblas_sgemv ( CblasColMajor, CblasTrans,
                        Nx, Mp, 1.0, (float*) this->V, Nx, (float*) F, 1, 0.0, (float*) cF, 1);
        }

      }
      else
      {
        std::cerr << "ERROR: Coefficient expansion is supported only for dim=2\n"; 
        exit(1);
      }
    }

    ~jPoly  ()
    {
      if (H) { free(H); H = 0; }
      if (V) { free(V); V = 0; }
      if (ydx) { free(ydx); ydx = 0; }
      if (x2m1) { free(x2m1); x2m1 = 0; } 
      if (mxp1) { free(mxp1); mxp1 = 0; }
      if (Pk) { free(Pk); Pk = 0; }
      if (Pnmk) { free(Pnmk); Pnmk = 0; }
    }

  private: 
    T *ydx = 0, *x2m1 = 0, *mxp1 = 0;
    T *Pk = 0, *Pnmk = 0;

};

template double jpoly<double> ( const double a, const double b, 
                                const unsigned int n,
                                const double x );

template float jpoly<float> ( const float a, const float b, 
                              const unsigned int n,
                              const float x );


template void jpoly<double> ( const double* x, unsigned int Nx, 
                              unsigned int Np, const double a, 
                              const double b, double* V );

template void jpoly<float> (  const float* x, unsigned int Nx, 
                              unsigned int Np, const float a, 
                              const float b, float* V );

template class jPoly<double>;

template class jPoly<float>;


#endif
