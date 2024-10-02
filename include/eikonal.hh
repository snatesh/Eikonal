#ifndef _EIKONAL_H
#define _EIKONAL_H

#include<cblas.h>
#include<nlopt.h>
#include<cstdlib>
#include<iostream>
#include<cstring>
struct optData
{

  unsigned int N;
  double *vecx, *vecy, *vecxx, *vecyy;
  double *vecf, *vecrhs;
  double *storl, *storb, *storh;
  double *Dx, *Dy, *rhs;
  double *vl, *vb, *vh;
  unsigned int nl, nb, nh;
  double *cu;
 
  
 
  optData ( unsigned int _N,
            unsigned int _nl,
            unsigned int _nb,
            unsigned int _nh,
            double* _Dx, double* _Dy,
            double* _rhs,
            double* _vl, double* _vb, 
            double* _vh, double* _cu ) 
    : N(_N), Dx(_Dx), Dy(_Dy), rhs(_rhs), 
      vl(_vl), vb(_vb), vh(_vh),
      nl(_nl), nb(_nb), nh(_nh),
      cu(_cu)
  {
    // Coeff derivs in x,y
    vecx = (double*) calloc(N, sizeof(double));
    vecy = (double*) calloc(N, sizeof(double));
    // outer product of coef derivs
    vecxx = (double*) calloc(N*N, sizeof(double));
    vecyy = (double*) calloc(N*N, sizeof(double));
    // objective function 
    // (combination of outer products, vectorized)
    vecf = (double*) calloc(N*N, sizeof(double));
    // outer product of coeffs of (1/f)
    // first is entry is 1 the case of f \equiv 1
    vecrhs = (double*) calloc(N*N, sizeof(double)); 
    cblas_dger  ( CblasColMajor, N, N, 
                  1.0, rhs, 1, rhs, 1, vecrhs, N );
    // storage for contraint result
    storl = (double*) calloc(nl, sizeof(double));
    storb = (double*) calloc(nb, sizeof(double));
    storh = (double*) calloc(nh, sizeof(double));
  }

  ~optData ()
  {
    free(vecx); free(vecy);
    free(vecxx); free(vecyy);
    free(vecf); free(vecrhs);
    free(storb); free(storh);
    free(storl);
  }

};
                
static unsigned long count = 0;
double F  ( unsigned int n, const double* cu, 
            double* grad, void* _data  )
{
  ++count;
  optData* data     = (optData*) _data;
  unsigned int N    = data->N;
  unsigned int NN   = N*N;
  double* vecx      = data->vecx; 
  double* vecy      = data->vecy;
  double* vecxx     = data->vecxx; 
  double* vecyy     = data->vecyy; 
  double* Dx        = data->Dx;
  double* Dy        = data->Dy;
  double* vecf      = data->vecf;
  double* vecrhs    = data->vecrhs;
 
  memset(vecx, 0, N*sizeof(double));
  memset(vecy, 0, N*sizeof(double));
  memset(vecf, 0, NN*sizeof(double)); 
  memset(vecxx, 0, NN*sizeof(double));
  memset(vecyy, 0, NN*sizeof(double));

  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                N, N, 1.0, Dx, N, cu, 1, 0.0, vecx, 1 ); 
  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                N, N, 1.0, Dy, N, cu, 1, 0.0, vecy, 1 ); 
  cblas_dger  ( CblasColMajor, N, N, 
                1.0, vecx, 1, vecx, 1, vecxx, N );
  cblas_dger  ( CblasColMajor, N, N, 
                1.0, vecy, 1, vecy, 1, vecyy, N );
  #pragma omp simd
  for (unsigned int i = 0; i < NN; ++i)
  {
    vecf[i] = vecxx[i] + vecyy[i] - vecrhs[i];
  } 

  // frobenius norm of A \equiv 2 norm vec(A)  
  double nrm  = cblas_dnrm2 (NN, vecf, 1);
  double nrm2 = nrm * nrm;
  
  if ( !(count % 100000) ) 
  { 
    std::cout << "Eval #" << count 
              << " : F = " << nrm2 
              << std::endl; 
  }
  
  return nrm2;
}



void cl ( unsigned int m, double* result, unsigned int n, 
          const double* cu, double* grad, void* f_data)
{
  optData* data = (optData*) f_data;
  double* vl = data->vl;
  double* storl = data->storl;
  memset(storl, 0, m*sizeof(double));
  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                m, n, 1.0, vl, n, cu, 1, 1.0, storl, 1 );

  for (unsigned int i = 0; i < m; ++i)
  {
    result[i] = storl[i];
  } 
}

void cb ( unsigned int m, double* result, unsigned int n, 
          const double* cu, double* grad, void* f_data)
{ 

  optData* data = (optData*) f_data;
  double* vb = data->vb;
  double* storb = data->storb;
  memset(storb, 0, m*sizeof(double));
  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                m, n, 1.0, vb , n, cu, 1, 0.0, storb, 1 );
  for (unsigned int i = 0; i < m; ++i)
  {
    result[i] = storb[i];
  } 
}

void ch ( unsigned int m, double* result, unsigned int n, 
          const double* cu, double* grad, void* f_data)
{
  optData* data = (optData*) f_data;
  double* vh = data->vh;
  double* storh = data->storh;
  memset(storh, 0, m*sizeof(double));
  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                m, n, 1.0, vh , n, cu, 1, 0.0, storh, 1 );
  for (unsigned int i = 0; i < m; ++i)
  {
    result[i] = storh[i];
  } 
}


#endif
