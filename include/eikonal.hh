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
  double *vecx, *vecy;
  double *storl, *storb, *storh;
  double *Dx, *Dy, *rhs;
  double *vl, *vb, *vh;
  unsigned int nl, nb, nh;
  double *cu, *ul, *ub, *uh;
 
  
 
  optData ( unsigned int _N,
            double* _Dx, double* _Dy,
            double* _rhs, double* _cu,
            unsigned int _nl = 0,
            unsigned int _nb = 0,
            unsigned int _nh = 0,
            double* _vl = nullptr, 
            double* _vb = nullptr, 
            double* _vh = nullptr, 
            double* _ul = nullptr, 
            double* _ub = nullptr,
            double* _uh = nullptr ) 
    : N(_N), Dx(_Dx), Dy(_Dy), rhs(_rhs), 
      vl(_vl), vb(_vb), vh(_vh),
      nl(_nl), nb(_nb), nh(_nh),
      cu(_cu), ul(_ul), ub(_ub), uh(_uh)
  {
    // Coeff derivs in x,y
    vecx = (double*) calloc(N, sizeof(double));
    vecy = (double*) calloc(N, sizeof(double));
    // storage for contraint result
    storl = (double*) calloc(nl, sizeof(double));
    storb = (double*) calloc(nb, sizeof(double));
    storh = (double*) calloc(nh, sizeof(double));
  }

  ~optData ()
  {
    free(vecx); free(vecy);
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
  double* vecx      = data->vecx; 
  double* vecy      = data->vecy;
  double* Dx        = data->Dx;
  double* Dy        = data->Dy;
  double* rhs       = data->rhs;
 
  memset(vecx, 0, N*sizeof(double));
  memset(vecy, 0, N*sizeof(double));


  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                N, N, 1.0, Dx, N, cu, 1, 0.0, vecx, 1 ); 
  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                N, N, 1.0, Dy, N, cu, 1, 0.0, vecy, 1 ); 

  double nrmx = cblas_dnrm2(N, vecx, 1);
  double nrmy = cblas_dnrm2(N, vecy, 1);
  double nrmrhs = cblas_dnrm2(N, rhs, 1);
  double nrm2 = nrmx * nrmx + nrmy * nrmy - nrmrhs * nrmrhs;
  
  if ( !(count % 100) ) 
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
  double* ul = data->ul;
  double* storl = data->storl;
  memset(storl, 0, m*sizeof(double));
  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                m, n, 1.0, vl, n, cu, 1, 1.0, storl, 1 );
  
  // u2 = q-ul, q=0 => c = u2_l +ul
  for (unsigned int i = 0; i < m; ++i)
  {
    result[i] = storl[i] + ul[i];
  } 
}

void cb ( unsigned int m, double* result, unsigned int n, 
          const double* cu, double* grad, void* f_data)
{ 

  optData* data = (optData*) f_data;
  double* vb = data->vb;
  double* ub = data->ub;
  double* storb = data->storb;
  memset(storb, 0, m*sizeof(double));
  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                m, n, 1.0, vb , n, cu, 1, 0.0, storb, 1 );
  for (unsigned int i = 0; i < m; ++i)
  {
    result[i] = storb[i] + ub[i];
  } 
}

void ch ( unsigned int m, double* result, unsigned int n, 
          const double* cu, double* grad, void* f_data)
{
  optData* data = (optData*) f_data;
  double* vh = data->vh;
  double* uh = data->uh;
  double* storh = data->storh;
  memset(storh, 0, m*sizeof(double));
  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                m, n, 1.0, vh , n, cu, 1, 0.0, storh, 1 );
  for (unsigned int i = 0; i < m; ++i)
  {
    result[i] = storh[i] + uh[i];
  } 
}


#endif
