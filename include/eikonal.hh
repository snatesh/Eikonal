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
  double *cu, *ul, *ub, *uh, *grad;
  double *cucopy; 
  double *clstorph, *cbstorph, *chstorph; 
  double *clstormh, *cbstormh, *chstormh; 
 
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
    // copy for gradient calculations
    cucopy = (double*) calloc(N, sizeof(double));
    // storage for constraint gradient calculations
    clstorph = (double*) calloc(nl, sizeof(double));
    cbstorph = (double*) calloc(nb, sizeof(double));
    chstorph = (double*) calloc(nh, sizeof(double));
    clstormh = (double*) calloc(nl, sizeof(double));
    cbstormh = (double*) calloc(nb, sizeof(double));
    chstormh = (double*) calloc(nh, sizeof(double));
  }

  ~optData ()
  {
    free(vecx); free(vecy);
    free(storb); free(storh);
    free(storl); free(cucopy);
    free(clstorph); free(clstormh);
    free(cbstorph); free(cbstormh);  
    free(chstorph); free(chstormh);
  }

};


static unsigned long count = 0;
double Fhlp  ( const double* cu, void* _data )
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

  double dxsq   = cblas_ddot(N, vecx, 1, vecx, 1);
  double dysq   = cblas_ddot(N, vecy, 1, vecy, 1);
  double rhssq  = cblas_ddot(N, rhs, 1, rhs, 1);
  double Fval   = dxsq + dysq - rhssq;
  
  if ( !(count % 100) ) 
  { 
    std::cout << "Eval #" << count 
              << " : F = " << Fval
              << std::endl; 
  }
  
  return Fval * Fval;

}


                

double F  ( unsigned int n, const double* cu, 
            double* grad, void* _data  )
{

  if (grad)
  {
    optData* data = (optData*) _data;
    double h = 1e-7;
    double fph, fmh, cuog;
    double* cucopy = data->cucopy;
    for (unsigned int i = 0; i < n; ++i) { cucopy[i] = cu[i]; }
    // compute finite difference approx to grad
    for (unsigned int i = 0; i < n; ++i)
    {
      // save current cu
      cuog = cucopy[i];
      // evaluate h in front
      cucopy[i] += h;
      fph = Fhlp(cucopy, _data);
      // evaluate h behind
      cucopy[i] -= 2*h;
      fmh = Fhlp(cucopy, _data);
      // reset cu
      cucopy[i] = cuog;
      // compute dFk/dx_j
      grad[i] = (fph - fmh) / (2.0 * h); 
    }
    
  }

  return Fhlp (cu, _data);

}

void cnstrntHlp ( unsigned int m, 
                  double* result, 
                  unsigned int n,
                  const double* cu,
                  void* f_data, 
                  unsigned int edge )
{
  optData* data = (optData*) f_data;
  if (edge == 0)
  {
    double* vl    = data->vl;
    double* ul    = data->ul;
    double* storl = data->storl;

    memset(storl, 0, m*sizeof(double));
    cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                  m, n, 1.0, vl, n, cu, 1, 1.0, storl, 1 );
    
    // u2 = q-ul, q=0 => c = u2_l +ul
    for (unsigned int i = 0; i < m; ++i)
    {
      result[i] = storl[i] + ul[i] ;
    } 
  }
  else if (edge == 1)
  {
    double* vb    = data->vb;
    double* ub    = data->ub;
    double* storb = data->storb;

    memset(storb, 0, m*sizeof(double));
    cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                  m, n, 1.0, vb , n, cu, 1, 0.0, storb, 1 );
    for (unsigned int i = 0; i < m; ++i)
    {
      result[i] = storb[i] + ub[i];
    } 
  }
  else if (edge == 2)
  {
    double* vh    = data->vh;
    double* uh    = data->uh;
    double* storh = data->storh;
    memset(storh, 0, m*sizeof(double));
    cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                  m, n, 1.0, vh , n, cu, 1, 0.0, storh, 1 );
    for (unsigned int i = 0; i < m; ++i)
    {
      result[i] = storh[i] + uh[i];
    } 
  }
}

void cnstrntGhlp (  unsigned int m, 
                    double* result, 
                    unsigned int n,
                    const double* cu,
                    double* grad, 
                    void* f_data,  
                    unsigned int edge )
{

  optData* data = (optData*) f_data;
  double* cucopy    = data->cucopy;
  for (unsigned int j = 0; j < n; ++j) { cucopy[j] = cu[j]; }
  double cuog, h = 1e-8;
  double *cstorph, *cstormh;
  if (edge == 0)
  {
    cstorph  = data->clstorph;
    cstormh  = data->clstormh;
  }
  else if (edge == 1)
  {
    cstorph  = data->cbstorph;
    cstormh  = data->cbstormh;
  }
  
  else if (edge == 2)
  {
    cstorph  = data->chstorph;
    cstormh  = data->chstormh;
  }

  // compute dc_i/dx_j
  for (unsigned int j = 0; j < n; ++j)
  {
    cuog = cucopy[j];
    // evaluate h in front
    cucopy[j] += h;
    cnstrntHlp(m, cstorph, n, cucopy, f_data, edge);
    // evaluate h behind
    cucopy[j] -= 2*h;
    cnstrntHlp(m, cstormh, n, cucopy, f_data, edge);
    // reset cu
    cucopy[j] = cuog;
    for (unsigned int i = 0; i < m; ++i)
    {
      grad[j + n*i] = (cstorph[i] - cstormh[i]) / (2.0 * h);
    }
  }

}


void cl ( unsigned int m, double* result, unsigned int n, 
          const double* cu, double* grad, void* f_data)
{
  if (grad)
  {
    cnstrntGhlp(m, result, n, cu, grad, f_data, 0);
  }
  cnstrntHlp(m, result, n, cu, f_data, 0);
}

void cb ( unsigned int m, double* result, unsigned int n, 
          const double* cu, double* grad, void* f_data)
{
  if (grad)
  {
    cnstrntGhlp(m, result, n, cu, grad, f_data, 1);
  } 
  cnstrntHlp(m, result, n, cu, f_data, 1);
}

void ch ( unsigned int m, double* result, unsigned int n, 
          const double* cu, double* grad, void* f_data)
{
  if (grad)
  {
    cnstrntGhlp(m, result, n, cu, grad, f_data, 2);
  } 
  cnstrntHlp(m, result, n, cu, f_data, 2);
}


#endif
