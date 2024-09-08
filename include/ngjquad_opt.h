#ifndef _NGJQUAD_OPT_H
#define _NGJQUAD_OPT_H
#include<string>
#include<iostream>
#include<iomanip>
#include<vector>
#include<complex.h>
#include<cblas.h>
#include<lapacke.h>
#include<nlopt.h>
#include<jacobi_mat_ON_tri.h>
#include<jPoly.h>

typedef double _Complex Complex;

typedef struct
{
  double *Vm, *Hm, *Fk;
  double a,b,c;
  unsigned int N;
  unsigned int m;
  unsigned int M = static_cast<unsigned int>(0.5 * m * (m + 1));
} nlopt_func_data;

unsigned long count = 0;

double nloptF(unsigned int n, const double* Zk, double* grad, void* data)
{
  ++count;
  nlopt_func_data* d = (nlopt_func_data*) data;
  unsigned int m = d->m;
  unsigned int M = d->M;
  unsigned int N = static_cast<unsigned int>(n / 3.0);
  double a, b, c; a = d->a; b = d->b; c = d->c;
  const double *Xk, *Yk, *Wk;
  double *Hm, *Vm, *Fk;
  Xk = Zk; Yk = Zk + N; Wk = Zk + 2*N;
  Hm = d->Hm; Vm = d->Vm; Fk = d->Fk;
  double* I0 = (double*) calloc(M, sizeof(double)); I0[0] = 1;  
  jPoly_tri<double>(Xk, Yk, Hm, N, m-1, a, b, c, Vm);
  cblas_dgemv(CblasColMajor,  CblasTrans,  N, M, 1.0, Vm, N, Wk, 1, -1.0, I0, 1);
  for (unsigned int i =0; i < M; ++i) { Fk[i] = I0[i]; }
  free(I0);
  double norm2 = pow(cblas_dnrm2(M, Fk, 1),2);
  if ( !(count % 1000) ) { std::cout << "Eval #" << count << " : F = " << norm2 << std::endl; }
  return norm2;
}

  
void nloptieqC( unsigned int m, double *result, unsigned int n, const double* Zk, 
                double* grad=nullptr, void* f_data=nullptr )
{
  unsigned int N = static_cast<unsigned int>(n / 3.0);
  for (unsigned int i = 0; i < m-1; ++i)
  {
    result[i] = Zk[i] + Zk[i + N] - 1.0; 
  }
}

double nlopteqC(  unsigned int n, const double* Zk, 
                  double* grad=nullptr, void* data=nullptr ) 
{
  double sumw = 0.0;
  unsigned int N = static_cast<unsigned int>(n / 3.0);
  for (unsigned int i = 2*N; i < n; ++i) { sumw += Zk[i]; }
  return sumw - 1.0;
} 


void F( double* Zk, double* Vm, unsigned int N, unsigned int m, 
        double* Hm, double a, double b, double c, double* Fk )
{
  const double *Xk, *Yk, *Wk;
  Xk = Zk; Yk = Zk + N; Wk = Zk + 2*N;
  unsigned int M = static_cast<unsigned int>(0.5 * m * (m + 1));
  double* I0 = (double*) calloc(M, sizeof(double)); I0[0] = 1;  
  jPoly_tri<double>(Xk, Yk, Hm, N, m-1, a, b, c, Vm);
  cblas_dgemv(CblasColMajor,  CblasTrans,  N, M, 1.0, Vm, N, Wk, 1, -1.0, I0, 1);
  for (unsigned int i =0; i < M; ++i) { Fk[i] = I0[i]; }
  free(I0);
}



void newton(double* Fk, double* Vm, double* Hm, unsigned int N, unsigned int m, double a, double b, double c, double* Zk, bool wolfe)
{
  double h = 1e-7; 
  double tol = 1e-7; 
  double tol_up = 1e5; 
  unsigned int maxiter = 1000;
  unsigned int M = static_cast<unsigned int>(0.5 * m * (m + 1));
  double pk = cblas_dnrm2(M, Fk, 1);
  double* Fkph    = (double*) calloc(M, sizeof(double));
  double* Fkmh    = (double*) calloc(M, sizeof(double));
  double* Zkph    = (double*) calloc(3*N, sizeof(double));
  double* Zkmh    = (double*) calloc(3*N, sizeof(double));
  double* gradFk  = (double*) calloc(M*3*N, sizeof(double)); 
  double* dZk     = (double*) calloc(3*N, sizeof(double)); 
  double* S       = (double*) calloc(M, sizeof(double));

  double rho = 0.9; double gam = 1e-4;
  unsigned int iter = 0; lapack_int rank[1];
  while (pk > tol && pk < tol_up && iter < maxiter)
  {
    iter += 1;
    /* dZk initialized to Fk, overwritten by dgelsd to dZk 
       which is lsq sol to gradFk*dZk = Fk */
    for (unsigned int i = 0; i < M; ++i)    { dZk[i] = Fk[i]; }
    for (unsigned int i = M; i < 3*N; ++i)  { dZk[i] = 0; }
    for (unsigned int i = 0; i < 3*N; ++i) { Zkph[i] = Zkmh[i] = Zk[i]; }
    // compute finite difference approx to grad
    for (unsigned int jj = 0; jj < 3*N; ++jj)
    {
      // eval above and below Zk
      Zkph[jj] += h; Zkmh[jj] -= h;
      F(Zkph, Vm, N, m, Hm, a, b, c, Fkph);
      F(Zkmh, Vm, N, m, Hm, a, b, c, Fkmh);
      // compute dFk/dx_j
      for (unsigned int ii = 0; ii < M; ++ii)
      {
        gradFk[ii + M*jj] = (Fkph[ii] - Fkmh[ii]) / (2.0 * h);
      }
      // revert to original Zk
      Zkph[jj] -= h;
      Zkmh[jj] += h;
    }
    // compute descent direction
    if (LAPACKE_dgelsd(LAPACK_COL_MAJOR, M, 3*N, 1, gradFk, M, dZk, 3*N, S, 1e-16, rank))
    {
      std::cerr << "ERROR: Lapack dgelsd: Pseudoinverse" << std::endl;
    }
    // linesearch with Wolfe conditions
    double alph = 0.0005; double wfnrm;
    for (unsigned int i = 0; i < 3*N; ++i) { Zk[i] -= alph*dZk[i]; }
    F(Zk, Vm, N, m, Hm, a, b, c, Fk);
    pk = cblas_dnrm2(M, Fk, 1);

    if (wolfe)
    {
      for (unsigned int iter_i = 0; iter_i < maxiter; ++iter_i)
      {
        cblas_dgemv(CblasColMajor,  CblasNoTrans, M, 3*N, gam*alph, gradFk, M, dZk, 1, 1.0, Fk, 1);
        wfnrm = cblas_dnrm2(M, Fk, 1);
        std::cout << wfnrm << std::endl;
        if (pk > wfnrm)
        {
          alph = rho*alph;
          for (unsigned int i = 0; i < 3*N; ++i) { Zk[i] -= alph*dZk[i]; }
          F(Zk, Vm, N, m, Hm, a, b, c, Fk);
          pk = cblas_dnrm2(M, Fk, 1);
        } 
      }
    }

    if ( !(iter%100) ) { std::cout << "objective norm : " << pk << std::endl; }
  }
  free(Fkph);
  free(Fkmh);
  free(Zkph);
  free(Zkmh);
  free(gradFk);
  free(dZk);
  free(S);
}

#endif
