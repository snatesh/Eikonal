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


typedef struct nlopt_func_data
{
  double *Vm, *Hm, *Fk, *Z0, *Jn1, *Jn2, *X0, *Y0;
  double *Vm_T, *W0, *S, *F0;
  Complex *Jn, *XY0;
  double a,b,c;
  unsigned int N, M, n, m;
  bool run0; 
  double minf, minf1;
 
  nlopt_func_data ( unsigned int _m, unsigned int _n, 
                    double _a, double _b, double _c )
    : n(_n), m(_m), a(_a), b(_b), c(_c), run0(true)
  {
    this->N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
    this->M = static_cast<unsigned int>(0.5 * m * (m + 1));
    // generate jacobi matrices for x,y 
    this->Hm = (double*) calloc(m*m, sizeof(double));
    this->Jn1 = (double*) calloc(N*N, sizeof(double));
    this->Jn2 = (double*) calloc(N*N, sizeof(double));
    this->Jn = (Complex*) calloc(N*N, sizeof(Complex));
    this->XY0 = (Complex*) calloc(N, sizeof(Complex));
    this->X0  = (double*) calloc(N, sizeof(double));
    this->Y0  = (double*) calloc(N, sizeof(double));
    this->Vm = (double*) calloc(N*M, sizeof(double));
    this->Vm_T = (double*) calloc(M*N, sizeof(double));
    this->W0 = (double*) calloc(M, sizeof(double)); W0[0] = 1.0;
    this->S = (double*) calloc(M, sizeof(double)); 
    this->Z0 = (double*) calloc(3*N, sizeof(double));
    this->F0 = (double*) calloc(3*N, sizeof(double)); 
    this->Fk = this->F0;
  }

  ~nlopt_func_data()
  {
    free(Jn1); free(Jn2); free(Jn);
    free(XY0); free(X0); free(Y0);
    free(Hm); free(Vm); free(Vm_T);
    free(W0); free(F0); free(Z0); free(S); 
    Jn1 = Jn2 = X0 = Y0 = 0; Jn = XY0 = 0;
    Hm = Vm = Vm_T = W0 = F0 = Z0 = S = 0; 
  }
};

inline void set_args  ( int argc, char* argv[], 
                        nlopt_algorithm& alg, int& n, 
                        int& m, double& tol, double& tolc, 
                        bool& use_newton, bool& use_wolfe, 
                        double& alph  )
{
  std::cout << "SEARCH FOR NEAR OPTIMAL GAUSS QUAD" << std::endl;
  switch(argc)
  {
    case 9:
    {
      unsigned int algtype = std::stoi(argv[1]);
      switch(algtype)
      {
        case 0:
        { 
          alg = NLOPT_LN_NELDERMEAD;
          std::cout << "ALGORITHM: NELDERMEAD\n";
          break;
        }
        case 1:
        {
          alg = NLOPT_LN_COBYLA;
          std::cout << "ALGORITHM: COBYLA\n";
          break;
        }
        case 2:
        {  
          alg = NLOPT_LN_SBPLX;
          std::cout << "ALGORITHM: SBPLX";
          break;
        }
        default:
        {
          std::cerr << "Algtype not supported" << std::endl;
          break;
        }
      }
      n = std::stoi(argv[2]);
      m = std::stoi(argv[3]);
      tol = std::stod(argv[4]); 
      tolc = std::stod(argv[5]);
      use_newton = (bool) std::stoi(argv[6]); 
      use_wolfe = (bool) std::stoi(argv[7]);
      alph = std::stod(argv[8]);
      break;
    }
    case 1:
    {
      alg = NLOPT_LN_NELDERMEAD;
      std::cout << "ALGORITHM: NELDERMEAD\n";
      n = 4; m = 6;
      tol = 1e-5;
      tolc = 1e-5;
      use_wolfe = false;
      use_newton = false;
      alph = 0;
      break;
    }
    default:
    {
      std::cerr << 
        "Usage: ./ngjquad_opt algtype, n m tol tolc use_newton use_wolfe alph\n";
      std::cerr << "  algtype (int) - 0-2\n";
      std::cerr << "  n,m (int) > 0\n";
      std::cerr << "  tol,tolc (double) > 0\n";
      std::cerr << "  use_newton, use_wolfe (int) - 0,1\n";
      std::cerr << "  alph (double) > 0\n";
      exit(1);
    }
  }
  std::cout << "\nPARAMETERS:" << std::endl;
  std::cout << std::setw(15) << "n          = " << n << std::endl;
  std::cout << std::setw(15) << "m          = " << m << std::endl;
  std::cout << std::setw(15) << "tol        = " << tol << std::endl;
  std::cout << std::setw(15) << "tolc       = " << tolc << std::endl;
  std::cout << std::setw(15) << "use_newton = " << use_newton << std::endl;
  std::cout << std::setw(15) << "use_wolfe  = " << use_wolfe << std::endl;
  std::cout << std::setw(15) << "alpha      = " << alph << std::endl;
}

unsigned long count = 0;
inline double nloptF(unsigned int n, const double* Zk, double* grad, void* data)
{
  ++count;
  nlopt_func_data* d = (nlopt_func_data*) data;
  unsigned int m = d->m;
  unsigned int M = d->M;
  unsigned int N = d->N;
  double a = d->a, b = d->b, c = d->c;
  const double *Xk, *Yk, *Wk;
  double *Hm, *Vm, *Fk;
  Xk = Zk; Yk = Zk + N; Wk = Zk + 2*N;
  Hm = d->Hm; Vm = d->Vm; Fk = d->Fk;
  double* I0 = (double*) calloc(M, sizeof(double)); I0[0] = 1;  
  jPoly_tri<double>(Xk, Yk, Hm, N, m-1, a, b, c, Vm);
  cblas_dgemv(CblasColMajor,  CblasTrans,  N, M, 1.0, Vm, N, Wk, 1, -1.0, I0, 1);
  #pragma omp parallel for
  for (unsigned int i =0; i < M; ++i) { Fk[i] = I0[i]; }
  free(I0);
  double norm2 = pow(cblas_dnrm2(M, Fk, 1),2);
  if ( !(count % 100000) && d->run0 ) 
  { 
    std::cout << "Eval #" << count << " : F = " << norm2 << std::endl; 
  }
  return norm2;
}

inline void cond(double* Vm, unsigned int N, unsigned int M, double* rcond, bool norm2=true)
{
  if (norm2)
  {
    unsigned int dimS = (N <= M) ? N : M;
    double* S = (double*) calloc(dimS, sizeof(double));
    double* superb = (double*) calloc(dimS, sizeof(double));
    LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'N','N', N, M, Vm, N, S, nullptr, 1, nullptr, 1, superb);
    rcond[0] = S[dimS-1] / S[0];
    free(S); free(superb);
  }
  else
  {
    //  inf norm of Vm
    double  normVm = LAPACKE_dlange(LAPACK_COL_MAJOR, '1', N, M, Vm, N);
    // LU of A
    int* ipiv = (int*) calloc(N, sizeof(int));
    LAPACKE_dgetrf (LAPACK_COL_MAJOR, N, M, Vm, N, ipiv); 
    free(ipiv);
    LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', N, Vm, M, normVm, rcond);
  }
}

unsigned long count1 = 0;
inline double nloptF1(unsigned int n, const double* Zk, double* grad, void* data)
{
  ++count1;
  nlopt_func_data* d = (nlopt_func_data*) data;
  const double *Xk = Zk; 
  const double *Yk = Zk + d->N;
  jPoly_tri<double>(Xk, Yk, d->Hm, d->N, d->m-1, d->a, d->b, d->c, d->Vm);
  double rcond[1];
  cond(d->Vm, d->N, d->M, rcond);
  if ( !(count1 % 100000) ) 
  { 
    std::cout << "Eval #" << count << " : F = " << 1.0 / rcond[0] << std::endl; 
  }
  return 1.0 / rcond[0];
}

inline void nloptieqC ( unsigned int m, double *result, unsigned int n, 
                        const double* Zk, double* grad=nullptr, 
                        void* f_data=nullptr )
{
  unsigned int N = static_cast<unsigned int>(n / 3.0);
  #pragma omp parallel for
  for (unsigned int i = 0; i < m-1; ++i)
  {
    result[i] = Zk[i] + Zk[i + N] - 1.0; 
  }
}

inline void nloptieqC1 (  unsigned int m, double *result, unsigned int n, 
                          const double* Zk, double* grad=nullptr, 
                          void* f_data=nullptr )
{
  unsigned int N = static_cast<unsigned int>(n / 2.0);
  #pragma omp parallel for
  for (unsigned int i = 0; i < m-1; ++i)
  {
    result[i] = Zk[i] + Zk[i + N] - 1.0; 
  }
}

inline double nlopteqC( unsigned int n, const double* Zk, 
                        double* grad=nullptr, void* data=nullptr ) 
{
  double sumw = 0.0;
  unsigned int N = static_cast<unsigned int>(n / 3.0);
  #pragma omp parallel for
  for (unsigned int i = 2*N; i < n; ++i) { sumw += Zk[i]; }
  return sumw - 1.0;
} 

inline double nlopteqC1 ( unsigned int n, const double* Zk,
                          double* grad, void* data) 
{
  nlopt_func_data* d = (nlopt_func_data*) data;
  d->run0 = false;
  return nloptF(d->n, Zk, nullptr, d) - d->minf;
  
}  

inline void F( double* Zk, double* Vm, unsigned int N, unsigned int m, 
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

inline void init_opt  ( nlopt_func_data* d )
{

  unsigned int m = d->m; unsigned int n = d->n; 
  unsigned int M = d->M; unsigned int N = d->N;
  double* Jn1 = d->Jn1; double* Jn2 = d->Jn2;
  Complex* Jn = d->Jn; Complex* XY0 = d->XY0;
  double* X0 = d->X0; double* Y0 = d->Y0;
  double* Hm = d->Hm; double* Vm = d->Vm;
  double* Vm_T = d->Vm_T; double* W0 = d->W0;
  double* S = d->S; double* Z0 = d->Z0; 
  double* F0 = d->F0; 
  double a = d->a, b = d->b, c = d->c;
  structure_factors_tri<double>(m, a, b, c, Hm);
  // generate jacobi matrices for x,y 
  jacobi_mat_ON_tri<double>(n, a, b, c, Jn1, Jn2);
  // compute initial nodes and weights from eigenvalues of Jn
  for (unsigned int i = 0; i < N*N; ++i) { Jn[i] = Jn1[i] + I*Jn2[i]; }
  if (LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'N', N, Jn, N, XY0, NULL, N, NULL, N))
  {
    std::cerr << "ERROR: Lapack ZGEEV: Eigenvalues" << std::endl;
  }
  for (unsigned int i = 0; i < N; ++i) 
  { 
    X0[i] = creal(XY0[i]); 
    Y0[i] = cimag(XY0[i]); 
  } 
  // evaluate Vandermonde on initial nodes
  jPoly_tri<double>(X0, Y0, Hm, N, m-1, a, b, c, Vm);
  unsigned int i = 0;
  for (unsigned int col = 0; col < M; ++col)
  {
    for (unsigned int row = 0; row < N; ++row)
    {
      Vm_T[col + row*M] = Vm[i];
      i += 1;
    }
  }
  lapack_int rank[1]; 
  // solve least squares system for initial weights
  if (LAPACKE_dgelsd(LAPACK_COL_MAJOR, M, N, 1, Vm_T, M, W0, M, S, -1.0, rank))
  {
    std::cerr << "ERROR: Lapack dgelsd: Pseudoinverse" << std::endl;
  }
  // copy into Z0 for opt routines
  for (unsigned int i = 0; i < N; ++i)
  {
    Z0[i]       = X0[i];
    Z0[i + N]   = Y0[i];
    Z0[i + 2*N] = W0[i];
  }
  F(Z0, Vm, N, m, Hm, a, b, c, F0);
}

inline void nlopt_run ( nlopt_algorithm alg, nlopt_func_data* d,
                        double tol, double tolc )
{
  std::cout <<"\n BEGIN NLOPT 0 \n"; 
  unsigned int N = d->N; 
  // x, y , w > 0
  double* lb = (double*) calloc(3*N, sizeof(double)); 
  // x, y , w < 1
  double* ub = (double*) calloc(3*N, sizeof(double));
  double* tolieq = (double*) calloc(2*N, sizeof(double));
  for (unsigned int i = 0; i < 3*N; ++i) { ub[i] = 1; }
  for (unsigned int i = 0; i < 2*N; ++i) { tolieq[i] = tolc; }

  nlopt_opt opt = nlopt_create(alg, 3*N); 
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);
  nlopt_set_min_objective(opt, nloptF, d);
  nlopt_add_equality_constraint(opt, nlopteqC, NULL, tolc);
  nlopt_add_inequality_mconstraint(opt, 2*N, nloptieqC, NULL, tolieq);
  nlopt_set_xtol_rel(opt, tol);
  nlopt_set_stopval(opt, tol);
  double minF;
  if (nlopt_optimize(opt, d->Z0, &minF) < 0) 
  {
    std::cerr << "NLOPT failed!" << std::endl;
  }
  else 
  {
    std::cout << "found minimum with objective val = " << minF << std::endl;
  }
  nlopt_destroy(opt);
  free(lb); free(ub);
  free(tolieq);
  d->minf = minF;
  std::cout <<"END NLOPT \n\n";
}

inline void nlopt_run1 ( nlopt_algorithm alg, nlopt_func_data* d,
                         double tol, double tolc )
{
  std::cout <<"\n BEGIN NLOPT 1 \n"; 
  unsigned int N = d->N; 
  // x, y > 0
  double* lb = (double*) calloc(2*N, sizeof(double)); 
  // x, y < 1
  double* ub = (double*) calloc(2*N, sizeof(double));
  double* tolieq = (double*) calloc(2*N, sizeof(double));
  for (unsigned int i = 0; i < 2*N; ++i) { ub[i] = 1; }
  for (unsigned int i = 0; i < 2*N; ++i) { tolieq[i] = tolc; }

  nlopt_opt opt = nlopt_create(alg, 2*N); 
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);
  nlopt_set_min_objective(opt, nloptF1, d);
  nlopt_add_equality_constraint(opt, nlopteqC1, NULL, tolc);
  nlopt_add_inequality_mconstraint(opt, 2*N, nloptieqC1, NULL, tolieq);
  nlopt_set_stopval(opt, tol);
  double minF;
  if (nlopt_optimize(opt, d->Z0, &minF) < 0) 
  {
    std::cerr << "NLOPT failed!" << std::endl;
  }
  else 
  {
    std::cout << "found minimum with objective val = " << minF << std::endl;
  }
  nlopt_destroy(opt);
  free(lb); free(ub);
  d->minf1 = minF;
  free(tolieq);
  std::cout <<"END NLOPT \n\n";
}

inline void newton( double* Fk, double* Vm, double* Hm, 
                    unsigned int N, unsigned int m, 
                    double a, double b, double c, double* Zk, 
                    bool use_wolfe, double alph = 0.01 )
{
  std::cout << "BEGIN NEWTON" << std::endl;
  double h = 1e-7; 
  double tol = 1e-7; 
  double tol_up = 1e5; 
  unsigned int maxiter = 100;
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
    double wfnrm;
    for (unsigned int i = 0; i < 3*N; ++i) { Zk[i] -= alph*dZk[i]; }
    F(Zk, Vm, N, m, Hm, a, b, c, Fk);
    pk = cblas_dnrm2(M, Fk, 1);
    if (use_wolfe)
    {
      for (unsigned int iter_i = 0; iter_i < maxiter; ++iter_i)
      {
        cblas_dgemv(CblasColMajor,  CblasNoTrans, M, 3*N, gam*alph, gradFk, M, dZk, 1, 1.0, Fk, 1);
        wfnrm = cblas_dnrm2(M, Fk, 1);
        if (pk > wfnrm)
        {
          alph = rho*alph;
          for (unsigned int i = 0; i < 3*N; ++i) { Zk[i] -= alph*dZk[i]; }
          F(Zk, Vm, N, m, Hm, a, b, c, Fk);
          pk = cblas_dnrm2(M, Fk, 1);
          if ( !(iter_i%10) ) { std::cout << "norm(F) (wolfe): " << pk << std::endl; }
        } 
      }
    }

    if ( !(iter%10) ) { std::cout << "norm(F) : " << pk << std::endl; }
  }

  free(Fkph);
  free(Fkmh);
  free(Zkph);
  free(Zkmh);
  free(gradFk);
  free(dZk);
  free(S);
  std::cout << "END NEWTON\n" << std::endl;
}

#endif
