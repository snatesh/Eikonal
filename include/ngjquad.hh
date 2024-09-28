#ifndef _NGJQUAD_H
#define _NGJQUAD_H

#include<algorithm>
#include<complex.h>
#include<cblas.h>
#include<lapacke.h>
#include<nlopt.h>
#include<jMat.hh>
#include<jPoly.hh>
#include<jevd.hh>

typedef double _Complex Complex;
static unsigned long count = 0;
static unsigned long count1 = 0;


struct optData
{

  double *Vm = 0, *Fk = 0 ;
  double *Jn1 = 0, *Jn2 = 0, *Jn3 = 0, *J = 0;
  Complex *Jnz = 0, *XY0 = 0;
  double *X0 = 0, *Y0 = 0;
  double *X10 = 0, *X20 = 0, *X30 = 0;
  double *Vm_T = 0, *W0 = 0, *S = 0, *F0 = 0;
  jMat<double>* Jn = 0;
  jPoly<double>* Pm = 0;
  double a, b, c, d;
  unsigned int N, M, n, m;
  bool run0; 
  double minf, minf1;
  unsigned int nthreads;
  unsigned int dim;
  double* Z0 = 0; // this holds the optimization variables
  double *legx, *legw; // needed for jmat d=3

  std::string dir;
  std::string splxT;
  
  void init()
  {
  
    initmsg(0);
    if (this->dim == 1)
    {
      this->N = n; this->M = m;
      this->Jn = new jMat<double>(n, a, b);
      this->Pm = new jPoly<double>(N, M, a, b, nthreads);
      this->Vm = Pm->V;
      this->Jn1 = Jn->J;
      this->X0 = (double*) calloc(N, sizeof(double));
      this->Vm_T = (double*) calloc(M*N, sizeof(double));
      this->W0  = (double*) calloc((N >= M ? N : M), sizeof(double));
      W0[0] = 1.0;
      this->S = (double*) calloc(M, sizeof(double)); 
      this->Z0 = (double*) calloc(2*N, sizeof(double));
      this->F0 = (double*) calloc(2*N, sizeof(double)); 
      this->Fk = this->F0;
    }
    else if (this->dim == 2)
    {
      this->N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
      this->M = static_cast<unsigned int>(0.5 * m * (m + 1));
      // generate jacobi matrices for x,y 
      this->Jn = new jMat<double>(n, a, b, c);
      this->Pm = new jPoly<double>(N, m-1, a, b, c, nthreads); 
      this->Vm = Pm->V;
      this->Jn1 = Jn->Jn1;
      this->Jn2 = Jn->Jn2;
      this->Jnz = (Complex*) calloc(N*N, sizeof(Complex));
      this->XY0 = (Complex*) calloc(N, sizeof(Complex));
      this->X0  = (double*) calloc(N, sizeof(double));
      this->Y0  = (double*) calloc(N, sizeof(double));
      this->Vm_T  = (double*) calloc(M*N, sizeof(double));
      this->W0  = (double*) calloc((N >= M ? N : M), sizeof(double));
      W0[0] = 1.0;
      this->S = (double*) calloc(M, sizeof(double)); 
      this->Z0 = (double*) calloc(3*N, sizeof(double));
      this->F0 = (double*) calloc(3*N, sizeof(double)); 
      this->Fk = this->F0;
    }
    else if (this->dim == 3)
    { 
      this->N = dimPI3(n-1); 
      this->M = dimPI3(m-1);
      // generate jacobi matrices for x,y 
      //this->Jn = new jMat<double>(n, a, b, c, d, dir);
      unsigned int nlg = 20; 
      this->Jn = new jMat<double>(n, a, b, c, d, nlg, legx, legw, 1);
      this->Pm = new jPoly<double>(N, m-1, a, b, c, d, nthreads); 
      this->Vm = Pm->V;
      this->Jn1 = Jn->Jn1;
      this->Jn2 = Jn->Jn2;
      this->Jn3 = Jn->Jn3;
      this->J = (double*) calloc(N*N*3, sizeof(double));
      this->X10  = (double*) calloc(N, sizeof(double));
      this->X20  = (double*) calloc(N, sizeof(double));
      this->X30  = (double*) calloc(N, sizeof(double));
      this->Vm_T  = (double*) calloc(M*N, sizeof(double));
      this->W0  = (double*) calloc((N >= M ? N : M), sizeof(double));
      W0[0] = 1.0;
      this->S = (double*) calloc(M, sizeof(double)); 
      this->Z0 = (double*) calloc(4*N, sizeof(double));
      this->F0 = (double*) calloc(4*N, sizeof(double)); 
      this->Fk = this->F0;
    }
    initmsg(1);
  }
   
  optData ( unsigned int _m, unsigned int _n, 
            double _a, double _b, 
            unsigned int _nthreads )
    : n(_n), m(_m), a(_a), b(_b), 
      run0(true), dim(1), 
      nthreads(_nthreads) 
    { 
      if ( m < (2 * n - 1) )
      {
        std::cout << "\nOPTIMAL GAUSS-JACOBI QUADRATURE IS OBTAINABLE in 1D\n";
        std::cout << "N nodes can integrate 2N-1 POLYNOMIALS\n";
        std::cout << "DEFAULTING M = 2N-1\n";
        this->m = 2 * n - 1; 
      } 
      this->init(); 
    }
  
  optData ( unsigned int _m, unsigned int _n, 
            double _a, double _b, double _c, 
            unsigned int _nthreads )
    : n(_n), m(_m), a(_a), b(_b), c(_c), 
      run0(true), dim(2), nthreads(_nthreads) { this->init(); }
  
  optData ( unsigned int _m, unsigned int _n, 
            double _a, double _b, double _c, 
            double _d, std::string _dir, 
            unsigned int _nthreads )
    : n(_n), m(_m), a(_a), b(_b), c(_c), d(_d), dir(_dir), 
      run0(true), dim(3), nthreads(_nthreads) { this->init(); }
  
  optData ( unsigned int _m, unsigned int _n, 
            double _a, double _b, double _c, 
            double _d, double* _legx, double* _legw, 
            unsigned int _nthreads )
    : n(_n), m(_m), a(_a), b(_b), c(_c), d(_d), 
      legx(_legx), legw(_legw), run0(true), dim(3),  
      nthreads(_nthreads) { this->init(); }
  
  ~optData()
  {
    if (Jn) { delete Jn; Jn = 0; }
    if (Pm) { delete Pm; Pm = 0; }
    if (X0) { free(X0); X0 = 0; }
    if (Y0) { free(Y0); Y0 = 0; }
    if (Z0) { free(Z0); Z0 = 0; }
    if (Vm_T) { free(Vm_T); Vm_T = 0; }
    if (S) { free(S); S = 0; }
    if (W0) { free(W0); W0 = 0; }
    if (F0) { free(F0); F0 = 0; }
    if (X10) { free(X10); X10 = 0; }
    if (X20) { free(X20); X20 = 0; }
    if (X30) { free(X30); X30 = 0; }
    if (J)  { free(J); J = 0; }
    if (XY0) { free(XY0); XY0 = 0; }
    if (Jnz) { free(Jnz); Jnz = 0; }
  }

  void initmsg(unsigned int step)
  {
    if (this->dim < 1 || this->dim > 3)
    {   
      std::cerr << "dim must be 1,2 or 3! Exiting ..\n";
      exit(1);
    }
    if (step == 0) 
    {
      splxT = (this->dim == 2 ? "(TRIANGLE)" : 
                (this->dim == 3 ? "(TETRAHEDRON)" : "(LINE)"));
      std::cout << "\nBEGIN NGJQUAD " << splxT << " INITIALIZATION\n"; 
    }
    else if (step == 1)
    {
      std::cout << "\nORDER OF SOURCE BASIS : " << n-1 << std::endl;
      std::cout << "ORDER OF TARGET BASIS : " << m-1 << std::endl;
      std::cout << "\nSEARCHING FOR QUADRATURE RULE OF SIZE N = " << N << std::endl;
      std::cout << "TO EXACTLY INTEGRATE M = " << M << " POLYNOMIALS\n"
                << "WITH TOTAL DEGREE m = " << 0 << " .. " << m - 1 << "\n\n";
      std::cout << "END NGJQUAD " << splxT << " INITIALIZATION\n";
    }
  }

};


struct ngjQuad
{
  optData* optdata = 0;
  unsigned int n, m;
  double tol, tolc, alph;
  bool use_newton, use_wolfe;
  unsigned int nthreads;
  nlopt_algorithm alg;
  unsigned int dim;
 

  ngjQuad ( unsigned int _n, unsigned int _m, 
            double a, double b, 
            double _tol, double _tolc, double _alph,
            bool _use_newton, bool _use_wolfe,
            nlopt_algorithm _alg, unsigned int _nthreads)
    : n(_n), m(_m), tol(_tol), tolc(_tolc), alph(_alph),
      use_newton(_use_newton), use_wolfe(_use_wolfe),
      alg(_alg), nthreads(_nthreads), dim(1)
  {
    this->optdata = new optData(m, n, a, b, nthreads);
  }

  
  ngjQuad ( unsigned int _n, unsigned int _m, 
            double a, double b, double c, 
            double _tol, double _tolc, double _alph,
            bool _use_newton, bool _use_wolfe,
            nlopt_algorithm _alg, unsigned int _nthreads)
    : n(_n), m(_m), tol(_tol), tolc(_tolc), alph(_alph),
      use_newton(_use_newton), use_wolfe(_use_wolfe),
      alg(_alg), nthreads(_nthreads), dim(2)
  {
    this->optdata = new optData(m, n, a, b, c, nthreads);
  }
  
  ngjQuad ( unsigned int _n, unsigned int _m, 
            double a, double b, double c, double d, 
            double _tol, double _tolc, double _alph,
            bool _use_newton, bool _use_wolfe,
            nlopt_algorithm _alg, unsigned int _nthreads, 
            std::string dir)
    : n(_n), m(_m), tol(_tol), tolc(_tolc), alph(_alph),
      use_newton(_use_newton), use_wolfe(_use_wolfe),
      alg(_alg), nthreads(_nthreads), dim(3)
  {
    this->optdata = new optData(m, n, a, b, c, d, dir, nthreads);
  }
  
  ngjQuad ( unsigned int _n, unsigned int _m, 
            double a, double b, double c, double d,
            double* legx, double* legw, 
            double _tol, double _tolc, double _alph,
            bool _use_newton, bool _use_wolfe,
            nlopt_algorithm _alg, unsigned int _nthreads )
    : n(_n), m(_m), tol(_tol), tolc(_tolc), alph(_alph),
      use_newton(_use_newton), use_wolfe(_use_wolfe),
      alg(_alg), nthreads(_nthreads), dim(3)
  {
    this->optdata = new optData(m, n, a, b, c, d, legx, legw, nthreads);
  }
  
  ~ngjQuad() 
  { 
    count = count1 = 0;
    if (optdata) { delete optdata; optdata = 0; }
  }
  
  
  static inline void cond ( double* Vm, unsigned int N, 
                            unsigned int M, double* rcond, 
                            bool norm2=true )
  {
    // 2-norm cond
    if (norm2)
    {
      unsigned int dimS = (N <= M) ? N : M;
      double* S = (double*) calloc(dimS, sizeof(double));
      double* superb = (double*) calloc(dimS, sizeof(double));
  
      LAPACKE_dgesvd( LAPACK_COL_MAJOR, 'N','N', N, M, Vm, N, S, 
                      nullptr, 1, nullptr, 1, superb  );
      rcond[0] = S[dimS-1] / S[0];
  
      free(S); free(superb);
    }
    // inf norm cond
    else
    {
      int* ipiv = (int*) calloc(N, sizeof(int));
      double normVm = LAPACKE_dlange(LAPACK_COL_MAJOR, '1', N, M, Vm, N);
      // LU of A
      LAPACKE_dgetrf(LAPACK_COL_MAJOR, N, M, Vm, N, ipiv); 
      LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', N, Vm, M, normVm, rcond);
      free(ipiv);
    }
  }
  
  static inline double optF ( unsigned int n, const double* Zk, 
                              double* grad, void* _data  )
  {
    ++count;
    optData* data = (optData*) _data;
    unsigned int m = data->m;
    unsigned int M = data->M;
    unsigned int N = data->N;
    unsigned int dim = data->dim;
    double *Vm = data->Pm->V, *Fk = data->Fk;
    const double *Wk, *Xk, *Yk; 
    const double *X1k, *X2k, *X3k;
    if (dim == 1)
    {
      Xk = Zk; Wk = Zk + N;
      Fk = data->Fk;
      data->Pm->computeV(Xk);
    } 
    else if (dim == 2)
    {
      Xk = Zk; Yk = Zk + N; 
      Wk = Zk + 2*N;
      Fk = data->Fk;
      data->Pm->computeV(Xk, Yk);
    }
    else if (dim == 3)
    {
      X1k = Zk; X2k = Zk + N; X3k = Zk + 2*N;
      Wk = Zk + 3*N;
      data->Pm->computeV(X1k, X2k, X3k);
    }
    double* I0 = (double*) calloc(M, sizeof(double)); I0[0] = 1;  
    cblas_dgemv( CblasColMajor, CblasTrans,  N, M, 1.0, Vm, N, Wk, 1, -1.0, I0, 1 );
    #pragma omp simd 
    for (unsigned int i = 0; i < M; ++i) { Fk[i] = I0[i]; }
    free(I0);
    double norm2;
    norm2 = pow(cblas_dnrm2(M, Fk, 1), 2);
    if ( !(count % 100000) && data->run0 ) 
    { 
      std::cout << "Eval #" << count << " : F = " << norm2 << std::endl; 
    }
    return norm2;
  }
  
  static inline double optF1  ( unsigned int n, const double* Zk, 
                                double* grad, void* _data  )
  {
    ++count1;
    optData* data = (optData*) _data;
    unsigned int dim = data->dim; 
    if (dim == 1)
    {
      const double *Xk = Zk;
      data->Pm->computeV(Xk);
    }
    else if (dim == 2) 
    { 
      const double *Xk = Zk; 
      const double *Yk = Zk + data->N;
      data->Pm->computeV(Xk, Yk); 
    }
    else if (dim == 3) 
    {
      const double *X1 = Zk;
      const double* X2 = Zk + data->N;
      const double* X3 = Zk + 2 * data->N; 
      data->Pm->computeV(X1, X2, X3); 
    }
    double rcond[1];
    cond(data->Pm->V, data->N, data->M, rcond);
    if ( !(count1 % 100000) ) 
    { 
      std::cout << "Eval #" << count1 << " : F = " << 1.0 / rcond[0] << std::endl; 
    }
    return 1.0 / rcond[0];
  }
  
  
  static inline void optieqC (  unsigned int m, double* result, unsigned int n, 
                                const double* Zk, double* grad, 
                                void* _data )
  {
   
    optData* data = (optData*) _data;
    unsigned int dim = data->dim;
    if (dim == 1)
    {
      // only bound constraints
    }  
    else if (dim == 2)
    {
      unsigned int N = static_cast<unsigned int>(n / 3.0);
      #pragma omp simd
      for (unsigned int i = 0; i < m-1; ++i)
      {
        result[i] = Zk[i] + Zk[i + N] - 1.0; 
      }
    }
    else if (dim == 3)
    {
      unsigned int N = static_cast<unsigned int>(n / 4.0);
      #pragma omp simd
      for (unsigned int i = 0; i < m-1; ++i)
      {
        result[i] = Zk[i] + Zk[i + N] + Zk[i + 2*N] - 1.0; 
      }
    }
  }
  
  static inline void optieqC1 ( unsigned int m, double* result, unsigned int n, 
                                const double* Zk, double* grad, void* f_data)
  {
    optData* data = (optData*) f_data;
    unsigned int dim = data->dim;
    if (dim == 1)
    {
      // only bound constraints
    }
    else if (dim == 2)
    {
      unsigned int N = static_cast<unsigned int>(n / 2.0);
      #pragma omp simd 
      for (unsigned int i = 0; i < m-1; ++i)
      {
        result[i] = Zk[i] + Zk[i + N] - 1.0; 
      }
    }
    else if (dim == 3)
    {
      unsigned int N = static_cast<unsigned int>(n / 3.0);
      #pragma omp simd 
      for (unsigned int i = 0; i < m-1; ++i)
      {
        result[i] = Zk[i] + Zk[i + N] + Zk[i + 2*N] - 1.0; 
      }
    }
  }
  
  static inline double opteqC ( unsigned int n, const double* Zk, 
                                double* grad, void* _data ) 
  {
    double sumw = 0.0;
    optData* data = (optData*) _data;
    unsigned int dim = data->dim;
    unsigned int N = static_cast<unsigned int>(n / (dim+1.0));
    #pragma omp simd reduction(+:sumw)
    for (unsigned int i = dim*N; i < n; ++i) { sumw += Zk[i]; }
    return sumw - 1.0;
  } 
  
  static inline double opteqC1  ( unsigned int n, const double* Zk,
                                  double* grad, void* _data) 
  {
    optData* data = (optData*) _data;
    data->run0 = false;
    return optF(data->n, Zk, nullptr, data) - data->minf;
  }  
  
  inline void F() 
  {
    unsigned int N = optdata->N;
    unsigned int M = optdata->M;
    unsigned int m = optdata->m;
    double* Zk = optdata->Z0;
    double* Fk = optdata->F0;
    double* Wk;
    double* I0 = (double*) calloc(M, sizeof(double)); I0[0] = 1; 
    
    if (this->dim == 1)
    {
      const double *Xk;
      Xk = Zk; Wk = Zk + N;
      optdata->Pm->computeV(Xk);
    }
    else if (this->dim == 2)
    {
      const double *Xk, *Yk;
      Xk = Zk; Yk = Zk + N; Wk = Zk + 2*N;
      optdata->Pm->computeV(Xk, Yk);
    }
    else if (this->dim == 3)
    {
      const double *X1k, *X2k, *X3k;
      X1k = Zk; X2k = Zk + N; X3k = Zk + 2*N; Wk = Zk + 3*N;
      optdata->Pm->computeV(X1k, X2k, X3k);
    }
    cblas_dgemv ( CblasColMajor, CblasTrans, N, M, 1.0, optdata->Pm->V, 
                  N, Wk, 1, -1.0, I0, 1);
    for (unsigned int i = 0; i < M; ++i) { Fk[i] = I0[i]; }
    free(I0);

  }
  
  inline void init ()
  {
    unsigned int m = optdata->m; unsigned int n = optdata->n; 
    unsigned int M = optdata->M; unsigned int N = optdata->N;
    double* Vm = optdata->Vm; double* Vm_T = optdata->Vm_T; 
    double* W0 = optdata->W0; double* S = optdata->S; 
    double* Z0 = optdata->Z0; double* F0 = optdata->F0; 
    if (this->dim == 3)
    {
      double* Jn1 = optdata->Jn1; 
      double* Jn2 = optdata->Jn2;
      double* Jn3 = optdata->Jn3;
      double* J = optdata->J;
      double a = optdata->a, b = optdata->b, c = optdata->c, d = optdata->d;
      for (unsigned int j = 0; j < N; ++j)
      {
        for (unsigned int i = 0; i < N; ++i)
        {
          J[i + N*j]          = Jn1[i + j*N];
          J[i + N*(j + N)]    = Jn2[i + j*N];
          J[i + N*(j + 2*N)]  = Jn3[i + j*N];
        }
      }
      // compute initial nodes and weights from jevd of J=[Jx,Jy,Jz] 
      jointDiag<double>* jevd = new jointDiag<double>(N, this->dim, 1e-10, J, 1); 
  
      for (unsigned int i = 0; i < N; ++i) 
      { 
        optdata->X10[i] = J[i + N*i] ; 
        optdata->X20[i] = J[i + N*(i + N)] ;
        optdata->X30[i] = J[i + N*(i + 2*N)] ; 
      }
      // evaluate Vandermonde on initial nodes
      optdata->Pm->computeV(optdata->X10, optdata->X20, optdata->X30); 
    }
    else if (this->dim == 2)
    {
      double* Jn1 = optdata->Jn1; double* Jn2 = optdata->Jn2;
      Complex* XY0 = optdata->XY0;
      // compute initial nodes and weights from eigenvalues of Jn
      for (unsigned int i = 0; i < N*N; ++i) { optdata->Jnz[i] = Jn1[i] + I*Jn2[i]; }
  
      if (LAPACKE_zgeev ( LAPACK_COL_MAJOR, 'N', 'N', N, optdata->Jnz, N, 
                          XY0, nullptr, N, nullptr, N ))
      {
        std::cerr << "ERROR: NGJQUAD INIT" << std::endl;
      }
  
      for (unsigned int i = 0; i < N; ++i) 
      { 
        optdata->X0[i] = creal(XY0[i]); 
        optdata->Y0[i] = cimag(XY0[i]); 
      } 
      // evaluate Vandermonde on initial nodes
      optdata->Pm->computeV(optdata->X0, optdata->Y0);
    }
    else if (this->dim == 1)
    {
      double* X0 = optdata->X0; 
      double* X0i = (double*) calloc(N, sizeof(double));
      // compute initial nodes and weights from eigenvalues of Jn
      if (LAPACKE_dgeev ( LAPACK_COL_MAJOR, 'N', 'N', N, optdata->Jn1, N, 
                          X0, X0i, nullptr, N, nullptr, N ))
      {
        std::cerr << "ERROR: NGJQUAD INIT" << std::endl;
      }
      free(X0i);
      // evaluate Vandermonde on initial nodes
      std::sort(optdata->X0, optdata->X0 + N);
      optdata->Pm->computeV(optdata->X0);
    }
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
    if (LAPACKE_dgelsd  ( LAPACK_COL_MAJOR, M, N, 1, Vm_T, 
                          M, W0, (N >= M ? N : M), S, 
                          -1.0, rank ))
    {
      std::cerr << "ERROR: Lapack *gelsd: Pseudoinverse" << std::endl;
    }
  
    if (this->dim == 3)
    {
      // copy into Z0 for opt routines
      for (unsigned int i = 0; i < N; ++i)
      {
        Z0[i]       = optdata->X10[i];
        Z0[i + N]   = optdata->X20[i];
        Z0[i + 2*N] = optdata->X30[i];
        Z0[i + 3*N] = W0[i];
      }
  
    }
    else if (this->dim == 2)
    {
      // copy into Z0 for opt routines
      for (unsigned int i = 0; i < N; ++i)
      {
        Z0[i]       = optdata->X0[i];
        Z0[i + N]   = optdata->Y0[i];
        Z0[i + 2*N] = W0[i];
      }
    }
    else if (this->dim == 1)
    {
      for (unsigned int i = 0; i < N; ++i)
      {
        Z0[i]       = optdata->X0[i];
        Z0[i + N]   = W0[i];
      }
    }
    F();
  }
  
  
  inline void runXW ()
  {
    double minF = cblas_dnrm2(optdata->M, optdata->F0, 1); 
    if (minF > tolc)
    {
      std::cout <<"\n BEGIN NLOPT (X,W) \n"; 
      unsigned int N = optdata->N; 
      nlopt_opt opt;
      double *lb, *ub, *tolieq;
      if (this->dim == 1)
      {
        // x > -1, w > 0
        lb = (double*) calloc(2*N, sizeof(double));
        for (unsigned int i = 0; i < N; ++i) { lb[i] = -1; } 
        // x, w < 1
        ub = (double*) calloc(2*N, sizeof(double));
        tolieq = (double*) calloc(1*N, sizeof(double));
        for (unsigned int i = 0; i < 2*N; ++i) { ub[i] = 1; }
        for (unsigned int i = 0; i < 1*N; ++i) { tolieq[i] = tolc; }
        opt = nlopt_create(alg, 2*N); 
      }
      else if (this->dim == 2)
      {
        // x, y , w > 0
        lb = (double*) calloc(3*N, sizeof(double)); 
        // x, y , w < 1
        ub = (double*) calloc(3*N, sizeof(double));
        tolieq = (double*) calloc(2*N, sizeof(double));
        for (unsigned int i = 0; i < 3*N; ++i) { ub[i] = 1; }
        for (unsigned int i = 0; i < 2*N; ++i) { tolieq[i] = tolc; }
        opt = nlopt_create(alg, 3*N); 
        nlopt_add_inequality_mconstraint(opt, 2*N, optieqC, optdata, tolieq);
      }
      else if (this->dim == 3)
      {
        // x, y, z, w > 0
        lb = (double*) calloc(4*N, sizeof(double)); 
        // x, y, z, w < 1
        ub = (double*) calloc(4*N, sizeof(double));
        tolieq = (double*) calloc(3*N, sizeof(double));
        for (unsigned int i = 0; i < 4*N; ++i) { ub[i] = 1; }
        for (unsigned int i = 0; i < 3*N; ++i) { tolieq[i] = tolc; }
        opt = nlopt_create(alg, 4*N); 
        nlopt_add_inequality_mconstraint(opt, 3*N, optieqC, optdata, tolieq);
  
      }
      nlopt_set_lower_bounds(opt, lb);
      nlopt_set_upper_bounds(opt, ub);
      nlopt_set_min_objective(opt, optF, optdata);
      nlopt_add_equality_constraint(opt, opteqC, optdata, tolc);
      nlopt_set_xtol_rel(opt, tol);
      nlopt_set_stopval(opt, tol);

      if (nlopt_optimize(opt, optdata->Z0, &minF) < 0) 
      {
        std::cerr << "NLOPT failed! Exiting .." << std::endl;
        exit(1);
      }
      nlopt_destroy(opt);
      free(lb); free(ub);
      free(tolieq);
      optdata->minf = minF;
    } 
    std::cout << "Quadrature exactness in norm up to specified order : " 
              << minF << std::endl;
    std::cout <<"END NLOPT \n\n";
  }
  
  inline void runX ()
  {
    std::cout <<"\n BEGIN NLOPT (X)\n"; 
    unsigned int N = optdata->N;
    nlopt_opt opt; 
    double *ub, *lb, *tolieq;
    if (this->dim == 1)
    {
      // x > -1
      lb = (double*) calloc(1*N, sizeof(double));
      // x < 1
      ub = (double*) calloc(1*N, sizeof(double));
      tolieq = (double*) calloc(1*N, sizeof(double));
      for (unsigned int i = 0; i < 1*N; ++i) { lb[i] = -1; ub[i] = 1; }
      for (unsigned int i = 0; i < 1*N; ++i) { tolieq[i] = tolc; }
      opt = nlopt_create(alg, 1*N); 
    }
    else if (this->dim == 2)
    { 
      // x, y > 0
      lb = (double*) calloc(2*N, sizeof(double)); 
      // x, y < 1
      ub = (double*) calloc(2*N, sizeof(double));
      tolieq = (double*) calloc(2*N, sizeof(double));
      for (unsigned int i = 0; i < 2*N; ++i) { ub[i] = 1; }
      for (unsigned int i = 0; i < 2*N; ++i) { tolieq[i] = tolc; }
      opt = nlopt_create(alg, 2*N); 
      nlopt_add_inequality_mconstraint(opt, 2*N, optieqC1, optdata, tolieq);
    }
    else if (this->dim == 3)
    {
      // x, y, z > 0
      lb = (double*) calloc(3*N, sizeof(double)); 
      // x, y ,z < 1
      ub = (double*) calloc(3*N, sizeof(double));
      tolieq = (double*) calloc(3*N, sizeof(double));
      for (unsigned int i = 0; i < 3*N; ++i) { ub[i] = 1; }
      for (unsigned int i = 0; i < 3*N; ++i) { tolieq[i] = tolc; }
      opt = nlopt_create(alg, 3*N); 
      nlopt_add_inequality_mconstraint(opt, 3*N, optieqC1, optdata, tolieq);
    } 
  
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, optF1, optdata);
    nlopt_add_equality_constraint(opt, opteqC1, optdata, tolc);
  
    nlopt_set_xtol_rel(opt, 1e-8);
    nlopt_set_stopval(opt, 1);
    double minF;
    if (nlopt_optimize(opt, optdata->Z0, &minF) < 0) 
    {
      std::cerr << "NLOPT failed!" << std::endl;
    }
    else 
    {
      std::cout << "Conditioning of interpolation operator on new abscissa : " 
                << minF << std::endl;
    }
    nlopt_destroy(opt);
    free(lb); free(ub);
    optdata->minf1 = minF;
    free(tolieq);
    std::cout <<"END NLOPT \n\n";
  }

  inline void newton() 
  {
    std::cout << "BEGIN NEWTON" << std::endl;
    double h = 1e-8, tol = 1e-7; 
    double tol_up = 1e5; 
    unsigned int maxiter = 10000;
    unsigned int M = optdata->M;
    unsigned int N = optdata->N;

    double pk = cblas_dnrm2(M, optdata->F0, 1);
    double* Fkph    = (double*) calloc(M, sizeof(double));
    double* Fkmh    = (double*) calloc(M, sizeof(double));
    double* Zkph    = (double*) calloc((dim+1)*N, sizeof(double));
    double* Zkmh    = (double*) calloc((dim+1)*N, sizeof(double));
    double* gradFk  = (double*) calloc(M*(dim+1)*N, sizeof(double)); 
    double* dZk     = (double*) calloc((dim+1)*N, sizeof(double)); 
    double* S       = (double*) calloc(M, sizeof(double));
  
    double rho = 0.9; double gam = 1e-4;
    unsigned int iter = 0; lapack_int rank[1];
    double* Zk = optdata->Z0;
    while (pk > tol && pk < tol_up && iter < maxiter)
    {
      iter += 1;
      /* dZk initialized to Fk, overwritten by dgelsd to dZk 
         which is lsq sol to gradFk*dZk = Fk */
      for (unsigned int i = 0; i < M; ++i)    { dZk[i] = optdata->Z0[i]; }
      for (unsigned int i = M; i < (dim+1)*N; ++i)  { dZk[i] = 0; }
      for (unsigned int i = 0; i < (dim+1)*N; ++i) { Zkph[i] = Zkmh[i] = Zk[i]; }
      // save current Z0 address
      Zk = optdata->Z0;
      // compute finite difference approx to grad
      for (unsigned int jj = 0; jj < (dim+1)*N; ++jj)
      {
        // eval above and below Zk
        Zkph[jj] += h; Zkmh[jj] -= h;
        // switch pointers, call F
        optdata->Z0 = Zkph; F();
        for (unsigned int ii = 0; ii < M; ++ii) { Fkph[ii] = optdata->F0[ii]; }
        // switch pointers, call F
        optdata->Z0 = Zkmh; F();
        for (unsigned int ii = 0; ii < M; ++ii) { Fkmh[ii] = optdata->F0[ii]; }
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
      if (  LAPACKE_dgelsd  ( LAPACK_COL_MAJOR, M, (dim+1)*N, 1, 
                              gradFk, M, dZk, (dim+1)*N, S, 1e-16, rank ) )
      {
        std::cerr << "ERROR: Lapack dgelsd: Pseudoinverse" << std::endl;
      }
      // linesearch with Wolfe conditions
      double wfnrm;
      for (unsigned int i = 0; i < (dim+1)*N; ++i) { Zk[i] -= alph*dZk[i]; }
      // restore Zk after derivative computations
      optdata->Z0 = Zk; F();
      pk = cblas_dnrm2(M, optdata->F0, 1);
      if (use_wolfe)
      {
        for (unsigned int iter_i = 0; iter_i < maxiter; ++iter_i)
        {
          cblas_dgemv ( CblasColMajor,  CblasNoTrans, M, (dim+1)*N, 
                        gam*alph, gradFk, M, dZk, 1, 1.0, optdata->F0, 1  );
          wfnrm = cblas_dnrm2(M, optdata->F0, 1);
          if (pk > wfnrm)
          {
            alph = rho*alph;
            for (unsigned int i = 0; i < (dim+1)*N; ++i) { optdata->Z0[i] -= alph*dZk[i]; }
            F();
            pk = cblas_dnrm2(M, optdata->F0, 1);
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

};



#endif
