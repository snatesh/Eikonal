#ifndef _NGJQUAD_REDOBJ_H
#define _NGJQUAD_REDOBJ_H

#include<algorithm>
#include<complex.h>
#include<cblas.h>
#include<lapacke.h>
#include<nlopt.h>
#include<jMat.hh>
#include<dMat.hh>
#include<jPoly.hh>
#include<jevd.hh>

typedef double _Complex Complex;
static unsigned long count = 0;
static unsigned long count1 = 0;


struct optData
{

  double *Vm = 0;
  double *Jn1 = 0, *Jn2 = 0;
  Complex *Jnz = 0, *XY0 = 0;
  double *X0 = 0, *Y0 = 0;
  double *Q = 0, *W0 = 0;
  double* one = 0;
  jMat<double>* Jn = 0;
  jPoly<double>* Pm = 0;
  jPoly<double>* Pmx = 0;
  jPoly<double>* Pmy = 0;
  double a, b, c, d;
  unsigned int N, M, n, m;
  bool run0; 
  double minf, minf1;
  unsigned int nthreads;
  unsigned int dim;
  double* Z0 = 0; 
  double* G = 0;
  double *Habc = 0, *Ha1bc1 = 0, *Hab1c1 = 0;
  double *Dx = 0, *Dy = 0;
  double *dQ = 0, *dqstor = 0, *dgQ = 0;
  double *dgQdq = 0;

  std::string dir;
  std::string splxT;
  
  void init()
  {
  
    initmsg(0);
    if (this->dim == 2)
    {
      this->N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
      this->M = static_cast<unsigned int>(0.5 * m * (m + 1));
      // generate jacobi matrices for x,y 
      this->Jn = new jMat<double>(n, a, b, c);
      this->Pm = new jPoly<double>(N, m-1, a, b, c, nthreads); 
      this->Pmx = new jPoly<double>(N, m-1, a+1, b, c+1, nthreads); 
      this->Pmy = new jPoly<double>(N, m-1, a, b+1, c+1, nthreads); 
      this->Habc  = (double*) calloc((m+1)*(m+1), sizeof(double));
      this->Ha1bc1  = (double*) calloc((m+1)*(m+1), sizeof(double));
      this->Hab1c1 = (double*) calloc((m+1)*(m+1), sizeof(double));
      this->Dx = (double*) calloc(M*M, sizeof(double));
      this->Dy = (double*) calloc(M*M, sizeof(double));
      sFactors(m+1, a, b, c, Habc);  
      sFactors(m+1, a+1, b, c+1, Ha1bc1);  
      sFactors(m+1, a, b+1, c+1, Hab1c1);  
      dMat(a, b, c, Habc, Ha1bc1, m-1, 0, Dx);
      dMat(a, b, c, Habc, Hab1c1, m-1, 1, Dy); 
      this->Vm = Pm->V;
      this->Jn1 = Jn->Jn1;
      this->Jn2 = Jn->Jn2;
      this->Jnz = (Complex*) calloc(N*N, sizeof(Complex));
      this->XY0 = (Complex*) calloc(N, sizeof(Complex));
      this->X0  = (double*) calloc(N, sizeof(double));
      this->Y0  = (double*) calloc(N, sizeof(double));
      this->Q   = (double*) calloc(N*N, sizeof(double));  
      this->W0  = (double*) calloc(N, sizeof(double));
      this->one = (double*) calloc(N, sizeof(double));
      this->dQ  = (double*) calloc(N*dim*N*N, sizeof(double));
      this->dgQ = (double*) calloc(N*N, sizeof(double));
      this->dqstor = (double*) calloc(M, sizeof(double));
      this->dgQdq = (double*) calloc(N*N, sizeof(double));
      for (unsigned int i = 0; i < N; ++i) { one[i] = 1.0; }
      this->Z0 = (double*) calloc(dim*N, sizeof(double));
      //this->G = (double*) calloc((dim+1)*N*(dim*N), sizeof(double));
      this->G = (double*) calloc(N*(dim*N), sizeof(double));
      for (unsigned int i = 0; i < N; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j) 
        {
          G[i + N*(i + j*N)] = 1.0;
        }
        //G[i + N*i] = 1.0;
        //G[i + N*(i+N)] = 1.0;
      }
      //for (unsigned int i = 0; i < dim*N; ++i)
      //{
      //  G[i + i*((dim+1)*N)] = -1.0;
      //}   
      //for (unsigned int i = dim*N; i < (dim+1)*N; ++i)
      //{
      //  for (unsigned int j = 0; j < dim; ++j)
      //  {   
      //    G[i + ((i-dim*N)+j*N)*((dim+1)*N)] = 1.0;
      //  }   
      //}
    }
    else
    {
      std::cerr << "only dim=2 is supported. Exiting..\n";
      exit(1);
    }
    initmsg(1);
  }
   
  
  optData ( unsigned int _m, unsigned int _n, 
            double _a, double _b, double _c, 
            unsigned int _nthreads )
    : n(_n), m(_m), a(_a), b(_b), c(_c), 
      run0(true), dim(2), nthreads(_nthreads) { this->init(); }
  
  ~optData()
  {
    if (Jn) { delete Jn; Jn = 0; }
    if (Pm) { delete Pm; Pm = 0; }
    if (Pmx) { delete Pmx; Pmx = 0; }
    if (Pmy) { delete Pmy; Pmy = 0; }
    if (X0) { free(X0); X0 = 0; }
    if (Y0) { free(Y0); Y0 = 0; }
    if (Z0) { free(Z0); Z0 = 0; }
    if (Q) { free(Q); Q = 0; }
    if (W0) { free(W0); W0 = 0; }
    if (XY0) { free(XY0); XY0 = 0; }
    if (Jnz) { free(Jnz); Jnz = 0; }
    if (G)  { free(G); G = 0; }
    if (Habc) { free(Habc); Habc = 0; }
    if (Ha1bc1) { free(Ha1bc1); Ha1bc1 = 0; }
    if (Hab1c1) { free(Hab1c1); Hab1c1 = 0; }
    if (Dx)   { free(Dx); Dx = 0; }
    if (Dy)   { free(Dy); Dy = 0; }
    if (one) { free(one); one = 0; }
    if (dQ) { free(dQ); dQ = 0; }
    if (dqstor) { free(dqstor); dqstor = 0; }
    if (dgQ)  { free(dgQ); dgQ = 0; }
    if (dgQdq) { free(dgQdq); dgQdq = 0; }
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
  
  inline void dPTP ()
  {
    double* X = Z0; double* Y = Z0 + N;
    if (this->dim == 2)
    {
      double *V[2], *D[2]; 
      V[0] = Pmx->V; V[1] = Pmy->V;
      D[0] = Dx; D[1] = Dy;
      Pm->computeV(X, Y);
      Pmx->computeV(X, Y);
      Pmy->computeV(X, Y);
      double dqkldzij;
      for (unsigned int i = 0; i < N; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
        {
          for (unsigned int k = 0; k < N; ++k)
          {
            for (unsigned int l = 0; l < N; ++l)
            {
              if (i != l && i != k) 
              { 
                dqkldzij = 0; 
              }
              else if (i == l && i != k) 
              { 
                cblas_dgemv ( CblasColMajor, CblasTrans, M, M, 1.0, 
                              D[j], M, &(V[j][l]), N, 
                              0.0, dqstor, 1 );
                dqkldzij = cblas_ddot(M, &Vm[k], N, dqstor, 1);
              }
              else if (i != l && i == k) 
              { 
                cblas_dgemv ( CblasColMajor, CblasTrans, M, M, 1.0, 
                              D[j], M, &(V[j][k]), N, 
                              0.0, dqstor, 1 );
                dqkldzij = cblas_ddot(M, &Vm[l], N, dqstor, 1);
              }
              else if (i == l && i == k) 
              { 
                cblas_dgemv ( CblasColMajor, CblasTrans, M, M, 1.0, 
                              D[j], M, &(V[j][k]), N, 
                              0.0, dqstor, 1 );
                dqkldzij = 2.0 * cblas_ddot(M, &Vm[k], N, dqstor, 1);
              }
              //dQ[i + N*(j + dim*(k + N*l))] = dqkldzij; 
              dQ[l + N*(k + N*(j + dim*i))] = dqkldzij;
            }
          }
        }
      }
    }
    else
    {
      std::cerr << "Only dim=2 is supported. Exiting..\n";
      exit(1);
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

  // optimize quadrature problem 
  static inline double optF_help (  unsigned int n, const double* Zk, 
                                    double* grad, void* _data  )
  {
    ++count;
    optData* data = (optData*) _data;
    unsigned int m = data->m;
    unsigned int M = data->M;
    unsigned int N = data->N;
    unsigned int dim = data->dim;
    double *Vm = data->Pm->V;
    const double *Xk, *Yk;
    if (dim == 2)
    {
      Xk = Zk; Yk = Zk + N; 
      data->Pm->computeV(Xk, Yk);
      // compute V*V'
      cblas_dgemm ( CblasColMajor, 
                    CblasNoTrans, 
                    CblasTrans, 
                    N, N, M, 
                    1.0, Vm, N, 
                    Vm, N,  
                    0.0, data->Q, N );
      // perform cholesky Q = U'*U
      // with upper triangular U stored
      // in upper triangle block of Q
      int info = LAPACKE_dpotrf (LAPACK_COL_MAJOR, 'U', N, data->Q, N);
      if ( !info )
      {
        // if cholesky successfull, use factors to compute Q\1
        LAPACKE_dpotrs  ( LAPACK_COL_MAJOR, 'U', N, 1,
                          data->Q, N, data->one, N);
      }
      else
      {
        // if cholesky not succesfull, we have non-spd Q
        // at the current iterate, so we use a 
        // general symmetric solver
        std::cerr << "Could not compute cholesky factorization of Q. \n";
        std::cerr << "Info code: " << info << std::endl;
        // recompute Q
        cblas_dgemm ( CblasColMajor, 
                      CblasNoTrans, 
                      CblasTrans, 
                      N, N, M, 
                      1.0, Vm, N, 
                      Vm, N,  
                      0.0, data->Q, N );
        lapack_int* ipiv = (lapack_int*) calloc(N, sizeof(lapack_int));
        int info1 = LAPACKE_dsysv ( LAPACK_COL_MAJOR, 'L', N, 1, 
                                    data->Q, N, ipiv, data->one, N );
        if (info1 > 0) { std::cerr << "Could not compute inv(Q)*1\n"; exit(1); }
        free(ipiv); 
      }
      // copy sol to proper location, reset rhs
      for (unsigned int i = 0; i < N; ++i)
      {
        data->W0[i] = data->one[i];
        data->one[i] = 1.0;
      }
      // compute 1*W (sum of weights)
      double sum = 0.0;
      #pragma omp simd reduction(+:sum)
      for (unsigned int i = 0; i < N; ++i) { sum += data->W0[i]; }
      double fval = -0.5 * (sum - 1.0);     
      if ( !(count % 5) && data->run0 ) 
      { 
        std::cout << "Eval #" << count << " : F = " << fval << std::endl; 
      }
 
      return fval;
    
    }
    else
    {
      std::cerr << "Only dim=2 is supported. Exiting..\n";
      exit(1);
    }
  
  }
 
  // optimize quadrature problem
  static inline double optF ( unsigned int n, const double* Zk,
                              double* grad, void* _data )
  {
    
    optData* data = (optData*) _data;
    double fval = optF_help(n, Zk, grad, _data);
    if (grad)
    {

      unsigned int N = data->N;
      unsigned int dim = data->dim;
      double* Q = data->Q; 
      double* dQ = data->dQ;
      double* dgQ = data->dgQ;
      // compute dQ
      data->dPTP();
      // compute dg(Q)
      double* Qinvone = data->W0;
      cblas_dger (CblasColMajor, N, N, 1.0, Qinvone, 1, Qinvone, 1, dgQ, N);
      double* dq;
      // compute dgdz
      double* dgQdq = data->dgQdq;
      for (unsigned int i = 0; i < N; ++i)
      {
        for (unsigned int j = 0; j < dim; ++j)
        {
          dq = &(dQ[N*N*(j + dim*i)]);  
          cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans, 
                        N, N, N, 1.0, dgQ, 
                        N, dq, N, 0.0, dgQdq, N );
          double sum = 0.0;
          #pragma omp simd reduction(+:sum)
          for (unsigned int i = 0; i < N; ++i)
          {
            sum += dgQdq[i + i*N];
          }
          grad[i + j*N] = 0.5 * sum; 
        } 
      } 
      // reset dgQ
      for (unsigned int i = 0; i < N*N; ++i) { dgQ[i] = 0; }
    }
    return fval;
  }
 
  // optimize condition number of interpolation matrix on abscissa 
  static inline double optF1  ( unsigned int n, const double* Zk, 
                                double* grad, void* _data  )
  {
    ++count1;
    optData* data = (optData*) _data;
    unsigned int dim = data->dim; 
    if (dim == 2) 
    { 
      const double *Xk = Zk; 
      const double *Yk = Zk + data->N;
      data->Pm->computeV(Xk, Yk); 
    }
    else
    {
      std::cerr << "Only dim=2 supported. Exiting..\n";
      exit(1);
    }
    double rcond[1];
    cond(data->Pm->V, data->N, data->M, rcond);
    if ( !(count1 % 100) ) 
    { 
      std::cout << "Eval #" << count1 << " : cond(V(Z)) = " << 1.0 / rcond[0] << std::endl; 
    }
    return 1.0 / rcond[0];
  }
  
  static inline void optieqC (  unsigned int m, double* result, unsigned int n, 
                                const double* Zk, double* grad, 
                                void* _data )
  {
    optData* data = (optData*) _data;
    unsigned int dim = data->dim;
    if (dim == 2)
    {
      unsigned int N = data->N;
      #pragma omp simd
      for (unsigned int i = 0; i < N; ++i)
      {
        result[i] = Zk[i] + Zk[i + N] - 1.0 + 1e-8; 
      }
      if (grad)
      {
        for (unsigned int j = 0; j < m; ++j)
        {
          for (unsigned int i = 0; i < n; ++i)
          {
            grad[i + n*j] = data->G[j + m*i];
          }
        }
      }
    }
  }
  
  static inline void optieqC1 ( unsigned int m, double* result, unsigned int n, 
                                const double* Zk, double* grad, void* f_data)
  {
    optData* data = (optData*) f_data;
    unsigned int dim = data->dim;
    if (dim == 2)
    {
      unsigned int N = data->N;
      #pragma omp simd 
      for (unsigned int i = 0; i < m; ++i)
      {
        result[i] = Zk[i] + Zk[i + N] - 1.0 + 1e-8; 
      }
    }
    else
    {
      std::cerr << "Only dim=2 is supported. Exiting..\n";
      exit(1);
    }
  }
  
  static inline double opteqC1  ( unsigned int n, const double* Zk,
                                  double* grad, void* _data) 
  {
    optData* data = (optData*) _data;
    data->run0 = false;
    return optF(n, Zk, nullptr, data) - data->minf;
  }  
  
  inline void init ()
  {
    unsigned int m = optdata->m; unsigned int n = optdata->n; 
    unsigned int M = optdata->M; unsigned int N = optdata->N;
    double* Vm = optdata->Vm; 
    double* Z0 = optdata->Z0; 
    if (this->dim == 2)
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
    else
    {
      std::cerr << "Only dim=2 is supported. Exiting..\n";
      exit(1);
    }

  
    if (this->dim == 2)
    {
      // copy into Z0 for opt routines
      for (unsigned int i = 0; i < N; ++i)
      {
        Z0[i]       = optdata->X0[i];
        Z0[i + N]   = optdata->Y0[i];
      }
    }
    optF ( 2*optdata->N, Z0, nullptr, optdata );
  }
  


 
  inline void runXW ()
  {
    double minF = optF(2*optdata->N, optdata->Z0, nullptr, optdata); 
    if (minF > tolc)
    {
      std::cout <<"\n BEGIN NLOPT (X,W) \n";
      std::cout << "Initial value of objective : " << minF << std::endl; 
      unsigned int N = optdata->N; 
      nlopt_opt opt, opt_i;
      double *lb, *ub, *tolieq;
      if (this->dim == 2)
      {
        // x, y > 0
        lb = (double*) calloc(2*N, sizeof(double)); 
        // x, y < 1
        ub = (double*) calloc(2*N, sizeof(double));
        tolieq = (double*) calloc(N, sizeof(double));
        for (unsigned int i = 0; i < 2*N; ++i) 
        { 
          ub[i] = 1.0 - 1e-8;
          lb[i] = 1e-8; 
        }
        for (unsigned int i = 0; i < N; ++i) { tolieq[i] = tolc; }
        //opt = nlopt_create(NLOPT_AUGLAG, 2*N);
        //opt_i = nlopt_create(alg, 2*N);
        //nlopt_set_local_optimizer(opt, opt_i);
        opt = nlopt_create(alg, 2*N);
        nlopt_set_lower_bounds(opt, lb);
        nlopt_set_upper_bounds(opt, ub);
        nlopt_add_inequality_mconstraint(opt, N, optieqC, optdata, tolieq);
        nlopt_set_min_objective(opt, optF, optdata);
      }
      else
      {
        std::cerr << "Only dim=2 is supported. Exiting..\n";
        exit(1);
      }

      nlopt_set_xtol_rel(opt, tol);
      nlopt_set_ftol_rel(opt, tol);
      //nlopt_set_stopval(opt, tol);

      nlopt_result result = nlopt_optimize(opt, optdata->Z0, &minF);
      if (result < 0 && result != -4) 
      {
        std::cerr << "NLOPT failed! \n"; 
        std::cerr << nlopt_result_to_string(result);
        std::cerr << "\n Exiting ..\n";
        exit(1);
      }
      else
      {
        if (result == -4)
        {
          std::cout << nlopt_result_to_string(result) << std::endl;
        }
        std::cout << "NLOPT converged with residual " << minF << std::endl;
      }
      nlopt_destroy(opt);
      nlopt_destroy(opt_i);
      free(lb); free(ub);
      free(tolieq);
      optdata->minf = minF;
    } 
    std::cout << "Quadrature exactness in norm up to specified order : " 
              << minF << std::endl;

    double rcond[1];
    cond(optdata->Pm->V, optdata->N, optdata->M, rcond);
    std::cout << "Conditioning of interpolation operator on new abscissa : "
              << 1.0 / rcond[0] << std::endl; 
    std::cout <<"END NLOPT \n\n";
  }
  
  inline void runX ()
  {
    std::cout <<"\n BEGIN NLOPT cond(V(Z))\n"; 
    unsigned int N = optdata->N;
    nlopt_opt opt; 
    double *ub, *lb, *tolieq;
    if (this->dim == 2)
    { 
      // x, y > 0
      lb = (double*) calloc(2*N, sizeof(double)); 
      // x, y < 1
      ub = (double*) calloc(2*N, sizeof(double));
      tolieq = (double*) calloc(2*N, sizeof(double));
      for (unsigned int i = 0; i < 2*N; ++i) 
      { 
        ub[i] = 1.0 - 1e-8;
        lb[i] = 1e-8; 
      }
      for (unsigned int i = 0; i < 2*N; ++i) { tolieq[i] = tolc; }
      opt = nlopt_create(NLOPT_LN_COBYLA, 2*N); 
      nlopt_add_inequality_mconstraint(opt, 2*N, optieqC1, optdata, tolieq);
    }
    else
    {
      std::cerr << "Only dim=2 is supported. Exiting..\n";
      exit(1);
    }
  
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, optF1, optdata);
    nlopt_add_equality_constraint(opt, opteqC1, optdata, tolc);
  
    nlopt_set_xtol_rel(opt, 1e-1);
    nlopt_set_stopval(opt, 3);
    if (nlopt_optimize(opt, optdata->Z0, &optdata->minf1) < 0) 
    {
      std::cerr << "NLOPT failed!" << std::endl;
    }
    else 
    {
      std::cout << "Conditioning of interpolation operator on new abscissa : " 
                << optdata->minf1 << std::endl;
    }
    nlopt_destroy(opt);
    free(lb); free(ub);
    free(tolieq);
    std::cout <<"END NLOPT \n\n";
  }
};



#endif
