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
#include<../include/matplotlibcpp.h>



typedef double _Complex Complex;
namespace plt = matplotlibcpp;

void plot_tri(double* tri)
{
  plt::plot((std::vector<double>) {tri[0],tri[1]}, (std::vector<double>) {tri[3],tri[4]},"b-");
  plt::plot((std::vector<double>) {tri[1],tri[2]}, (std::vector<double>) {tri[4],tri[5]},"b-");
  plt::plot((std::vector<double>) {tri[2],tri[0]}, (std::vector<double>) {tri[5],tri[3]},"b-");
}

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



void newton(double* Fk, double* Vm, double* Hm, unsigned int N, unsigned int m, double a, double b, double c, double* Zk)
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

/*  
  - n is max poly degree in source basis 
  - N = (n+1) choose (n-1) is number of nodes in quad rule
    and also num polys in source basis
  - there are N polys up to total degree n in source basis
  - m is max homogeneous poly degree attempted to integrate
    exactly by quad rule of size N
  - there are M polys up to total degree m in target basis 
  - we seek quad rules of size N < M 
*/

int main(int argc, char* argv[])
{
  unsigned int n, m;
  double tol, tolc;
  if (argc == 5)
  {
    n = std::stoi(argv[1]);
    m = std::stoi(argv[2]);
    tol = std::stod(argv[3]); 
    tolc = std::stod(argv[4]); 
  }
  else if (argc == 1)
  {
    n = 4; m = 6;
    tol = 1e-5; double tolc = 1e-5;
  }
  else
  {
    std::cerr << "Incorrect number of command line arguments (n, m, tol, tolc)\n";
    std::cerr << "Only " << argc << " provided" << std::endl;
  }
  double a, b, c, kap; a = b = c = 0.5; kap = abs(a+b+c);
  double wabc = tgamma(kap+1.5) / ( tgamma(a+0.5) * tgamma(b+0.5) * tgamma(c+0.5) );
  double* Hm = (double*) calloc(m*m, sizeof(double));
  structure_factors_tri<double>(m, a, b, c, Hm);
  unsigned int N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
  unsigned int M = static_cast<unsigned int>(0.5 * m * (m + 1));


  /**************************** BEGIN OPT INITIALIZATION *******************/
  // generate jacobi matrices for x,y 
  double* Jn1 = (double*) calloc(N*N, sizeof(double));
  double* Jn2 = (double*) calloc(N*N, sizeof(double));
  Complex* Jn = (Complex*) calloc(N*N, sizeof(Complex));
  jacobi_mat_ON_tri<double>(n, a, b, c, Jn1, Jn2);
  // compute initial nodes and weights from eigenvalues of Jn
  for (unsigned int i = 0; i < N*N; ++i) { Jn[i] = Jn1[i] + I*Jn2[i]; }
  Complex* XY0 = (Complex*) calloc(N, sizeof(Complex));
  double* X0  = (double*) calloc(N, sizeof(double));
  double* Y0  = (double*) calloc(N, sizeof(double));
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
  double* Vm = (double*) calloc(N*M, sizeof(double));
  double* Vm_T = (double*) calloc(M*N, sizeof(double));
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

  // solve least squares system for initial weights
  double* W0 = (double*) calloc(M, sizeof(double)); W0[0] = 1.0;
  double* S = (double*) calloc(N, sizeof(double)); lapack_int rank[1];
  if (LAPACKE_dgelsd(LAPACK_COL_MAJOR, M, N, 1, Vm_T, M, W0, M, S, -1.0, rank))
  {
    std::cerr << "ERROR: Lapack dgelsd: Pseudoinverse" << std::endl;
  }
  double* Z0 = (double*) calloc(3*N, sizeof(double));
  for (unsigned int i = 0; i < N; ++i)
  {
    Z0[i]       = X0[i];
    Z0[i + N]   = Y0[i];
    Z0[i + 2*N] = W0[i];
  }
  double* F0 = (double*) calloc(3*N, sizeof(double));
  F(Z0, Vm, N, m, Hm, a, b, c, F0);
  nlopt_func_data* d = (nlopt_func_data*) malloc(sizeof(nlopt_func_data));
  d->Vm = Vm; d->Hm = Hm; d->Fk = F0;
  d->a = a; d->b = b; d->c = c;
  d->N = N; d->m = m; d->M = M;
  
  /*************************** END INITIALIZATION **********************/


  /***************************** BEGIN NLOPT **************************/
  
  std::cout <<"\n\n BEGIN NLOPT \n\n";

  // x, y , w > 0
  double* lb = (double*) calloc(3*N, sizeof(double)); 
  // x, y , w < 1
  double* ub = (double*) calloc(3*N, sizeof(double));
  double* tolieq = (double*) calloc(2*N, sizeof(double));
  double* ineqres = (double*) calloc(2*N, sizeof(double));
  for (unsigned int i = 0; i < 3*N; ++i) { ub[i] = 1; }
  for (unsigned int i = 0; i < 2*N; ++i) { tolieq[i] = tolc; }

  nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, 3*N); 
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);
  nlopt_set_min_objective(opt, nloptF, d);
  nlopt_add_equality_constraint(opt, nlopteqC, NULL, tolc);
  nlopt_add_inequality_mconstraint(opt, 2*N, nloptieqC, NULL, tolieq);
  nlopt_set_stopval(opt, tol);
  double minF;
  if (nlopt_optimize(opt, Z0, &minF) < 0) 
  {
    std::cerr << "NLOPT failed!" << std::endl;
  }
  else 
  {
    std::cout << "found minimum with objective val = " << minF << std::endl;
  }
  nlopt_destroy(opt);

  std::cout <<"\n END NLOPT \n\n";
 
  /************************** END NLOPT  *********************************/

  /************************ BEGIN NEWTON *********************************/
  newton(F0, Vm, Hm, N, m, a, b, c, Z0);
  double* Zk = Z0; double* Fk = F0; 
  double sum = 0;
  for (unsigned int i = 0; i < N; ++i) { sum += Zk[2*N+i]; }
  std::cout << "final sum of weights : " << sum << std::endl;  
  std::cout << "final value of objective : " << nloptF(3*N, Z0, nullptr, d)  << std::endl;
  /************************ END NEWTON *********************************/

  
  /******************************** BEGIN INTEGRATION TEST **************/
  double* Ftest = (double*) calloc(N, sizeof(double));
  for (unsigned int i = 0; i < N; ++i) { Ftest[i] = std::sin( pow(Zk[i], 2) + pow(Zk[i+N], 2) ); }
  double Ival = cblas_ddot(N, Zk + 2*N, 1, Ftest, 1) / wabc;
  std::cout << "Integral val : "; 
  std::cout << std::setw(10) << Ival << std::endl;
  /******************* END INTEGRATION TEST ************************/


  /********************* BEGIN PLOT ****************************/
  std::vector<double> X(Zk, Zk+N);
  std::vector<double> Y(Zk+N, Zk+2*N);
  std::vector<double> X_0(X0, X0+N);
  std::vector<double> Y_0(Y0, Y0+N);
  
  double tri[6] = {0, 1, 0, 0, 0, 1};
  plot_tri(tri);
  plt::plot(X, Y, "ro");
  plt::plot(X_0, Y_0, "b.");
  plt::show(); 
  /********************* END PLOT ****************************/

  
  free(Jn1); free(Jn2); free(Jn);
  free(XY0); free(X0); free(Y0);
  free(Hm); free(Vm); free(Vm_T);
  free(W0); free(F0); free(Z0);
  free(S); free(d); free(ub); 
  free(lb); free(tolieq);
  return 0;
}


    // TODO: add below Wolf condition-based backtracking line search to newton()
    //for (unsigned int iter_i = 0; iter_i < maxiter; ++iter_i)
    //{
    //  cblas_dgemv(CblasColMajor,  CblasNoTrans, M, 3*N, gam*alph, gradFk, M, dZk, 1, 1.0, Fk, 1);
    //  wfnrm = cblas_dnrm2(M, Fk, 1);
    //  std::cout << "wfnrm " << wfnrm << std::endl;
    //  if (pk > wfnrm)
    //  {
    //    alph = rho*alph;
    //    for (unsigned int i = 0; i < 3*N; ++i) { Zk[i] -= alph*dZk[i]; }
    //    F(Zk, N, m, Hm, a, b, c, Fk);
    //    pk = cblas_dnrm2(M, Fk, 1);
    //  } 
    //}
    
