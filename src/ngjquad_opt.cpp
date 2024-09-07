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




double F(double* Zk, unsigned int N, unsigned int m, double* Hm, double a, double b, double c)
{
  double *Xk, *Yk, *Wk;
  Xk = Zk; Yk = Zk + N; Wk = Zk + 2*N;
  unsigned int M = static_cast<unsigned int>(0.5 * m * (m + 1));
  double* Vm = (double*) calloc(N*M, sizeof(double));
  double* I0 = (double*) calloc(M, sizeof(double)); I0[0] = 1;  
  jPoly_tri<double>(Xk, Yk, Hm, N, m-1, a, b, c, Vm);
  cblas_dgemv(CblasColMajor,  CblasTrans,  N, M, 1.0, Vm, N, Wk, 1, -1.0, I0, 1);
  double norm = cblas_dnrm2(M, I0, 1);
  free(Vm);
  free(I0);
  return norm * norm;
}

int main(int argc, char* argv[])
{
  // jacobi poly params
  double a, b, c, kap; a = b = c = 0.5; kap = abs(a+b+c);
  // normalization for weight
  double wabc 
    = tgamma(kap+1.5) / 
      ( tgamma(a+0.5) * tgamma(b+0.5) * tgamma(c+0.5) );
  unsigned int n, m;
  
  /*
  - n is max poly degree in source basis 
  - N = (n+1) choose (n-1) is number of nodes in quad rule
  - m is max homogeneous poly degree attempted to integrate
    exactly by quad rule of size N
  - there are M polys up to total degree m
  */
  
  n = 4; m = 7;
  // normalizations for polys in target basis
  double* Hm = (double*) calloc(m*m, sizeof(double));
  structure_factors_tri<double>(m, a, b, c, Hm);
  // number of polynomials in source basis
  unsigned int N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
  // number of polynomials in target basis
  unsigned int M = static_cast<unsigned int>(0.5 * m * (m + 1));
  // generate jacobi matrices for x,y 
  double* Jn1 = (double*) calloc(N*N, sizeof(double));
  double* Jn2 = (double*) calloc(N*N, sizeof(double));
  Complex* Jn = (Complex*) calloc(N*N, sizeof(Complex));
  jacobi_mat_ON_tri<double>(n, a, b, c, Jn1, Jn2);
  // compute initial nodes and weights from eigenvalues of Jn
  for (unsigned int i = 0; i < N*N; ++i) { Jn[i] = Jn1[i] + I*Jn2[i]; }
  Complex* X0 = (Complex*) calloc(N, sizeof(Complex));
  double* Xk  = (double*) calloc(N, sizeof(double));
  double* Yk  = (double*) calloc(N, sizeof(double));
  if (LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'N', N, Jn, N, X0, NULL, N, NULL, N))
  {
    std::cerr << "ERROR: Lapack ZGEEV: Eigenvalues" << std::endl;
  }
  for (unsigned int i = 0; i < N; ++i) 
  { 
    Xk[i] = creal(X0[i]); 
    Yk[i] = cimag(X0[i]); 
  } 
  // evaluate Vandermonde on initial nodes
  double* Vm = (double*) calloc(N*M, sizeof(double));
  double* Vm_T = (double*) calloc(M*N, sizeof(double));
  jPoly_tri<double>(Xk, Yk, Hm, N, m-1, a, b, c, Vm);
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

  double* Wk = (double*) calloc(N, sizeof(double)); 
  double* I0 = (double*) calloc(M, sizeof(double)); I0[0] = 1.0;
  double* S = (double*) calloc(N, sizeof(double));
  lapack_int rank[1];
  if (LAPACKE_dgelsd(LAPACK_COL_MAJOR, M, N, 1, Vm_T, M, I0, M, S, -1.0, rank))
  {
    std::cerr << "ERROR: Lapack dgelsd: Pseudoinverse" << std::endl;
  }

  std::vector<double> X(Xk, Xk+N);
  std::vector<double> Y(Yk, Yk+N);

  double tri[6] = {0, 1, 0, 0, 0, 1};
  plot_tri(tri);
  plt::plot(X, Y, "r.");
  plt::show(); 
 
  double sum = 0;
  for (unsigned int i = 0; i < N; ++i)
  {
    Wk[i] = I0[i];
    sum += Wk[i];

  } 

  double* Zk = (double*) calloc(3*N, sizeof(double));
  for (unsigned int i = 0; i < N; ++i)
  {
    Zk[i] = Xk[i];
    Zk[i + N] = Yk[i];
    Zk[i + 2*N] = Wk[i];
  }



  std::cout << "sum of initial weights : " << sum << std::endl;  
  std::cout << "initial value of objective : " <<  F(Zk, N, m, Hm, a, b, c) << std::endl;

  
  free(Jn1);
  free(Jn2);
  free(Jn);
  free(X0);
  free(Xk);
  free(Yk);
  free(Hm);
  free(Vm);
  free(Vm_T);
  free(Wk);
  return 0;

  //std::stringstream ss;
  //ss << "triquadLeg_" << n << "_" << m << ".txt";
  //std::string fname(ss.str()); ss = "";
  

}

