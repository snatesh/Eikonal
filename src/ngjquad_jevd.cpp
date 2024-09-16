#include<ngjquad_jevd.h>
#include<fstream>
#include<../include/matplotlibcpp.h>
namespace plt = matplotlibcpp;
void plot_tri(double* tri, std::string col);

/*

  Usage: ./ngjquad_opt algtype, n m tol tolc use_newton use_wolfe alph
    -  algtype (int) - 0-2 
        (0 - NELDERMEAD, 1 - COBYLA, 2 - SBPLX)
    -  n,m (int) > 0
    -  tol,tolc (double) > 0
    -  use_newton, use_wolfe (int) - 0,1
    -  alph (double) > 0
  
   n is max homogeneous poly degree in source basis 
   m is max homogeneous poly degree in target basis

   N = (n+1) choose (n-1) is number of nodes in quad rule
   and also num polys up to total degree n in source basis
   
   m is max homogeneous poly degree attempted to integrate
   exactly by quad rule of size N
   
   There are M polys up to total degree m in target basis 
   we seek quad rules of size N < M 

   This routine calls LAPACK complex eigenvalue routine
   ZGEEV to compute initial nodes for optimization
 
*/
int main(int argc, char* argv[])
{

  unsigned int nthreads;

  // parse command line args
  nlopt_algorithm alg; 
  int n, m; 
  double tol, tolc, alph;
  bool use_newton, use_wolfe; 
  set_args( argc, argv, alg, n, m, tol, tolc, 
            use_newton, use_wolfe, alph, nthreads);
  // initialize opt routines
  double a, b, c, kap; a = b = c = 0.5; kap = abs(a+b+c);
  double wabc = tgamma(kap+1.5) / ( tgamma(a+0.5) * tgamma(b+0.5) * tgamma(c+0.5) );
  unsigned int N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
  unsigned int M = static_cast<unsigned int>(0.5 * m * (m + 1));
  nlopt_func_data* d = new nlopt_func_data(m, n, a, b, c, nthreads);
  init_opt(d);
  //std::vector<double> X(d->Z0, d->Z0+N); std::vector<double> Y(d->Z0+N, d->Z0+2*N);
  //double tri[6] = {0, 1, 0, 0, 0, 1}; 
  //plot_tri(tri, "k-");
  //plt::plot(X, Y, "ro"); 
  //plt::show(); 
  // run nlopt
  nlopt_run ( alg, d, tol, tolc );
  // newton relaxation
  if (use_newton)
  {
    newton(d->F0, d->Vm, d->Hm, d->N, d->m, a, b, c, d->Z0, use_wolfe, alph);
  }
  // integration test on result d->Z0
  double rcond; double sum = 0;
  for (unsigned int i = 0; i < N; ++i) { sum += d->Z0[2*N+i]; }
  std::cout << "Sum of weights : " << sum << std::endl;  
  std::cout << "Objective value at argmin: " << nloptF(3*N, d->Z0, nullptr, d)  << std::endl;
  cond(d->Vm, N, M, &rcond);
  std::cout << "Conditioning of interpolation operator: " << 1.0 / rcond << std::endl;
  double* Ftest = (double*) calloc(N, sizeof(double));
  for (unsigned int i = 0; i < N; ++i) { Ftest[i] = std::sin( pow(d->Z0[i], 2) + pow(d->Z0[i+N], 2) ); }
  double Ival = cblas_ddot(N, d->Z0 + 2*N, 1, Ftest, 1) / wabc;
  printf("Integral_T sin(x^2+y^2) : %5.16f \n", Ival);
  for (unsigned int i = 0; i < N; ++i)
  {
    if (d->Z0[i] < 0 || d->Z0[i+N] < 0 || d->Z0[i] + d->Z0[i+N] - 1 > 0)
    {
      std::cerr << "ERROR: nodes outside of triangle!" << std::endl; 
      std::stringstream ss; ss << "triquadLeg_" << n << "_" << m << ".txt";
      std::ofstream ofile(ss.str());
      for (unsigned int i = 0; i < 3*N; ++i) 
      { 
        ofile << std::setprecision(std::numeric_limits<double>::max_digits10) 
              << d->Z0[i] << std::endl; 
      }
      delete d;
      free(Ftest);    
      exit(1); 
    }
  }
  
  // run nlopt on next objective
  nlopt_run1 (alg, d, tol, tolc ); 
  // newton relaxation
  if (use_newton)
  {
    newton(d->F0, d->Vm, d->Hm, d->N, d->m, a, b, c, d->Z0, use_wolfe, alph);
  }

  
  std::stringstream ss; ss << "triquadLeg_" << n << "_" << m << ".txt";
  std::ofstream ofile(ss.str());
  for (unsigned int i = 0; i < 3*N; ++i) 
  { 
    ofile << std::setprecision(std::numeric_limits<double>::max_digits10) 
          << d->Z0[i] << std::endl; 
  }


  // integration test on result d->Z0
  sum = 0.0;
  for (unsigned int i = 0; i < N; ++i) { sum += d->Z0[2*N+i]; }
  std::cout << "Sum of weights : " << sum << std::endl;  
  std::cout << "Conditioning of interpolation operator : " << nloptF1(2*N, d->Z0, nullptr, d)  << std::endl;
  Ival = cblas_ddot(N, d->Z0 + 2*N, 1, Ftest, 1) / wabc;
  printf("Integral_T sin(x^2+y^2) : %5.16f \n", Ival);
 
  // plot resulting quadrature nodes
  //std::vector<double> X(d->Z0, d->Z0+N); std::vector<double> Y(d->Z0+N, d->Z0+2*N);
  //std::vector<double> X_0(d->X0, d->X0+N); std::vector<double> Y_0(d->Y0, d->Y0+N);
  //double tri[6] = {0, 1, 0, 0, 0, 1}; 
  //plot_tri(tri, "k-");
  //plt::plot(X, Y, "ro"); 
  //plt::plot(X_0, Y_0, "ks"); 
  //plt::show(); 
  // cleanup
  delete d; 
  free(Ftest);
  return 0;
}

void plot_tri(double* tri, std::string col)
{
  std::vector<double> x1, y1, x2, y2, x3, y3;
  x1 = {tri[0], tri[1]}; y1 = {tri[3], tri[4]};
  x2 = {tri[1], tri[2]}; y2 = {tri[4], tri[5]};
  x3 = {tri[2], tri[0]}; y3 = {tri[5], tri[3]};
  plt::plot(x1, y1, col);
  plt::plot(x2, y2, col);
  plt::plot(x3, y3, col);
}
