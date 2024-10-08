#include<eikonal.hh>
#include<cmath>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<matplotlibcpp.h>
#include<vector>
#include<map>

namespace plt = matplotlibcpp;

double Frhs(double x, double y)
{
  return 1.0;
}

double Fu(double x, double y)
{
  return -1.0 * x * y * (x * x + y * y);
}


int main(int argc, char* argv[])
{
 
  std::ifstream Zfile("../testing/testdata/triquadLeg_16_27.txt");
  
  unsigned int N = 136;
  unsigned int nthreads = 6;
  double* Z = (double*) calloc(3*N, sizeof(double));
  double* frhs = (double*) calloc(N, sizeof(double));
  double* fu = (double*) calloc(N, sizeof(double));
  for (unsigned int i = 0; i < 3*N; ++i) { Zfile >> Z[i]; } 
  double* X = Z;
  double* Y = &Z[N];
  double* W = &Z[2*N];

  for (unsigned int i = 0; i < N; ++i)
  {
    frhs[i] = Frhs(X[i], Y[i]);
    fu[i] = Fu(X[i], Y[i]);
  }

  unsigned int n = 15; 
 
  eikonal* solver = new eikonal ( n, N, N, N, N, frhs, fu, 
                                  X, Y, W, nthreads );
  printMat(solver->Dx, solver->Np, solver->Np);
  printMat(solver->Dy, solver->Np, solver->Np);
  printMat(solver->G, solver->Np, solver->Np);
  printMat(solver->cu, solver->Np, 1);
  printMat(solver->crhs, solver->Np, 1);
  
  unsigned int Np = solver->Np;
  double *lob, *upb, *toleq; 
  lob   = (double*) malloc(Np * sizeof(double));
  upb   = (double*) malloc(Np * sizeof(double));
  toleq = (double*) malloc(solver->Ne * sizeof(double)); 
  double tol = 1e-14;
  for (unsigned int i = 0; i < Np; ++i)
  {
    upb[i]     = 1000 ; 
    lob[i]     = -1000; 
  }
  for (unsigned int i = 0; i < 3*N; ++i) { toleq[i] = tol; }
  
  optDataEik* data = new optDataEik ( solver->Np, solver->Dx, solver->Dy, 
                                      solver->G, solver->crhs, solver->cu, 
                                      solver->Ne, solver->polyedge->V); 
  nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, Np); 
  //nlopt_opt local_opt = nlopt_create(NLOPT_LD_CCSAQ, N);
  nlopt_set_lower_bounds(opt, lob);
  nlopt_set_upper_bounds(opt, upb);
  nlopt_set_min_objective(opt, F, data);
  nlopt_add_equality_mconstraint(opt, solver->Ne, cl, data, toleq);
  nlopt_set_xtol_rel(opt, tol);
  nlopt_set_ftol_rel(opt, tol);
  //nlopt_set_xtol_rel(local_opt, tol);
  //nlopt_set_ftol_rel(local_opt, tol);
  //nlopt_set_local_optimizer(opt, local_opt);
  double minF;
  nlopt_result result = nlopt_optimize(opt, data->cu, &minF);
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
  counteik = 0;
  printMat(solver->cu, solver->Np, 1);
  

  double* U = (double*) calloc(N, sizeof(double));
  solver->btransform(U);
  std::ofstream ofile("Usol.txt");
  for (unsigned int i = 0; i < N; ++i)
  {
    ofile << U[i] << std::endl;
  }
  ofile.close();
  std::vector<double> x(X, X+N);
  std::vector<double> y(Y, Y+N);
  std::vector<double> u(U, U+N); 

  std::map<std::string, std::string> kwargs;
  kwargs["marker"] = "o";
  kwargs["linestyle"] = "-";
  kwargs["linewidth"] = "1";
  kwargs["markersize"] = "12";

 //const long fg = plt::figure();
  plt::figure_size(1200, 780);
  plt::plot3(x, y, u);
  plt::xlim(0.0,1.0);
  plt::ylim(0.0,1.0);
  plt::show();

  delete data; 
  delete solver;
  nlopt_destroy(opt);
  free(lob);
  free(upb);
  free(Z);
  free(U);
  //nlopt_destroy(local_opt);
  
}
