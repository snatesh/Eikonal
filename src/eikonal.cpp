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

double dtoh(double x, double y)
{
  double X[2] = {x, y};
  double Z[2] = {(x-y+1.0)/2.0, (y-x+1.0)/2.0};
  double d = std::sqrt( std::pow(X[0]-Z[0], 2.0) +
                        std::pow(X[1]-Z[1], 2.0)  );
  return d;
}


double Fu(double x, double y)
{
  double z = dtoh(x, y);
  double min = (x < y ? x : y);
  double min1 = (min < z ? min : z);
  return min1;
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
    std::cout << X[i] << " " << Y[i] << std::endl;
    frhs[i] = Frhs(X[i], Y[i]);
    fu[i] = Fu(X[i], Y[i]);
  }

  unsigned int n = 15; 
 
  eikonal* solver = new eikonal ( n, N, N, N, N, frhs, fu, 
                                  X, Y, W, nthreads );


  
  unsigned int Np = solver->Np;
  double *lob, *upb, *toleq; 
  lob   = (double*) malloc(Np * sizeof(double));
  upb   = (double*) malloc(Np * sizeof(double));
  toleq = (double*) malloc(solver->Ne * sizeof(double)); 
  double tol = 1e-13;
  for (unsigned int i = 0; i < Np; ++i)
  {
    upb[i]     = HUGE_VAL/2.0 ; 
    lob[i]     = -HUGE_VAL/3.0; 
  }
  for (unsigned int i = 0; i < solver->Ne; ++i) { toleq[i] = tol; }
  
  optDataEik* data = new optDataEik ( solver->Np, solver->Dx, solver->Dy, 
                                      solver->G, solver->crhs, solver->cu, 
                                      solver->Ne, solver->polyedge->V); 
  nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, Np); 
  //nlopt_opt local_opt = nlopt_create(NLOPT_LD_CCSAQ, N);
  nlopt_set_lower_bounds(opt, lob);
  nlopt_set_upper_bounds(opt, upb);
  nlopt_set_min_objective(opt, F, data);
  nlopt_add_equality_mconstraint(opt, solver->Ne, cl, data, toleq);
  nlopt_set_xtol_rel(opt, 1e-7);
  nlopt_set_ftol_rel(opt, 1e-7);
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

  
  delete data; 
  delete solver;
  nlopt_destroy(opt);
  free(lob);
  free(upb);
  free(Z);
  free(U);
  //nlopt_destroy(local_opt);
  
}
