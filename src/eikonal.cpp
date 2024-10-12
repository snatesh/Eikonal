#include<eikonal.hh>
#include<cmath>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<map>
#include<random>

using std::sin;
using std::cos;
using std::pow;
double Frhs(double x, double y)
{

//return y*y*pow((7*pow(x,6) + 6*pow(x,5)*(-1 + y) + 2*x*pow(y,7) + (-1 + y)*pow(y,7)),2) + 
//  x*x*pow((pow(x,6) + 8*x*pow(y,7) + pow(x,5)*(-1 + 2*y) + pow(y,7)*(-8 + 9*y)),2); 
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
  //return x*(1-x-y)*y;//*(pow(x,5) + pow(y,7));
  double z = dtoh(x, y);
  double min = (x < y ? x : y);
  double min1 = (min < z ? min : z);
  return min1;
}



double L2err(eikonal* solver, double* fu)
{
  double* U = (double*) calloc(solver->N, sizeof(double));
  double* diff = (double*) calloc(solver->N, sizeof(double));
  solver->btransform(U);
  for (unsigned int i = 0; i < solver->N; ++i)
  {
    diff[i] = std::abs(U[i]-fu[i]);
  }
  
  double Ifu = cblas_ddot(solver->N, solver->W, 1, fu, 1);
  std::cout << Ifu << std::endl;
  std::cout << cblas_ddot(solver->N, solver->W, 1, U, 1);
  double abserr = cblas_ddot(solver->N, solver->W, 1, diff, 1);
  free(U); 
  free(diff);
  return abserr / std::abs(Ifu);
}

int main(int argc, char* argv[])
{
  bool unconstrained = false; 

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
  
  unsigned int n = std::stoi(argv[1]); 
  eikonal* solver = new eikonal ( n, N, 2*n, 2*n, 2*n, frhs, fu, 
                                  X, Y, W, nthreads, false, unconstrained);

  unsigned int Nopt = solver->Nopt;
  double *lob, *upb, *toleq; 
  lob   = (double*) malloc(Nopt * sizeof(double));
  upb   = (double*) malloc(Nopt * sizeof(double));
  toleq = (double*) malloc(solver->Ne * sizeof(double)); 
  double tol = 1e-13;
  for (unsigned int i = 0; i < Nopt; ++i)
  {
    upb[i]     = HUGE_VAL; 
    lob[i]     = -HUGE_VAL; 
  }
  for (unsigned int i = 0; i < solver->Ne; ++i) { toleq[i] = tol; }

  double* cu_sol;
  if (unconstrained)  { cu_sol = solver->cu_eq; }
  else                { cu_sol = solver->cu; }
  
  optDataEik* data = new optDataEik ( solver->Nopt, solver->Dx, solver->Dy, 
                                      solver->G, solver->crhs, cu_sol, 
                                      solver->Ne, solver->polyedge->V); 

  nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, solver->Nopt); 
  nlopt_set_lower_bounds(opt, lob);
  nlopt_set_upper_bounds(opt, upb);
  nlopt_set_min_objective(opt, F, data);
  if (not unconstrained)
  {
    nlopt_add_inequality_mconstraint(opt, solver->Ne, cl, data, toleq);
  }
  nlopt_set_xtol_rel(opt, sqrt(tol));
  nlopt_set_ftol_rel(opt, sqrt(tol));
  nlopt_set_stopval(opt, tol);
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
  free(toleq);
  free(Z);
  free(frhs);
  free(fu);
  free(U);
}
