#include<ngjquad_opt.h>
#include<../include/matplotlibcpp.h>

namespace plt = matplotlibcpp;
void plot_tri(double* tri, std::string col);

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
  unsigned int n, m; double tol, tolc;
  bool use_newton = true; bool use_wolfe = true;
  double alph; nlopt_algorithm alg;
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
      return 1;
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


  /**************************** BEGIN OPT INITIALIZATION *******************/
  double a, b, c, kap; a = b = c = 0.5; kap = abs(a+b+c);
  double wabc = tgamma(kap+1.5) / ( tgamma(a+0.5) * tgamma(b+0.5) * tgamma(c+0.5) );
  unsigned int N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
  unsigned int M = static_cast<unsigned int>(0.5 * m * (m + 1));
  // generate jacobi matrices for x,y 
  double* Hm = (double*) calloc(m*m, sizeof(double));
  double* Jn1 = (double*) calloc(N*N, sizeof(double));
  double* Jn2 = (double*) calloc(N*N, sizeof(double));
  Complex* Jn = (Complex*) calloc(N*N, sizeof(Complex));
  Complex* XY0 = (Complex*) calloc(N, sizeof(Complex));
  double* X0  = (double*) calloc(N, sizeof(double));
  double* Y0  = (double*) calloc(N, sizeof(double));
  double* Vm = (double*) calloc(N*M, sizeof(double));
  double* Vm_T = (double*) calloc(M*N, sizeof(double));
  double* W0 = (double*) calloc(M, sizeof(double)); W0[0] = 1.0;
  double* S = (double*) calloc(M, sizeof(double)); 
  double* Z0 = (double*) calloc(3*N, sizeof(double));
  double* F0 = (double*) calloc(3*N, sizeof(double));
  nlopt_func_data* d = (nlopt_func_data*) malloc(sizeof(nlopt_func_data));

  initialize(
              m, n, M, N, a, b, c,
              Jn1, Jn2, Jn, XY0,
              X0, Y0, Hm, Vm, Vm_T, 
              W0, S, Z0, F0, d 
            );
  
  /*************************** END INITIALIZATION **********************/


  /***************************** BEGIN NLOPT **************************/
  
  std::cout <<"\n BEGIN NLOPT \n";

  // x, y , w > 0
  double* lb = (double*) calloc(3*N, sizeof(double)); 
  // x, y , w < 1
  double* ub = (double*) calloc(3*N, sizeof(double));
  double* tolieq = (double*) calloc(2*N, sizeof(double));
  double* ineqres = (double*) calloc(2*N, sizeof(double));
  for (unsigned int i = 0; i < 3*N; ++i) { ub[i] = 1; }
  for (unsigned int i = 0; i < 2*N; ++i) { tolieq[i] = tolc; }

  nlopt_opt opt = nlopt_create(alg, 3*N); 
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

  std::cout <<"END NLOPT \n\n";
 
  /************************** END NLOPT  *********************************/
  double* Zk = Z0; double* Fk = F0; 
  /************************ BEGIN NEWTON *********************************/
  if (use_newton)
  {
    newton(F0, Vm, Hm, N, m, a, b, c, Z0, use_wolfe, alph);
  }
  /************************ END NEWTON *********************************/
  double rcond;
  double sum = 0;
  for (unsigned int i = 0; i < N; ++i) { sum += Zk[2*N+i]; }
  std::cout << "final sum of weights : " << sum << std::endl;  
  std::cout << "final value of objective : " << nloptF(3*N, Z0, nullptr, d)  << std::endl;
  jPoly_tri<double>(Zk, Zk + N, Hm, N, m-1, a, b, c, Vm);
  cond(Vm, N, M, &rcond);
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
  plot_tri(tri, "k-");
  plt::plot(X, Y, "ro");
  plt::plot(X_0, Y_0, "ks");
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

namespace plt = matplotlibcpp;

void plot_tri(double* tri, std::string col)
{
  std::vector<double> x1( {tri[0], tri[1]} );
  std::vector<double> y1( {tri[3], tri[4]} );
  std::vector<double> x2( {tri[1], tri[2]} );
  std::vector<double> y2( {tri[4], tri[5]} );
  std::vector<double> x3( {tri[2], tri[0]} );
  std::vector<double> y3( {tri[5], tri[3]} );
  plt::plot(x1, y1, col);
  plt::plot(x2, y2, col);
  plt::plot(x3, y3, col);
}
