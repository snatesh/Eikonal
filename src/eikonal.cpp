#include<cblas.h>
#include<jPoly.h>
#include<ngjquad.h>
#include<kMat.h>
#include<dMat.h>

inline void set_args  ( int argc, char* argv[], 
                        nlopt_algorithm& alg, int& dim, int& n, 
                        int& m, double& tol, double& tolc, 
                        bool& use_newton, bool& use_wolfe, 
                        double& alph, unsigned int& nthreads,
                        std::string& dir) 
{
  std::cout << "\nSEARCH FOR NEAR OPTIMAL GAUSSIAN QUADRATURE\n";
  switch(argc)
  {
    case 11:
    {
      nthreads = std::stoi(argv[10]);
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
      dim = std::stoi(argv[2]);
      n = std::stoi(argv[3]);
      m = std::stoi(argv[4]);
      tol = std::stod(argv[5]); 
      tolc = std::stod(argv[6]);
      use_newton = (bool) std::stoi(argv[7]); 
      use_wolfe = (bool) std::stoi(argv[8]);
      alph = std::stod(argv[9]);
      if (dim == 3 && dir == "") { dir = "../testing/testdata/"; }
      break;
    }
    case 1:
    {
      nthreads = 1;
      alg = NLOPT_LN_NELDERMEAD;
      std::cout << "ALGORITHM: NELDERMEAD\n";
      dim = 2; n = 4; m = 6;
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
        "Usage: ./ngjquad_opt algtype dim n m tol tolc use_newton use_wolfe alph\n";
      std::cerr << "  algtype (int) - 0-2\n";
      std::cerr << "  dim (2,3)\n";
      std::cerr << "  n,m (int) > 0\n";
      std::cerr << "  tol,tolc (double) > 0\n";
      std::cerr << "  use_newton, use_wolfe (int) - 0,1\n";
      std::cerr << "  alph (double) > 0\n";
      std::cerr << "  nthreads (int) > 0\n";
      exit(1);
    }
  }
  std::cout << "\nPARAMETERS:" << std::endl;
  std::cout << std::setw(15) << "n          = " << n << std::endl;
  std::cout << std::setw(15) << "m          = " << m << std::endl;
  std::cout << std::setw(15) << "tol        = " << tol << std::endl;
  std::cout << std::setw(15) << "tolc       = " << tolc << std::endl;
  std::cout << std::setw(15) << "use_newton = " << 
                                (use_newton ? "true\n" : "false\n");
  std::cout << std::setw(15) << "use_wolfe  = " << 
                                (use_newton && use_wolfe ? "true\n" : "false\n");
  std::cout << std::setw(15) << "alpha      = " << alph << std::endl;
  std::cout << std::setw(15) << "nthreads   = " << nthreads << std::endl;


}


int main(int argc, char* argv[])
{

  unsigned int nthreads;

  // parse command line args
  nlopt_algorithm alg; 
  int n, m; 
  double tol, tolc, alph;
  std::string dir;
  int dim;
  bool use_newton, use_wolfe; 
  set_args( argc, argv, alg, dim, n, m, tol, tolc, 
            use_newton, use_wolfe, alph, nthreads, dir);
  // initialize opt routines
  double a, b, c, d, kap; a = b = c = d = 0.5; kap = abs(a+b+c);
  double wabc = tgamma(kap+1.5) / ( tgamma(a+0.5) * tgamma(b+0.5) * tgamma(c+0.5) );
  double wabcd = 6;
  ngjQuad* gjquad = new ngjQuad ( n, m, a, b, c, tol, tolc, 
                                  alph, use_newton, use_wolfe, alg,
                                  nthreads  );

  gjquad->init();
  gjquad->runXW();
  double* Z0 = gjquad->optdata->Z0;
  unsigned int n_k = n - 1; 
  unsigned int N = rn3(n_k);

  // (a,b,c)->(a+1,b,c)->(a+1,b+1,c)->(a+1,b+1,c+1)
  double** Ks = (double**) calloc(3, sizeof(double*));
  double** H = (double**) calloc(6, sizeof(double*));
  double** D = (double**) calloc(2, sizeof(double*));
  double* K0 = (double*) calloc(N*N, sizeof(double));
  double* K = (double*) calloc(N*N, sizeof(double));
  double* Dx = (double*) calloc(N*N, sizeof(double));
  double* Dy = (double*) calloc(N*N, sizeof(double));
  for (unsigned int i = 0; i < 6; ++i)
  {
    H[i] = (double*) calloc((n+1)*(n+1), sizeof(double)); 
  }
  sFactors(n+1, a, b, c, H[0]);
  sFactors(n+1, a+1, b, c, H[1]);
  sFactors(n+1, a+1, b+1, c, H[2]);
  sFactors(n+1, a+1, b+1, c+1, H[3]);
  sFactors(n+1, a+1, b, c+1, H[4]);
  sFactors(n+1, a, b+1, c+1, H[5]);
  for (unsigned int i = 0; i < 3; ++i)
  {
    Ks[i] = (double*) calloc(N*N, sizeof(double));
  }
  kMat(a, b, c, H[0], H[1], n_k, 0, Ks[0]);
  kMat(a+1, b, c, H[1], H[2], n_k, 1, Ks[1]);
  kMat(a+1,b+1, c, H[2], H[3], n_k, 2, Ks[2]);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
              N, N, N, 1.0, Ks[1], N, Ks[0], N, 0.0, K0, N);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
              N, N, N, 1.0, Ks[2], N, K0, N, 0.0, K, N);
  for (unsigned int i = 0; i < 2; ++i)
  {
    D[i] = (double*) calloc(N*N, sizeof(double));
  }
  dMat(a, b, c, H[0], H[4], n_k, 0, D[0]);
  dMat(a, b, c, H[0], H[5], n_k, 1, D[1]);
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
              N, N, N, 1.0, K, N, D[0], N, 0.0, Dx, N); 
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
              N, N, N, 1.0, K, N, D[1], N, 0.0, Dy, N); 
  
  printMat(Dx, N, N);
  printMat(Dy, N, N);

  for (unsigned int i = 0; i < 6; ++i)
  {
    free(H[i]);
  }
  free(H);
  for (unsigned int i = 0; i < 3; ++i)
  {
    free(Ks[i]);
  }
  free(Ks); free(K); free(K0);
  for (unsigned int i = 0; i < 2; ++i)
  {
    free(D[i]);
  }
  free(D); free(Dx); free(Dy);
  return 0;
}
