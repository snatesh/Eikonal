#include<ngjquad.hh>
#include<legQuad.hh>
#include<matplotlibcpp.h>
#include<map>
#include<vector>
namespace plt = matplotlibcpp;
void plot_tri(double* tri, std::string col);



/*

  Usage: ./ngjquad algtype, n m tol tolc use_newton use_wolfe alph
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
          alg = NLOPT_LD_SLSQP;
          std::cout << "ALGORITHM: SLSQP\n";
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
      alg = NLOPT_LD_SLSQP;
      std::cout << "ALGORITHM: SLSQP\n";
      dim = 2; n = 4; m = 6;
      tol = 1e-13;
      tolc = 1e-13;
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
  ngjQuad* gjquad;
  legQuad<double>* legq;
  if (dim == 1)
  {
    gjquad = new ngjQuad  ( n, m, 0, 0, tol, tolc, alph,
                            use_newton, use_wolfe, alg, nthreads );
  }
  if (dim == 2)
  {
    gjquad = new ngjQuad  ( n, m, a, b, c, tol, tolc, 
                            alph, use_newton, use_wolfe, alg,
                            nthreads  );
  }
  else if (dim == 3)
  {
    unsigned int nlg = 20;
    legq = new legQuad<double>(nlg); 
    gjquad = new ngjQuad( n, m, a, b, c, d, legq->x, legq->w,
                          tol, tolc, alph, use_newton, use_wolfe, 
                          alg, nthreads );
  }

  gjquad->init();
  double* Z0 = gjquad->optdata->Z0;
  unsigned int N = gjquad->optdata->N;
  double sumw = 0;
  for (unsigned int i = dim*N; i < (dim+1)*N; ++i) { sumw += Z0[i]; } 
  std::cout << "(initial) Sum of weights : " << sumw << std::endl;
  
  if (use_newton) gjquad->newton();
  
  gjquad->runXW();
  
  //if (dim > 1) { gjquad->runX(); }

  sumw = 0;
  for (unsigned int i = dim*N; i < (dim+1)*N; ++i) { sumw += Z0[i]; } 
  std::cout << "(final ) Sum of weights : " << sumw << std::endl;
  
  if (dim == 3)
  {
    std::stringstream ss;
    ss  << "_N" << gjquad->optdata->N << "_n" << gjquad->optdata->n-1 
        << "_M" << gjquad->optdata->N << "_m" << gjquad->optdata->m-1;
    std::stringstream ssx, ssy, ssz, ssw;
    ssx << "xtet" << ss.str() << ".txt";
    ssy << "ytet" << ss.str() << ".txt";
    ssz << "ztet" << ss.str() << ".txt";
    ssw << "wtet" << ss.str() << ".txt";
    std::ofstream xfile(ssx.str());
    std::ofstream yfile(ssy.str());
    std::ofstream zfile(ssz.str());
    std::ofstream wfile(ssw.str());
    for (unsigned int i = 0; i < gjquad->optdata->N; ++i)
    {
      xfile << gjquad->optdata->Z0[i] << std::endl; 
      yfile << gjquad->optdata->Z0[i+N] << std::endl; 
      zfile << gjquad->optdata->Z0[i+2*N] << std::endl; 
      wfile << gjquad->optdata->Z0[i+3*N] << std::endl; 
    }
    xfile.close();
    yfile.close();
    zfile.close();
    wfile.close();
    // integrate test function on the tet
    double* Ftest = (double*) calloc(N, sizeof(double));
    for (unsigned int i = 0; i < N; ++i)
    {
      Ftest[i] = std::sin(  Z0[i]*Z0[i] + Z0[i+N]*Z0[i+N] + 
                            std::pow(2.0, Z0[i+2*N])  );
    }
    double Ival = cblas_ddot(N, Z0 + 3*N, 1, Ftest, 1);
    free(Ftest); 
    printf("Integral : %5.16f \n", Ival / 6.0);
  }  

  if (dim == 1)
  {
      // integrate test function on line
      double* Ftest = (double*) calloc(N, sizeof(double));
      for (unsigned int i = 0; i < N; ++i) 
      {
        Z0[i]  = (Z0[i]+1)/2.0;
        Z0[i+N] = 2 * Z0[i+N] / 2.0; 
        Ftest[i] = std::exp(std::sin(Z0[i]*Z0[i]) + std::pow(2.0, Z0[i]) - 5*Z0[i]);
      }
      double Ival = cblas_ddot(N, Z0 + N, 1, Ftest, 1);
      free(Ftest);

      printf("Integral : %5.16f \n", Ival);
      std::vector<double> x, y, z;
      x = {0, 1, 0, 0, 1};
      y = {0, 0, 1, 0, 0};
      z = {0, 0, 0, 1, 0};
      plt::plot3(x, y, z);//, kwargs);
      plt::xlim(0.0,1.0);
      plt::ylim(0.0,1.0);
      plt::show();
  }

  if (dim == 2)
  { 
    // integrate test function on tri
    double* Ftest = (double*) calloc(N, sizeof(double));
    for (unsigned int i = 0; i < N; ++i) 
    { 
      Ftest[i] = std::sin( pow(Z0[i], 2) + pow(Z0[i+N], 2) ); 
    }
    double Ival = cblas_ddot(N, Z0 + 2*N, 1, Ftest, 1) / wabc;
    free(Ftest);

    printf("Integral : %5.16f \n", Ival);
    // plot resulting quadrature nodes
    std::vector<double> X(Z0, Z0+N); std::vector<double> Y(Z0+N, Z0+2*N);
    double tri[6] = {0, 1, 0, 0, 0, 1}; 
    plot_tri(tri, "k-");
    plt::plot(X, Y, "ro"); 
    plt::show(); 
  }
  delete gjquad;
  if (dim == 3) { delete legq; }
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

void plot_tet()
{
  std::map<std::string, std::string> kwargs;
  kwargs["marker"] = "o";
  kwargs["linestyle"] = "-";
  kwargs["linewidth"] = "1";
  kwargs["markersize"] = "12";
  std::vector<double> x, y, z;
  x = {0, 1, 0, 0, 1};
  y = {0, 0, 1, 0, 0};
  z = {0, 0, 0, 1, 0};
  plt::plot3(x, y, z);//, kwargs);
  plt::xlim(0.0,1.0);
  plt::ylim(0.0,1.0);
  plt::show();
}



