#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#include<cblas.h>
#include<gtest/gtest.h>
#include<structure_factors.h>
#include<jPoly.h>
#include<promotion_mat_tri.h>
#include<dxdy_mat_tri.h>
#include<jacobi_mat_ON_tri.h>

#ifdef PLOT
#include<../include/matplotlibcpp.h>
#endif

#include<nlopt.h>
#include<stdio.h>

/* 
  Top of file is reserved for test helper functions 
  and small-sized data definition
  
  gtest routines and main are at bottom of file 
*/

void printMat(const double* A, const unsigned int m, const unsigned int n)
{
  for (unsigned int i = 0; i < m; ++i)
  {
    for (unsigned int j = 0; j < n; ++j)
    {
      std::cout << std::setw(10);
      std::cout << A[i + m*j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}


double twonorm(const double* A, const double* B, const unsigned int N)
{
  double norm = 0.0;
  double absdiff;
  for (unsigned int i = 0; i < N; ++i) 
  {
    absdiff = std::abs(A[i]-B[i]);
    norm += pow(absdiff,2);
  }
  return sqrt(norm);
}

double vecdot(const double* A, const double* B, const unsigned int N)
{
  return cblas_ddot(N, A, 1, B, 1);
}

double infnorm(const double* A, const double* Aref, const unsigned int N)
{
  double norm = -1.0;
  double maxref = Aref[0];
  double absdiff;
  for (unsigned int i = 0; i < N; ++i) 
  {
    absdiff = std::abs(A[i]-Aref[i]);
    maxref = (std::abs(Aref[i]) > maxref) ? std::abs(Aref[i]) : maxref;
    norm = (absdiff > norm) ? absdiff : norm;
  }
  return norm / maxref;
}

double ftest1(double x, double y)
{
  double b = -2.258846861003648;
  double a1 = 0.5376671395461;
  double a2 = 1.833885014595086;
  return std::cos(2.0 * M_PI * b + a1 * x + a2 * y); 
}

double ftest2(double x, double y)
{
  double b1 = 0.726885133383238;
  double b2 = -0.303440924786016;
  double a1 = 0.488893770311789;
  double a2 = 1.034693009917860;
  return 1.0/((pow(1.0/a1,2)+pow(x-b1,2))*(pow(1.0/a2,2)+pow(y-b2,2)));
}

double ftest3(double x, double y)
{
  double b1 = -2.258846861003648;
  double b2 = 0.862173320368121;
  double a1 = 0.5376671395461;
  double a2 = 1.833885014595086;
  return std::exp(-(pow(a1,2)*pow(x-b1,2)+pow(a2,2)*pow(y-b2,2)));
}


double ftest4(double x, double y)
{
  double b1 = -2.258846861003648;
  double b2 = 0.862173320368121;
  double a1 = 0.5376671395461;
  double a2 = 1.833885014595086;
  return pow(a1 * x + a2 * y, 9);
}

double ftest5(double x, double y)
{
  double b1 = -2.258846861003648;
  double b2 = 0.862173320368121;
  double a1 = 0.5376671395461;
  double a2 = 1.833885014595086;
  return pow(a1 * x + a2 * y, 23);
}

double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
    if (grad) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}

typedef struct {
    double a, b;
} my_constraint_data;

double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data *d = (my_constraint_data *) data;
    double a = d->a, b = d->b;
    if (grad) {
        grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
        grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
 }

namespace
{



TEST(structureFactorTest, TolCheck)
{
  double tol  = 1e-15;
  unsigned int Np = 7;
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* H = (double*) calloc(Np*Np, sizeof(double));
  double* Href = (double*) calloc(Np*Np, sizeof(double));
  std::ifstream Hfile("../testdata/Href.txt");
  for (unsigned int i = 0; i < Np*Np; ++i)
  {
    Hfile >> Href[i];

  }
  structure_factors_tri<double>(Np, a, b, c, H);
  double diff = infnorm(H,Href,Np*Np);
  EXPECT_LT(diff, tol);
  //printMat(H,Np,Np);
  std::cout << "\nrelative infinity norm of diff = " << diff << "\n\n";
  free(H);
}

TEST(jPolyTest, TolCheck)
{
  double tol  = 1e-15;
  unsigned int Np = 7;
  unsigned int Nx = 5;
  double a = 0.5; double b = 0.5; 
  double* x = (double*) calloc(Nx, sizeof(double));
  double* V = (double*) malloc(Nx*Np*sizeof(double));
  double* Vref = (double*) calloc(Nx*Np, sizeof(double));
  std::ifstream Vfile("../testdata/Vref.txt");
  for (unsigned int i = 0; i < Nx*Np; ++i)
  {
    Vfile >> Vref[i];
  } 
  double h = 2.0/(Nx-1);
  for (unsigned int i = 0; i < Nx; ++i)
  {
    x[i] = -1.0 + h*i;
  }
  jPoly<double>(x, Nx, Np, 0.5, 0.5, V);  
  double diff = infnorm(V,Vref,Nx*Np);
  EXPECT_LT(diff, tol);
  //printMat(V,Nx,Np);
  std::cout << "\nrelative infinity norm of diff = " << diff << "\n\n";
  free(x);
  free(V);
}

TEST(jPolyTriTest, TolCheck)
{
  double tol  = 1e-12;
  double a = 0.5; double b = 0.5; double c = 0.5; 
  unsigned int N = 136;
  unsigned int m = 27;
  unsigned int Np = static_cast<unsigned int>(0.5 * m * (m + 1));
  double* Xk = (double*) calloc(N, sizeof(double));
  double* Yk = (double*) calloc(N, sizeof(double));
  double* Vm = (double*) calloc(N*Np, sizeof(double));
  double* Vmref = (double*) calloc(N*Np, sizeof(double));
  double* Hm = (double*) calloc(m*m, sizeof(double));
  // read in test data
  std::ifstream Zkfile("../testdata/Zkref.txt");
  std::ifstream Vmfile("../testdata/Vmref.txt");
  for (unsigned int i = 0; i < N; ++i) { Zkfile >> Xk[i]; }
  for (unsigned int i = 0; i < N; ++i) { Zkfile >> Yk[i]; }
  for (unsigned int i = 0; i < N*Np; ++i) { Vmfile >> Vmref[i]; }
  structure_factors_tri(m, a, b, c, Hm);
  jPoly_tri<double>(Xk, Yk, Hm, N, m-1, a, b, c, Vm);
  double diff = infnorm(Vm, Vmref, N*Np);
  EXPECT_LT(diff, tol);
  std::cout << "\nrelative infinity norm of diff = " << diff << "\n\n";
  free(Hm);
  free(Vm);
  free(Xk);
  free(Yk);
  free(Vmref);
}

#ifdef PLOT
namespace plt = matplotlibcpp;
#endif
TEST(jPolyTriConvTest, ConvPlotCheck)
{
  double tol = 1e-14;
  double (*farr[5])(double,double) = {ftest1, ftest2, ftest3, ftest4, ftest5};
  for (unsigned int iF = 0; iF < 5; ++iF)
  {
    unsigned int n, m, N;
    double I, Iref;
    N = 136;
    double* Xkref = (double*) malloc(N * sizeof(double));
    double* Ykref = (double*) malloc(N * sizeof(double));
    double* Wkref = (double*) malloc(N * sizeof(double));
    double* fref  = (double*) malloc(N * sizeof(double));
    double* errs  = (double*) malloc(14 * sizeof(double));  
    unsigned int * Ns = (unsigned int*) malloc(14 * sizeof(unsigned int));
    std::ifstream Zkreffile("../testdata/Zkref.txt");
    for (unsigned int i = 0; i < N; ++i) { Zkreffile >> Xkref[i]; }
    for (unsigned int i = 0; i < N; ++i) { Zkreffile >> Ykref[i]; }
    for (unsigned int i = 0; i < N; ++i) { Zkreffile >> Wkref[i]; }


    for (unsigned int i = 0; i < N; ++i)
    {
      fref[i] = farr[iF](Xkref[i],Ykref[i]);
    }
    
    Iref = vecdot(fref, Wkref, N);
     
    std::stringstream ss; 
     double Nms[14][2] = 
      { 
        {2, 3}, {3, 5}, 
        {4, 6}, {5,8},
        {6,10}, {7,12},
        {8,13}, {9,15},
        {10,17}, {11,18},
        {12,20}, {13,22},
        {14,23}, {15,24}
      };
    for (unsigned int j = 0; j < 14; ++j)
    {
      n = Nms[j][0];
      m = Nms[j][1];
      N = static_cast<unsigned int>(0.5 * n * (n+1)); Ns[j] = N;
      double* Xk = (double*) malloc(N * sizeof(double));
      double* Yk = (double*) malloc(N * sizeof(double));
      double* Wk = (double*) malloc(N * sizeof(double));
      double* f  = (double*) malloc(N * sizeof(double));
      ss << "../testdata/triquadLeg_" << n << "_" << m << ".txt";
      std::ifstream Zkfile(ss.str());
      for (unsigned int i = 0; i < N; ++i) { Zkfile >> Xk[i]; }
      for (unsigned int i = 0; i < N; ++i) { Zkfile >> Yk[i]; }
      for (unsigned int i = 0; i < N; ++i) { Zkfile >> Wk[i]; }

      for (unsigned int i = 0; i < N; ++i)
      {
        f[i] = farr[iF](Xk[i], Yk[i]);
      }

      I = vecdot(f, Wk, N);
      errs[j] = std::abs(I - Iref) / std::abs(Iref);
      if (j == 13)
      {
        EXPECT_LT(errs[j], tol);
      }
      ss.str(""); 
      free(Xk);
      free(Yk);
      free(Wk);
      free(f);
    }

    #ifdef PLOT
    std::vector<unsigned int> Nsvec(Ns,Ns+14);
    std::vector<double> errsvec(errs,errs+14);
    
    plt::semilogy(Nsvec, errsvec, "o--");
    #endif
    ss.clear();
    free(Xkref);
    free(Ykref);
    free(Wkref);
    free(fref);
    free(errs);
  }
  #ifdef PLOT
  plt::save("../testdata/intconv.png");
  #endif
}

TEST(promMatA1Test, TolCheck)
{
  double tol  = 1e-14;
  unsigned int n = 14;
  unsigned int n_k = n - 1; 
  unsigned int N = static_cast<unsigned int>(0.5 * (n_k + 1) * (n_k + 2)); 
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* H_abc = (double*) calloc((n+1)*(n+1), sizeof(double));
  double* H_a1bc = (double*) calloc((n+1)*(n+1), sizeof(double));
  double* K_a1bc = (double*) calloc(N*N, sizeof(double)); 
  double* K_a1bc_ref = (double*) calloc(N*N, sizeof(double)); 
 
  std::ifstream K_a1bc_ref_file("../testdata/K_a1bc_ref.txt");
  for (unsigned int i = 0; i < N*N; ++i) { K_a1bc_ref_file >> K_a1bc_ref[i]; }
  structure_factors_tri<double>(n+1, a, b, c, H_abc);
  structure_factors_tri<double>(n+1, a+1, b, c, H_a1bc);
  promotion_mat_tri(a, b, c, H_abc, H_a1bc, n_k, 0, K_a1bc);
  double diff = infnorm(K_a1bc, K_a1bc_ref, N*N);
  EXPECT_LT(diff, tol);
  std::cout << "\nrelative infinity norm of diff = " << diff << "\n\n";

  free(H_abc);
  free(H_a1bc);
  free(K_a1bc);
  free(K_a1bc_ref);
}

TEST(promMatB1Test, TolCheck)
{
  double tol  = 1e-14;
  unsigned int n = 14;
  unsigned int n_k = n - 1; 
  unsigned int N = static_cast<unsigned int>(0.5 * (n_k + 1) * (n_k + 2)); 
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* H_abc = (double*) calloc((n+1)*(n+1), sizeof(double));
  double* H_ab1c = (double*) calloc((n+1)*(n+1), sizeof(double));
  double* K_ab1c = (double*) calloc(N*N, sizeof(double)); 
  double* K_ab1c_ref = (double*) calloc(N*N, sizeof(double)); 
 
  std::ifstream K_ab1c_ref_file("../testdata/K_ab1c_ref.txt");
  for (unsigned int i = 0; i < N*N; ++i) { K_ab1c_ref_file >> K_ab1c_ref[i]; } 
  structure_factors_tri<double>(n+1, a, b, c, H_abc);
  structure_factors_tri<double>(n+1, a, b+1.0, c, H_ab1c);
  promotion_mat_tri(a, b, c, H_abc, H_ab1c, n_k, 1, K_ab1c);
  
  double diff = infnorm(K_ab1c, K_ab1c_ref, N*N);
  EXPECT_LT(diff, tol);
  std::cout << "\nrelative infinity norm of diff = " << diff << "\n\n";

  free(H_abc);
  free(H_ab1c);
  free(K_ab1c);
  free(K_ab1c_ref);
}

TEST(promMatC1Test, TolCheck)
{
  double tol  = 1e-14;
  unsigned int n = 14;
  unsigned int n_k = n - 1; 
  unsigned int N = static_cast<unsigned int>(0.5 * (n_k + 1) * (n_k + 2)); 
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* H_abc = (double*) calloc((n+1)*(n+1), sizeof(double));
  double* H_abc1 = (double*) calloc((n+1)*(n+1), sizeof(double));
  double* K_abc1 = (double*) calloc(N*N, sizeof(double)); 
  double* K_abc1_ref = (double*) calloc(N*N, sizeof(double)); 
 
  std::ifstream K_abc1_ref_file("../testdata/K_abc1_ref.txt");
  for (unsigned int i = 0; i < N*N; ++i) { K_abc1_ref_file >> K_abc1_ref[i]; } 
  structure_factors_tri<double>(n+1, a, b, c, H_abc);
  structure_factors_tri<double>(n+1, a, b+1.0, c, H_abc1);
  promotion_mat_tri(a, b, c, H_abc, H_abc1, n_k, 2, K_abc1);
  
  double diff = infnorm(K_abc1, K_abc1_ref, N*N);
  EXPECT_LT(diff, tol);
  std::cout << "\nrelative infinity norm of diff = " << diff << "\n\n";

  free(H_abc);
  free(H_abc1);
  free(K_abc1);
  free(K_abc1_ref);
}

TEST(dxTriTest, TolCheck)
{
  double tol  = 1e-14;
  unsigned int n = 14;
  unsigned int n_k = n - 1; 
  unsigned int N = static_cast<unsigned int>(0.5 * (n_k + 1) * (n_k + 2)); 
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* H_abc = (double*) calloc((n+1)*(n+1), sizeof(double));
  double* H_a1bc1 = (double*) calloc((n+1)*(n+1), sizeof(double));
  double* D_a1bc1 = (double*) calloc(N*N, sizeof(double)); 
  double* D_a1bc1_ref = (double*) calloc(N*N, sizeof(double)); 
 
  std::ifstream D_a1bc1_ref_file("../testdata/Dx_a1bc1_ref.txt");
  for (unsigned int i = 0; i < N*N; ++i) { D_a1bc1_ref_file >> D_a1bc1_ref[i]; } 
  structure_factors_tri<double>(n+1, a, b, c, H_abc);
  structure_factors_tri<double>(n+1, a+1.0, b, c+1.0, H_a1bc1);
  dxdy_mat_tri(a, b, c, H_abc, H_a1bc1, n_k, 0, D_a1bc1);
  
  double diff = infnorm(D_a1bc1, D_a1bc1_ref, N*N);
  EXPECT_LT(diff, tol);
  std::cout << "\nrelative infinity norm of diff = " << diff << "\n\n";

  free(H_abc);
  free(H_a1bc1);
  free(D_a1bc1);
  free(D_a1bc1_ref);
}

TEST(dyTriTest, TolCheck)
{
  double tol  = 1e-14;
  unsigned int n = 14;
  unsigned int n_k = n - 1; 
  unsigned int N = static_cast<unsigned int>(0.5 * (n_k + 1) * (n_k + 2)); 
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* H_abc = (double*) calloc((n+1)*(n+1), sizeof(double));
  double* H_ab1c1 = (double*) calloc((n+1)*(n+1), sizeof(double));
  double* D_ab1c1 = (double*) calloc(N*N, sizeof(double)); 
  double* D_ab1c1_ref = (double*) calloc(N*N, sizeof(double)); 
 
  std::ifstream D_ab1c1_ref_file("../testdata/Dy_ab1c1_ref.txt");
  for (unsigned int i = 0; i < N*N; ++i) { D_ab1c1_ref_file >> D_ab1c1_ref[i]; } 
  structure_factors_tri<double>(n+1, a, b, c, H_abc);
  structure_factors_tri<double>(n+1, a, b+1.0, c+1.0, H_ab1c1);
  dxdy_mat_tri(a, b, c, H_abc, H_ab1c1, n_k, 1, D_ab1c1);
  
  double diff = infnorm(D_ab1c1, D_ab1c1_ref, N*N);
  EXPECT_LT(diff, tol);
  std::cout << "\nrelative infinity norm of diff = " << diff << "\n\n";

  free(H_abc);
  free(H_ab1c1);
  free(D_ab1c1);
  free(D_ab1c1_ref);
}

TEST(jacobiMatTriTest, TolCheck)
{
  double tol  = 1e-14;
  unsigned int n = 20;
  unsigned int N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* Jn1 = (double*) calloc(N*N, sizeof(double));
  double* Jn2 = (double*) calloc(N*N, sizeof(double));
  double* Jn1_ref = (double*) calloc(N*N, sizeof(double));
  double* Jn2_ref = (double*) calloc(N*N, sizeof(double));
  std::ifstream Jn1_ref_file("../testdata/Jn1_ref.txt");
  std::ifstream Jn2_ref_file("../testdata/Jn2_ref.txt");
  for (unsigned int i = 0; i < N*N; ++i) 
  { 
    Jn1_ref_file >> Jn1_ref[i];
    Jn2_ref_file >> Jn2_ref[i];
  } 

  jacobi_mat_ON_tri<double>(n, a, b, c, Jn1, Jn2);
  
  double diff1 = infnorm(Jn1, Jn1_ref, N*N);  
  double diff2 = infnorm(Jn2, Jn2_ref, N*N);  
  
  EXPECT_LT(diff1, tol);
  EXPECT_LT(diff2, tol);
  std::cout << "\nrelative infinity norm of diff for Jx = " << diff1 << "\n\n";
  std::cout << "\nrelative infinity norm of diff for Jy = " << diff2 << "\n\n";
  
  free(Jn1);
  free(Jn2);
  free(Jn1_ref);
  free(Jn2_ref);

}

TEST(nloptExampleTest, runCheck)
{
  double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
  nlopt_opt opt;
  opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_min_objective(opt, myfunc, NULL);
  my_constraint_data data[2] = { {2,0}, {-1,1} };
  nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-8);
  nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-8);
  nlopt_set_xtol_rel(opt, 1e-4);
  double x[2] = { 1.234, 5.678 };  /* `*`some` `initial` `guess`*` */
  double minf; /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */
  if (nlopt_optimize(opt, x, &minf) < 0) {
      printf("nlopt failed!\n");
  }
  else {
      printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
  }
  nlopt_destroy(opt);
}

TEST(jacobigenTEST, tolCheck)
{

  double tol  = 1e-14;
  unsigned int n = 12;
  unsigned int N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* Jn1 = (double*) calloc(N*N, sizeof(double));
  double* Jn2 = (double*) calloc(N*N, sizeof(double));
  double* Jn1_ref = (double*) calloc(N*N, sizeof(double));
  double* Jn2_ref = (double*) calloc(N*N, sizeof(double));
  double* Z = (double*) calloc(3*N, sizeof(double));
  double* D = (double*) calloc(N*N, sizeof(double)); 
  double* V = (double*) calloc(N*N, sizeof(double));
  // need to figure out N pts and num polys needed for Pn (not Pn-1) 
  std::cout << 3*N << std::endl; 
  std::ifstream Z_file("../testdata/triquadLeg_15_24.txt");
  for (unsigned int i = 0; i < 3*N; ++i) 
  { 
    Z_file >> Z[i]; std::cout << Z[i] << std::endl; 
  }
 
  jacobi_mat_ON_tri<double>(n, a, b, c, Jn1, Jn2);

  //double diff1 = infnorm(Jn1, Jn1_ref, N*N);  
  //double diff2 = infnorm(Jn2, Jn2_ref, N*N);  
  //
  //EXPECT_LT(diff1, tol);
  //EXPECT_LT(diff2, tol);

  free(Jn1); free(Jn2);
  free(Jn1_ref); free(Jn2_ref);

}

} // end gtest namespace

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  std::cout << "\n\n<<<<<<< RUNNING ALL TESTS >>>>>>>\n\n";
  int ret{RUN_ALL_TESTS()};
  if (!ret)
      std::cout << "\n\n<<<<<<<<<<<< SUCCESS >>>>>>>>>>>>\n\n";
  else
      std::cout << "\n\n<<<<<<<<<<<< FAILURE >>>>>>>>>>>>\n\n";
  return 0;
}
