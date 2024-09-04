#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#include<cblas.h>
#include"gtest/gtest.h"
#include"structure_factors.h"
#include"jPoly.h"
#include"promotion_mat_tri.h"
#include"../include/matplotlibcpp.h"


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
  //double sum = 0.0;
  //for (unsigned int i = 0; i < N; ++i)
  //{
  //  sum += A[i] * B[i];
  //} 
  //return sum;
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

double Nms[14][2] = {{2, 3}, {3, 5}, {4, 6},{5,8},{6,10},{7,12},{8,13},{9,15},{10,17},{11,18},{12,20},{13,22},{14,23},{15,24}};

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

namespace plt = matplotlibcpp;
TEST(jPolyTriConvTest, PlotCheck)
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

    std::vector<unsigned int> Nsvec(Ns,Ns+14);
    std::vector<double> errsvec(errs,errs+14);

    plt::semilogy(Nsvec, errsvec, "o--");

    ss.clear();
    free(Xkref);
    free(Ykref);
    free(Wkref);
    free(fref);
    free(errs);
  }
  plt::save("../testdata/intconv.png");
}

TEST(promMatA1Test, TolCheck)
{
  double tol  = 1e-14;
  unsigned int n = 13;
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
