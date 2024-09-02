#include<iomanip>
#include"gtest/gtest.h"
#include"structure_factors.h"
#include"jPoly.h"

/* 
  Top of file is reserved for test helper functions and data definition

  gtest routines and main are at bottom of file 
*/

void printMat(const double* A, const unsigned int m, const unsigned int n)
{

  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      std::cout << std::setw(10);
      std::cout << A[i + m*j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

double diffMat(const double* A, const double* B, 
             const unsigned int m, const unsigned int n)
{

  double diff = 0.0;
  for (int i = 0; i < m; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      diff += std::abs(A[i + m*j] - B[i + m*j]);

    }
  }
  return diff;
}

double Href[49] =
{
1.000000000000000,
                   0,
                   0,
                   0,
                   0,
                   0,
                   0,
   0.707106781186548,
   0.408248290463863,
                   0,
                   0,
                   0,
                   0,
                   0,
   0.577350269189626,
   0.333333333333333,
   0.258198889747161,
                   0,
                   0,
                   0,
                   0,
   0.500000000000000,
   0.288675134594813,
   0.223606797749979,
   0.188982236504614,
                   0,
                   0,
                   0,
   0.447213595499958,
   0.258198889747161,
   0.200000000000000,
   0.169030850945703,
   0.149071198499986,
                   0,
                   0,
   0.408248290463863,
   0.235702260395516,
   0.182574185835055,
   0.154303349962092,
   0.136082763487954,
   0.123091490979333,
                   0,
   0.377964473009227,
   0.218217890235992,
   0.169030850945703,
   0.142857142857143,
   0.125988157669742,
   0.113960576459638,
   0.104828483672192,
};

double Vref[35] =
{
   1.000000000000000,
   1.000000000000000,
   1.000000000000000,
   1.000000000000000,
   1.000000000000000,
  -1.500000000000000,
  -0.750000000000000,
                   0,
   0.750000000000000,
   1.500000000000000,
   1.875000000000000,
                   0,
  -0.625000000000000,
                   0,
   1.875000000000000,
  -2.187500000000000,
   0.546875000000000,
                   0,
  -0.546875000000000,
   2.187500000000000,
   2.460937500000000,
  -0.492187500000000,
   0.492187500000000,
  -0.492187500000000,
   2.460937500000000,
  -2.707031250000000,
                   0,
                   0,
                   0,
   2.707031250000000,
   2.932617187500000,
   0.418945312500000,
  -0.418945312500000,
   0.418945312500000,
   2.932617187500000,

};


namespace
{

TEST(structureFactorTest, TolCheck)
{
  double tol  = 1e-14;
  unsigned int Np = 7;
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* H = (double*) calloc(Np*Np, sizeof(double));
  structure_factors_tri(Np, a, b, c, H);
  EXPECT_LT(diffMat(H, Href, Np, Np), tol);
  free(H);
}

TEST(jPolyTest, TolCheck)
{
  double tol  = 1e-14;
  unsigned int Np = 7;
  unsigned int Nx = 5;
  double a = 0.5; double b = 0.5; 
  double* x = (double*) calloc(Nx, sizeof(double));
  double* V = (double*) malloc(Nx*Np*sizeof(double));
  double h = 2.0/(Nx-1);
  for (unsigned int i = 0; i < Nx; ++i)
  {
    x[i] = -1.0 + h*i;
  }
  jPoly<double>(x, Nx, Np, 0.5, 0.5, V);  
  EXPECT_LT(diffMat(V, Vref, Nx, Np), tol);
  free(x);
  free(V);
}

}

int main(int argc, char* argv[])
{

  ::testing::InitGoogleTest(&argc, argv);
  std::cout << "\n\nRUNNING ALL TESTS ..." << std::endl;
  int ret{RUN_ALL_TESTS()};
  if (!ret)
      std::cout << "<<<SUCCESS>>>" << std::endl;
  else
      std::cout << "FAILED" << std::endl;
  return 0;
}
