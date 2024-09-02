#include<iomanip>
#include"structure_factors.h"
#include"jPoly.h"




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

int main(int argc, char* argv[])
{
  unsigned int Np = 7;
  unsigned int Nx = 5;
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* H = (double*) calloc(Np*Np, sizeof(double));
  double* x = (double*) calloc(Nx, sizeof(double));
  double h = 2.0/(Nx-1);
  for (unsigned int i = 0; i < Nx; ++i)
  {
    x[i] = -1.0 + h*i;
    std::cout << x[i] << std::endl;
  }
  double* V = (double*) malloc(Nx*Np*sizeof(double));

  structure_factors_tri(Np, a, b, c, H);
  printMat(H,Np,Np);
  jPoly<double>(x,Nx,Np,0.5,0.5,V);  
  printMat(V,Nx,Np);
  
  unsigned int tmp = static_cast<unsigned int>(0.5 * (Np + 1) * (Np + 2));
  std::cout << tmp << std::endl;

  free(H);
  free(x);
  free(V);

  return 0;
}
