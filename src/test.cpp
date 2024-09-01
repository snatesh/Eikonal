#include<iostream>
#include"structure_factors.h"

int main(int argc, char* argv[])
{
  unsigned int N = 5;
  double a = 0.5; double b = 0.5; double c = 0.5;
  double* H = (double*) calloc(N*N, sizeof(double));

  structure_factors_tri(N, a, b, c, H);
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      std::cout << H[i + N*j] << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
