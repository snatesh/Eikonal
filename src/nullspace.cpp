#include<cblas.h>
#include<cstdlib>
#include<random>
#include<iostream>
#include<cmath>
#include<lapacke.h>
#include<fstream>
#include<iomanip>

void printMat(const double* A, const unsigned int m, const unsigned int n)
{
  std::cout << std::endl;
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

void transpose  ( const double* A, 
                  const unsigned int m, 
                  const unsigned int n,
                  double* AT )
{
  for (unsigned int j = 0; j < n; ++j)
  {
    for (unsigned int i = 0; i < m; ++i)
    {
      AT[j + n*i] = A[i + m*j];
    }
  }
}
              

int main(int argc, char* argv[])
{
  unsigned int M = 7;
  unsigned int N = 10;
  unsigned int dimS = M;

  double* A = (double*) calloc(M*N, sizeof(double));
  double* AT = (double*) calloc(M*N, sizeof(double));
  double* U = (double*) calloc(M*M, sizeof(double));
  double* UT = (double*) calloc(M*M, sizeof(double));
  double* S = (double*) calloc(dimS, sizeof(double));
  double* superb = (double*) calloc(dimS, sizeof(double));
  double* VT = (double*) calloc(N*N, sizeof(double));  
  double* V = (double*) calloc(N*N, sizeof(double));
  std::ifstream Afile("../testing/testdata/A.txt");
  for (unsigned int i = 0; i < M*N; ++i)
  {
    Afile >> A[i];
  }
  
  //printMat(A, M, N);
  transpose(A, M, N, AT);
  //printMat(AT, N, M); 
  //LAPACKE_dgesvd ( LAPACK_COL_MAJOR, 'A', 'A', 
  //                 M, N, A, M, S, 
  //                 U, M, VT, N, superb );
  LAPACKE_dgesvd ( LAPACK_COL_MAJOR, 'A', 'A', 
                   N, M, AT, N, S, 
                   V, N, UT, M, superb );
  printMat(V, N, N);
  //transpose(VT, N, N, V);
  //printMat(V, N, N);
  
  double* nullA = &V[N*M];
  printMat(nullA, N, N-M);

  double* test = (double*) calloc(M*(N-M), sizeof(double));

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
              M, N-M, N, 1.0, A, M, nullA, N, 0.0, test, M);
  

  for (unsigned int i = 0; i < M*(N-M); ++i)
  {
    std::cout << test[i] << std::endl;
  }
  free(test);
  free(A);
  free(AT);
  free(U);
  free(UT);
  free(S);
  free(superb);
  free(VT);
  free(V);
  

  return 0;

}
