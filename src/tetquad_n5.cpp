#include<jevd.h>
#include<fstream>

int main(int argc, char* argv[])
{

  unsigned int N = 35;
  double a = 0.5, b = 0.5, c = 0.5, d = 0.5; 
  double* Jx = (double*) calloc(N*N, sizeof(double));
  double* Jy = (double*) calloc(N*N, sizeof(double));
  double* Jz = (double*) calloc(N*N, sizeof(double));
  double* J = (double*) calloc(N*N*3, sizeof(double));
  double* X0 = (double*) calloc(N, sizeof(double)); 
  double* Y0 = (double*) calloc(N, sizeof(double)); 
  double* Z0 = (double*) calloc(N, sizeof(double)); 

  std::ifstream Jxfile("../testing/testdata/J5x_tet.txt");
  std::ifstream Jyfile("../testing/testdata/J5y_tet.txt");
  std::ifstream Jzfile("../testing/testdata/J5z_tet.txt");
  for (unsigned int j = 0; j < N*N; ++j)
  {
    Jxfile >> Jx[j];
    Jyfile >> Jy[j];
    Jzfile >> Jz[j];
  }
  for (unsigned int j = 0; j < N; ++j)
  {
    for (unsigned int i = 0; i < N; ++i)
    {
      J[i + N*j]          = Jx[i + j*N];
      J[i + N*(j + N)]    = Jy[i + j*N];
      J[i + N*(j + 2*N)]  = Jz[i + j*N];
    }
  } 

  jointDiag<double>* jevd = new jointDiag(N, 3, 1e-10, J, 1); 
  std::ofstream Xfile("xtet.txt");
  std::ofstream Yfile("ytet.txt");
  std::ofstream Zfile("ztet.txt");

  for (unsigned int i = 0; i < N; ++i) 
  { 
    Xfile << J[i + N*i] << std::endl; 
    Yfile << J[i + N*(i + N)] << std::endl;
    Zfile << J[i + N*(i + 2*N)] << std::endl; 
  } 

  Xfile.close();
  Yfile.close();
  Zfile.close();
 
  free(Jx);
  free(Jy);
  free(Jz);
  free(J);
  free(X0);
  free(Y0);
  free(Z0);
 


}
