#include<dMat.hh>

extern "C"
{

  void dMat ( double a, double b, double c, 
              const double* H, const double* H1, 
              unsigned int n, unsigned int mode, 
              double* D, bool weighted)
  {
    dMat<double>(a, b, c, H, H1, n, mode, D, weighted);
  }
}
