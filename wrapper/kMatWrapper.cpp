#include<kMat.hh>

extern "C"
{
  void kMat ( double a, double b, double c, const double* H,
              const double* H1, unsigned int n, 
              unsigned int mode, double* K )
  {
    kMat<double>(a, b, c, H, H1, n, mode, K);
  } 
}
