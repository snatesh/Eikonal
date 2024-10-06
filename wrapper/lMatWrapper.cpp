#include<lMat.hh>

extern "C"
{

  void lMat ( double a, double b, double c, 
              const double* H, const double* H1, 
              unsigned int n, unsigned int mode, 
              double* L)
  {
    lMat<double>(a, b, c, H, H1, n, mode, L);
  }

}
