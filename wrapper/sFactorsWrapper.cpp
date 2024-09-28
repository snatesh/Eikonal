#include<sFactors.hh>

extern "C"
{
  double Pochhammer(const double x, const unsigned int n)
  {
    return pochhammer<double>(x, n); 
  }

  double sFactor  ( const unsigned int n, 
                    const unsigned int k,
                    const unsigned int j,
                    const double a, 
                    const double b,
                    const double c, 
                    const double d  )
  {
    return sFactors<double>(n, k, j, a, b, c, d); 
  }

  void sFactors  ( const unsigned int N,
                   const double a, const double b, 
                   const double c, double* H )
  {
    sFactors<double>(N, a, b, c, H);
  }

}
