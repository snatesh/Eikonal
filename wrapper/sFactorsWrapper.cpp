#include<sFactors.hh>

extern "C"
{

  double sFactors_T3  ( const unsigned int n, 
                        const unsigned int k,
                        const unsigned int j,
                        const double a, 
                        const double b,
                        const double c, 
                        const double d )
  {
    double fac = sFactors<double>(n, k, j, a, b, c, d);
    return fac;
  }

  void sFactors_T2  ( const unsigned int N,
                      const double a, const double b, 
                      const double c, double* H )
  {
    sFactors<double>(N, a, b, c, H);
  }

  void sFactors_T1  ( const unsigned int N,
                      const double a, const double b, 
                      double* H)
  {
    structure_factors<double>(N, a, b, H);
  }

}
