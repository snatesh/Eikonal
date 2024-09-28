#include<jPoly.hh>

extern "C"
{

  jPoly<double>* jPoly_T3 ( unsigned int Nx,
                            unsigned int n,
                            double a,
                            double b,
                            double c,
                            double d,
                            unsigned int nthreads )
  {
    jPoly<double>* poly = new jPoly<double>(Nx, n, a, b, c, d, nthreads);
    return poly;                                                
  } 

  jPoly<double>* jPoly_T2 ( unsigned int Nx,
                            unsigned int n,
                            double a,
                            double b,
                            double c,
                            unsigned int nthreads )
  {
    jPoly<double>* poly = new jPoly<double>(Nx, n, a, b, c, nthreads);
    return poly;                                                
  } 

  jPoly<double>* jPoly_T1 ( unsigned int Nx,
                            unsigned int n,
                            double a,
                            double b,
                            unsigned int nthreads )
  {
    jPoly<double>* poly = new jPoly<double>(Nx, n, a, b, nthreads);
    return poly;                                                
  } 
}

