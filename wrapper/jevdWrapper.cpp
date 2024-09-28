#include<jevd.hh>

extern "C"
{
  jointDiag<double>* jevd ( unsigned int m,
                            unsigned int n,
                            double* J,
                            unsigned int nthreads,
                            double thresh = 1e-10,
                            bool hasV = false )
  {
    jointDiag<double>* jdiag = new jointDiag<double>( m, 
                                                      n, 
                                                      thresh, 
                                                      J,
                                                      nthreads,
                                                      hasV );
    return jdiag;
  }
} 
