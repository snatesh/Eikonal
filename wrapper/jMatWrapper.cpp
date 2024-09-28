#include<jMat.hh>


extern "C"
{

  jMat<double>* jMat_T1 ( unsigned int n, 
                          double a, 
                          double b )
  {
    jMat<double>* mat = new jMat<double>(n, a, b); 
    return mat;
  }
  
  jMat<double>* jMat_T2 ( unsigned int n, 
                          double a, 
                          double b,
                          double c )
  {
    jMat<double>* mat = new jMat<double>(n, a, b, c);
    return mat;
  }

  jMat<double>* jMat_T3 ( unsigned int n, 
                          double a, 
                          double b,
                          double c, 
                          double d,
                          unsigned int nlg,
                          double* x, 
                          double* w,
                          unsigned int nthreads )
  {
    jMat<double>* mat = new jMat<double>( n, a, b, c, d, 
                                          nlg, x, w, nthreads);
    return mat;
  }

}

