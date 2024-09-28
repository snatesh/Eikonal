#include<jWeight.hh>

extern "C"
{

  jWeight<double, 1>* jWeight_T1( double a, 
                                  double b,
                                  double* absc,
                                  double* wght )
  {
    double params[2] = {a, b};
    jWeight<double, 1>* jweight = new jWeight<double, 1>(params, absc, wght);
    return jweight;
  }
 
  jWeight<double, 2>* jWeight_T2( double a, 
                                  double b,
                                  double c,
                                  double* absc,
                                  double* wght )
  {
    double params[3] = {a, b, c};
    jWeight<double, 2>* jweight = new jWeight<double, 2>(params, absc, wght);
    return jweight;
  }

  jWeight<double, 3>* jWeight_T3( double a, 
                                  double b,
                                  double c,
                                  double d,
                                  double* absc,
                                  double* wght )
  {
    double params[4] = {a, b, c, d};
    jWeight<double, 3>* jweight = new jWeight<double, 3>(params, absc, wght);
    return jweight;
  }

}
