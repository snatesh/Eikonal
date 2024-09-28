#include<mapTensorQuad.hh>

extern "C"
{

  mapTensorQuad<double>* tensorQuad ( unsigned int nlg,
                                      double* x,
                                      double* w )
  {
    mapTensorQuad<double>* tquad = new mapTensorQuad<double>(nlg, x, w);
    return tquad; 
  }  

  void copy ( mapTensorQuad<double>* tquad,
              double* X,
              double* Y,
              double* Z,
              double* W )
  {
    tquad->copy(X, Y, Z, W);
  }
}
