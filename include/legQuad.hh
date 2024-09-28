#ifndef _LEGQUAD_H
#define _LEGQUAD_H
#include<ngjquad.hh>

/* map integral from tetrahedron to cube and 
   use a tensor product quadrature rule on the cube
   to numerically approximate integral on T */

/* generate 1D legendre quadrature nodes
   and weights on [-1,1] */
template<typename T>
struct legQuad
{
  static_assert(std::is_same_v<T, double> == true || 
                std::is_same_v<T, float> == true, 
                "ERROR: Only double and float template\
                 instantiations are allowed");
  
  unsigned int nlg;
  T *x = 0, *w = 0; 
  
  legQuad ( unsigned int _n ) : nlg(_n)  
  {
    ngjQuad* gjquad = new ngjQuad  (  nlg, nlg, 0, 0, 1e-16, 1e-16, 0,
                                      0, 0, NLOPT_LN_SBPLX, 1 );

    gjquad->init();
    gjquad->runXW();
    x = (T*) calloc(nlg, sizeof(T));
    w = (T*) calloc(nlg, sizeof(T));
    // copy over nodes/weights and delete opt object
    #pragma omp simd
    for (unsigned int i = 0; i < nlg; ++i) 
    { 
      x[i] = (T) gjquad->optdata->Z0[i]; 
      // add factor of 2 for weight normalization
      w[i] = (T) gjquad->optdata->Z0[i + nlg] * 2.0;
    }
    delete gjquad;
  } 
  ~legQuad() 
  { 
    if (x) { free(x); x = 0; }
    if (w) { free(w); w = 0; }
  }
};

template struct legQuad<double>;
template struct legQuad<float>;
 
#endif
