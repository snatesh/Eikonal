#ifndef _MAPTENSORQUAD_H
#define _MAPTENSORQUAD_H
#include<algorithm>
#include<omp.h>

/* map integral from tetrahedron to cube and 
   use a tensor product quadrature rule on the cube
   to numerically approximate integral on T 

  Inputs: 
    nlg   - number of legendre nodes and weights
    x, w  - legendre nodes and weights on [-1,1]

  Outputs:
    X,Y,Z,W - nlg^3 nodes and weights for the 
              Tetrahedron (x>=0,y>=0,z>=0,1-x-y-z>=0) 
            - generated on call to map()

  TODO: Unify for d=1,2. Currently, only d=3 is supported
*/


template<typename T>
struct mapTensorQuad
{
  unsigned int nlg;
  T *x = 0, *w = 0;
  T *X = 0, *Y = 0, *Z = 0, *W = 0;
  
  mapTensorQuad(unsigned int _nlg, T* _x, T* _w)
    : nlg(_nlg)
  {
    x = (T*) calloc(nlg, sizeof(T));
    w = (T*) calloc(nlg, sizeof(T));
    X = (T*) calloc(nlg*nlg*nlg, sizeof(T)); 
    Y = (T*) calloc(nlg*nlg*nlg, sizeof(T)); 
    Z = (T*) calloc(nlg*nlg*nlg, sizeof(T)); 
    W = (T*) calloc(nlg*nlg*nlg, sizeof(T)); 
    // copy nodes and weights and scale to (0,1)
    for (unsigned int i = 0; i < nlg; ++i)
    {
      x[i] = (_x[i] + 1.0) / 2.0;
      w[i] = (_w[i] / 2.0);
    }
    map();
  } 

  ~mapTensorQuad()
  {
    if (x) { free(x); x = 0; }
    if (w) { free(w); w = 0; }
    if (X) { free(X); X = 0; }
    if (Y) { free(Y); Y = 0; }
    if (Z) { free(Z); Z = 0; }
    if (W) { free(W); W = 0; }
  }

  inline void map()
  {
    T xc[3], xt[3], dF[9], c2tJ; 
    unsigned int ind = 0;
    #pragma omp simd collapse(3) 
    for (unsigned int k = 0; k < nlg; ++k)
    {
      for (unsigned int j = 0; j < nlg; ++j)
      {
        for (unsigned int i = 0; i < nlg; ++i)
        {
          // point in cube
          xc[0] = x[i]; xc[1] = x[j]; xc[2] = x[k];
          // xt is point in tet
          f(xc, xt); 
          // Jacobian at xc 
          c2tJ = jDet(xc, dF);
          // save combined weight at xt
          W[ind] = w[i] * w[j] * w[k] * c2tJ;
          // save xt points
          X[ind] = xt[0];
          Y[ind] = xt[1];
          Z[ind] = xt[2];
          ind += 1;
        }
      }
    }
  } 
 
  /* deformation map between standard tetrahedron 
     and unit cube in the first quadrant of R3 with
     standard basis */
  #pragma omp declare simd
  inline void f  ( T* x, T* F  )
  {
    T linf = *std::max_element(x, x + 3);
    T l1 = x[0] + x[1] + x[2];
    F[0] = x[0] * linf / l1;
    F[1] = x[1] * linf / l1;
    F[2] = x[2] * linf / l1;
  }
 
  #pragma omp declare simd 
  inline void gradF  ( T* x, T* dF )
  {
    T h = 1e-7;
    T xph[3] = {x[0] + h, x[1], x[2]};
    T yph[3] = {x[0], x[1] + h, x[2]};
    T zph[3] = {x[0], x[1], x[2] + h};
    T xmh[3] = {x[0] - h, x[1], x[2]};
    T ymh[3] = {x[0], x[1] - h, x[2]};
    T zmh[3] = {x[0], x[1], x[2] - h};
    T Fxph[3], Fyph[3], Fzph[3];
    T Fxmh[3], Fymh[3], Fzmh[3];
    // f(x+h,y,z), f(x-h,y,z), f(x,y+h,z)
    f(xph, Fxph); f(xmh, Fxmh); f(yph, Fyph);
    // f(x,y-h,z), f(x,y,z+h), f(x,y,z-h)
    f(ymh, Fymh); f(zph, Fzph); f(zmh, Fzmh);
    // dxF
    dF[0] = (Fxph[0]-Fxmh[0]) / (2 * h);
    dF[1] = (Fxph[1]-Fxmh[1]) / (2 * h);
    dF[2] = (Fxph[2]-Fxmh[2]) / (2 * h);
    // dyF
    dF[3] = (Fyph[0]-Fymh[0]) / (2 * h);   
    dF[4] = (Fyph[1]-Fymh[1]) / (2 * h);   
    dF[5] = (Fyph[2]-Fymh[2]) / (2 * h);   
    // dzF
    dF[6] = (Fzph[0]-Fzmh[0]) / (2 * h);
    dF[7] = (Fzph[1]-Fzmh[1]) / (2 * h);
    dF[8] = (Fzph[2]-Fzmh[2]) / (2 * h);
  }
  
  #pragma omp declare simd
  inline T jDet  ( T* x, T* dF )
  {
    // limiting value through origin (line alpha*(1,1,1))
    if (x[0] == x[1] && x[0] == x[2]) { return 1.0/27.0; }
    // otherwise compute by finite difference 
    gradF(x, dF);
    
    /*      | df0 df3 df6 |
       det[ | df1 df4 df7 | ]
            | df2 df5 df8 |     */ 
    
    T det1 = dF[0] * (dF[4] * dF[8] - dF[7] * dF[5]);
    T det2 = dF[3] * (dF[1] * dF[8] - dF[7] * dF[2]);
    T det3 = dF[6] * (dF[1] * dF[5] - dF[4] * dF[2]);
     
    return det1 - det2 + det3;
  
  } 

};

template struct mapTensorQuad<double>;
template struct mapTensorQuad<float>;


#endif
