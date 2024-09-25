#ifndef _MAPQUAD_H
#define _MAPQUAD_H
#include<ngjquad.h>

/* define deformation map between standard tetrahedron 
   and unit cube in the first quadrant of (x0,x1,x2) cartesian system */
template<typename T>
void f  ( T* x, T* F  )
{
  T linf = *std::max_element(x, x + 3);
  T l1 = x[0] + x[1] + x[2];
  F[0] = x[0] * linf / l1;
  F[1] = x[1] * linf / l1;
  F[2] = x[2] * linf / l1;
}

template<typename T>
void gradF  ( T* x, T* dF )
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
template<typename T>
T jDet  ( T* x, T* dF )
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

template void f<double>(double*, double*);
template void gradF<double>(double*, double*);
template double jDet<double>(double*, double*);


/* map integral from tetrahedron to cube and 
   use a tensor product quadrature rule on the cube
   to numerically approximate integral on T */
struct intTonC
{
  unsigned int n;
  double sum = 0.0;
  intTonC  ( unsigned int _n) : n(_n)
  {
    ngjQuad* gjquad = new ngjQuad  (  n, n, 0, 0, 1e-16, 1e-16, 0,
                                      0, 0, NLOPT_LN_SBPLX, 1 );

    gjquad->init();
    gjquad->runXW();
    double* x = gjquad->optdata->Z0;
    // factor of two missing below, 
    // but is removed due to rescaling anyway
    double* w = x + n; 
    // scale rule to [0,1]
    for (unsigned int i = 0; i < n; ++i) 
    {
      x[i]  = (x[i]+1)/2.0;
    } 

    double W, Xc[3], Xt[3], c2tJ, dF[9]; sum = 0;
    #pragma omp simd collapse(3) reduction(+:sum)
    for (unsigned int k = 0; k < n; ++k)
    {
      for (unsigned int j = 0; j < n; ++j)
      {
        for (unsigned int i = 0; i < n; ++i)
        {
          
          Xc[0] = x[i]; Xc[1] = x[j]; Xc[2] = x[k];
          f(Xc, Xt);
          W = w[i] * w[j] * w[k];
          c2tJ = jDet(Xc, dF);
          sum += integrand(Xt) * W * c2tJ;      
        }
      } 
    }
    delete gjquad;
  }

  double integrand(double* x)
  {
    return std::sin(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  }
};




#endif
