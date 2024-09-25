#ifndef _JWEIGHT_H
#define _JWEIGHT_H

#include<iostream>
#include<cmath>
#include<mapQuad.h>
#include<type_traits>



using std::tgamma;
using std::pow;

/* Weight function for L2(T) with Jacobi basis, where
   T is the standard simplex */

template<typename T, unsigned int dim>
class jWeight
{
  static_assert(std::is_same_v<T, double> == true || 
                std::is_same_v<T, float> == true, 
                "ERROR: Only double and float template\
                 instantiations are allowed");
  
  static_assert(dim >= 1 && dim <= 3, "ERROR: dim unsupported");  

  public:
    T wnorm = 0; 
    T *absc = 0, *wght = 0;
    T jparams[dim+1]; 
    
    T w(T* x)
    {
      return W(x, jparams) / wnorm;
    }    

    jWeight(T* _jparams)
    {
      for (unsigned int i = 0; i < dim+1; ++i) 
      {
        jparams[i] = _jparams[i];
      }
      absc = (T*) calloc(n, sizeof(T));
      wght = (T*) calloc(n, sizeof(T));
      ngjQuad* gjquad = new ngjQuad  (  n, n, 0, 0, 1e-16, 1e-16, 0,
                                        0, 0, NLOPT_LN_SBPLX, 1 );
      gjquad->init();
      gjquad->runXW();
      for (unsigned int i = 0; i < n; ++i)
      {
        absc[i] = gjquad->optdata->Z0[i];
        wght[i] = gjquad->optdata->Z0[i + n] * 2.0;
      }
      delete gjquad;
      wnorm = 0.0;
      switch(dim)
      {
        case 1:
        {
          T X[1];
          #pragma omp simd reduction(+:wnorm)
          for (unsigned int i = 0; i < n; ++i)
          {
            X[0] = (T) absc[i];
            wnorm +=  W(X, jparams) * ( (T) wght[i]);
          } 
          break;
        }
        case 2:
        {
          T X[2];
          #pragma omp simd collapse(2) reduction(+:wnorm)
          for (unsigned int j = 0; j < n; ++j)
          {
            for (unsigned int i = 0; i < n; ++i)
            {
              X[0] = (T) absc[i];
              X[1] = (T) absc[j];
              wnorm += W(X, jparams) * ( wght[i] * wght[j] );
            }
          }
          break;
        }
        case 3:
        {
          T X[3];
          #pragma omp simd collapse(3) reduction(+:wnorm)
          for (unsigned int k = 0; k < n; ++k)
          {

            for (unsigned int j = 0; j < n; ++j)
            {

              for (unsigned int i = 0; i < n; ++i)
              {
                X[0] = (T) absc[i];
                X[1] = (T) absc[j];
                X[2] = (T) absc[k];
                wnorm += W(X, jparams) * ( wght[i] * wght[j] * wght[k]);
              }
            }
          }
          break;
        }
      }
    }
 

  ~jWeight() 
  { 
    if (absc) { free(absc); absc = 0;}
    if (wght) { free(wght); wght = 0;}
  }
         
  private:
    unsigned int n = 20;

    T W(T* x, T* jparams) 
    { 
      switch(dim)
      {
        case 1:
        {
          T a = jparams[0], b = jparams[1];
          if (a < -1 || b < -1)
          {
            std::cerr << "ERROR: Jacobi 1D params must be > -1\n";
            exit(1);
          }
          return  pow(1-x[0], jparams[0]) * 
                  pow(1+x[0], jparams[1]);
        }
        case 2:
        {
          T a = jparams[0], b = jparams[1], c = jparams[2];
          if (a < -0.5 || b < -0.5 || c < -0.5)
          {
            std::cerr << "ERROR: Jacobi d>2 params must be > -1/2\n";
            exit(1);
          }
          return  pow(x[0], a - 0.5) * 
                  pow(x[1], b - 0.5) * 
                  pow(1 - x[0] - x[1], c - 0.5);
        }
        case 3:
        {
          T a = jparams[0], b = jparams[1], c = jparams[2], d = jparams[3];
          if (a < -0.5 || b < -0.5 || c < -0.5 || d < -0.5)
          {
            std::cerr << "ERROR: Jacobi d>2 params must be > -1/2\n";
            exit(1);
          }
          return  pow(x[0], a - 0.5) * 
                  pow(x[1], b - 0.5) * 
                  pow(x[2], c - 0.5) * 
                  pow(1 - x[0] - x[1] - x[2], d - 0.5);
        }
      }
    }
};

template class jWeight<double,1>;
template class jWeight<double,2>;
template class jWeight<double,3>;
template class jWeight<float,1>;
template class jWeight<float,2>;
template class jWeight<float,3>;

#endif
