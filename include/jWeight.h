#ifndef _JWEIGHT_H
#define _JWEIGHT_H

#include<iostream>
#include<cmath>
#include<type_traits>
#include<mapTensorQuad.h>


using std::tgamma;
using std::pow;

template<typename T, unsigned int dim> 
struct pVec 
{ 
  static_assert(dim >= 1 && dim <= 3, "ERROR: dim unsupported");  
  T pvec[dim+1]; 
  pVec(T* vec) 
  {
    if (dim == 1) 
    { 
      pvec[0] = vec[0]; pvec[1] = vec[1]; 
    }
    if (dim == 2) 
    { 
      pvec[0] = vec[0]; pvec[1] = vec[1];
      pvec[2] = vec[2];
    }
    if (dim == 3) 
    {
      pvec[0] = vec[0]; pvec[1] = vec[1];
      pvec[2] = vec[2]; pvec[3] = vec[3];
    }
  }
};


/*  Weight function for L2(T) with Jacobi basis, where
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
    pVec<T, dim> jparams; 
    
    T w(T* x)
    {
      return W(x, jparams) / wnorm;
    }    
    
    jWeight(pVec<T, dim> _jparams, T* _absc, T* _wght) 
      : jparams(_jparams), absc(_absc), wght(_wght)
    {
      bool supported = 1;
      if (dim > 1)
      {
        for (unsigned int i = 0; i <= dim; ++i)
        {
          supported *= (jparams.pvec[i] == 0.5);
        } 
        if (not supported)
        {
          std::cerr << "ERROR: Only legendre analog weight supported for dim > 1\n";
          exit(1);
        }
      }
      wnorm = 0.0;
      switch(dim)
      {
        case 1:
        {
          T X[1];
          #pragma omp simd reduction(+:wnorm)
          for (unsigned int i = 0; i < n; ++i)
          {
            X[0] = absc[i];
            wnorm +=  W(X, jparams) * (wght[i]);
          } 
          break;
        }
        case 2:
        {
          wnorm = 1.0 / 2.0;
          //T X[2];
          //#pragma omp simd collapse(2) reduction(+:wnorm)
          //for (unsigned int j = 0; j < n; ++j)
          //{
          //  for (unsigned int i = 0; i < n; ++i)
          //  {
          //    X[0] = absc[i];
          //    X[1] = absc[j];
          //    wnorm += W(X, jparams) * ( wght[i] * wght[j] );
          //  }
          //}
          break;
        }
        case 3:
        {
          wnorm = 1.0 / 6.0;
          //T X[3];
          //#pragma omp simd collapse(3) reduction(+:wnorm)
          //for (unsigned int k = 0; k < n; ++k)
          //{

          //  for (unsigned int j = 0; j < n; ++j)
          //  {

          //    for (unsigned int i = 0; i < n; ++i)
          //    {
          //      X[0] = absc[i];
          //      X[1] = absc[j];
          //      X[2] = absc[k];
          //      wnorm += W(X, jparams) * ( wght[i] * wght[j] * wght[k]);
          //    }
          //  }
          //}
          break;
        }
      }
    }
 
         
  private:
    unsigned int n = 20;

    T W(T* x, pVec<T, dim> jparams) 
    { 
      switch(dim)
      {
        case 1:
        {
          T a = jparams.pvec[0], b = jparams.pvec[1];
          if (a < -1 || b < -1)
          {
            std::cerr << "ERROR: Jacobi 1D params must be > -1\n";
            exit(1);
          }
          return  pow(1-x[0], jparams.pvec[0]) * 
                  pow(1+x[0], jparams.pvec[1]);
        }
        case 2:
        {
          return 1.0;
          //T a = jparams.pvec[0], b = jparams.pvec[1], c = jparams.pvec[2];
          //if (a < -0.5 || b < -0.5 || c < -0.5)
          //{
          //  std::cerr << "ERROR: Jacobi d>2 params must be > -1/2\n";
          //  exit(1);
          //}
          //return  pow(x[0], a - 0.5) * 
          //        pow(x[1], b - 0.5) * 
          //        pow(1 - x[0] - x[1], c - 0.5);
        }
        case 3:
        {
          return 1.0;
          //T a = jparams.pvec[0], b = jparams.pvec[1];
          //T c = jparams.pvec[2], d = jparams.pvec[3];
          //if (a < -0.5 || b < -0.5 || c < -0.5 || d < -0.5)
          //{
          //  std::cerr << "ERROR: Jacobi d>2 params must be > -1/2\n";
          //  exit(1);
          //}
          //return  pow(x[0], a - 0.5) * 
          //        pow(x[1], b - 0.5) * 
          //        pow(x[2], c - 0.5) * 
          //        pow(1 - x[0] - x[1] - x[2], d - 0.5);
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
