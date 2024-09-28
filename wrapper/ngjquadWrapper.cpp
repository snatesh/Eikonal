#include<ngjquad.hh>

extern "C"
{

  ngjQuad* ngjquad_T1 ( unsigned int n, 
                        unsigned int m, 
                        double a, 
                        double b, 
                        double tol, 
                        double tolc, 
                        double alpha,
                        bool use_newton, 
                        bool use_wolfe,
                        nlopt_algorithm alg, 
                        unsigned int nthreads )
  {
    ngjQuad* gjquad = new ngjQuad ( n, m, a, b, tol, tolc,
                                    alpha, use_newton, use_wolfe,
                                    alg, nthreads );
    return gjquad; 
  }
  
  ngjQuad* ngjquad_T2 ( unsigned int n, 
                        unsigned int m, 
                        double a, 
                        double b, 
                        double c, 
                        double tol, 
                        double tolc, 
                        double alpha,
                        bool use_newton, 
                        bool use_wolfe,
                        nlopt_algorithm alg, 
                        unsigned int nthreads )
  {
    ngjQuad* gjquad = new ngjQuad ( n, m, a, b, c, tol, tolc,
                                    alpha, use_newton, use_wolfe,
                                    alg, nthreads );
    return gjquad; 
  }

  ngjQuad* ngjquad_T3 ( unsigned int n, 
                        unsigned int m, 
                        double a, 
                        double b, 
                        double c, 
                        double d,
                        double* legx,
                        double* legw,
                        double tol, 
                        double tolc, 
                        double alpha,
                        bool use_newton, 
                        bool use_wolfe,
                        nlopt_algorithm alg, 
                        unsigned int nthreads )
  {
    ngjQuad* gjquad = new ngjQuad ( n, m, a, b, c, d, legx, legw, 
                                    tol, tolc, alpha, use_newton, 
                                    use_wolfe, alg, nthreads );
    return gjquad; 
  }

  void copyGJQuad (ngjQuad* gjquad, double* quad)
  {
    double* Z0 = gjquad->optdata->Z0;
    unsigned int N = gjquad->optdata->N;
    unsigned int dim = gjquad->optdata->dim;
    for (unsigned int i = 0; i < N*dim; ++i)
    {
      quad[i] = Z0[i];
    }
  }
}
