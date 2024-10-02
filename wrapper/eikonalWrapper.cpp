#include<eikonal.hh>
#include<cmath>
#include<iomanip>
extern "C"
{

  void printMat(const double* A, const unsigned int m, const unsigned int n)
  {
    for (unsigned int i = 0; i < m; ++i)
    {
      for (unsigned int j = 0; j < n; ++j)
      {
        std::cout << std::setw(10);
        std::cout << A[i + m*j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  
  void eikonalSolve ( unsigned int N,
                      unsigned int nl,
                      unsigned int nb,
                      unsigned int nh,
                      double* Dx, 
                      double* Dy,
                      double* rhs,
                      double* vl, 
                      double* vb, 
                      double* vh, 
                      double* cu_opt )
  {
    double *lb, *ub, *toleql, *toleqb, *toleqh; 
    lb      = (double*) malloc(N * sizeof(double));
    ub      = (double*) malloc(N * sizeof(double));
    toleql  = (double*) malloc(nl * sizeof(double)); 
    toleqb  = (double*) malloc(nb * sizeof(double)); 
    toleqh  = (double*) malloc(nh * sizeof(double)); 
    for (unsigned int i = 0; i < N; ++i)
    {
      ub[i]     = 10  ; 
      lb[i]     = -10 ; 
    }
    for (unsigned int i = 0; i < nl; ++i) { toleql[i] = 1e-14; }
    for (unsigned int i = 0; i < nb; ++i) { toleqb[i] = 1e-14; }
    for (unsigned int i = 0; i < nh; ++i) { toleqh[i] = 1e-14; }
    double tol = 1e-3;
    optData* data = new optData ( N, nl, nb, nh,
                                  Dx, Dy, rhs, 
                                  vl, vb, vh, cu_opt );  
    nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, N); 
    //nlopt_opt local_opt = nlopt_create(NLOPT_LN_SBPLX, N);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, F, data);
    nlopt_add_equality_mconstraint(opt, nl, cl, data, toleql);
    nlopt_add_equality_mconstraint(opt, nb, cb, data, toleqb);
    nlopt_add_equality_mconstraint(opt, nh, ch, data, toleqh);
    nlopt_set_xtol_rel(opt, tol);
    nlopt_set_stopval(opt, tol);
    //nlopt_set_xtol_rel(local_opt, tol);
    //nlopt_set_stopval(local_opt, tol);
    //nlopt_set_local_optimizer(opt, local_opt);
    double minF;
    if (nlopt_optimize(opt, data->cu, &minF) < 0) 
    {
      std::cerr << "NLOPT failed! Exiting .." << std::endl;
      exit(1);
    }
    count = 0;
    delete data; 
    nlopt_destroy(opt);
    //nlopt_destroy(local_opt);
  } 
}
