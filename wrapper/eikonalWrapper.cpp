#include<eikonal.hh>
#include<cmath>
#include<iomanip>
extern "C"
{

  void eikonalSolveH  ( unsigned int N,
                        double* Dx, 
                        double* Dy,
                        double* rhs,
                        double* cu_opt )
  {
    double *lb, *ub, *toleql, *toleqb, *toleqh; 
    lb      = (double*) malloc(N * sizeof(double));
    ub      = (double*) malloc(N * sizeof(double));
    for (unsigned int i = 0; i < N; ++i)
    {
      ub[i]     = HUGE_VAL  ; 
      lb[i]     = -HUGE_VAL ; 
    }
    double tol = 1e-6;
    optData* data = new optData ( N, Dx, Dy, rhs, cu_opt );  
    nlopt_opt opt = nlopt_create(NLOPT_LN_SBPLX, N); 
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, F, data);
    nlopt_set_xtol_rel(opt, tol);
    nlopt_set_stopval(opt, tol);
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

  void eikonalSolveP  ( unsigned int N,
                        unsigned int nl,
                        unsigned int nb,
                        unsigned int nh,
                        double* Dx, 
                        double* Dy,
                        double* rhs,
                        double* vl, double* ul, 
                        double* vb, double* ub,
                        double* vh, double* uh,
                        double* cu_opt )
  {
    double *lob, *upb, *toleql, *toleqb, *toleqh; 
    lob    = (double*) malloc(N * sizeof(double));
    upb     = (double*) malloc(N * sizeof(double));
    toleql  = (double*) malloc(nl * sizeof(double)); 
    toleqb  = (double*) malloc(nb * sizeof(double)); 
    toleqh  = (double*) malloc(nh * sizeof(double)); 
    for (unsigned int i = 0; i < N; ++i)
    {
      upb[i]     = HUGE_VAL  ; 
      lob[i]     = -HUGE_VAL ; 
    }
    for (unsigned int i = 0; i < nl; ++i) { toleql[i] = 1e-14; }
    for (unsigned int i = 0; i < nb; ++i) { toleqb[i] = 1e-14; }
    for (unsigned int i = 0; i < nh; ++i) { toleqh[i] = 1e-14; }
    double tol = 1e-8;
    optData* data = new optData ( N, Dx, Dy, rhs, cu_opt, 
                                  nl, nb, nh, 
                                  vl, vb, vh, 
                                  ul, ub, uh );  
    nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, N); 
    //nlopt_opt local_opt = nlopt_create(NLOPT_LN_SBPLX, N);
    nlopt_set_lower_bounds(opt, lob);
    nlopt_set_upper_bounds(opt, upb);
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
