#include<eikonal.hh>
#include<cmath>

extern "C"
{
  void eikonalSolve ( unsigned int N,
                      double* Dx, 
                      double* Dy,
                      double* vl, 
                      double* vb, 
                      double* vh, 
                      double* cu_opt )
  {
    double *lb, *ub, *toleq, tol = 1e-14;
    lb    = (double*) malloc(N * sizeof(double));
    ub    = (double*) malloc(N * sizeof(double));
    toleq = (double*) malloc(N * sizeof(double)); 
    for (unsigned int i = 0; i < N; ++i)
    {
      ub[i]     = HUGE_VAL  ;
      lb[i]     = -HUGE_VAL ; 
      toleq[i]  = 1e-6       ;
    }
    optData* data = new optData ( N, Dx, Dy, vl, vb, vh, cu_opt );  
    
    nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, N); 
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, F, data);
    nlopt_add_equality_mconstraint(opt, N, cl, data, toleq);
    nlopt_add_equality_mconstraint(opt, N, cb, data, toleq);
    nlopt_add_equality_mconstraint(opt, N, ch, data, toleq);
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
  } 
}
