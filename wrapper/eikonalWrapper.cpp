#include<eikonal.hh>
#include<cmath>
#include<iomanip>
extern "C"
{

  eikonal* createSolver ( unsigned int _n,
                          unsigned int _N,
                          unsigned int _nl,
                          unsigned int _nb,
                          unsigned int _nh,
                          double* frhs,
                          double* fu,
                          double* _X,
                          double* _Y, 
                          double* _W,
                          unsigned int _nthreads )
  {
    return new eikonal  ( _n, _n, _nl, _nb, _nh,
                          frhs, fu, _X, _Y, _W,
                          _nthreads );
  }
  
  void getCoeffs(eikonal* solver, double* cu)
  {
    for (unsigned int i = 0; i < solver->Np; ++i)
    {
      cu[i] = solver->cu[i]; 
    }
  }


  void eikonalSolveH  ( unsigned int N,
                        double* Dx, 
                        double* Dy,
                        double* Hess,
                        double* rhs,
                        double* cu_opt )
  {
    double *lb, *ub;
    lb      = (double*) malloc(N * sizeof(double));
    ub      = (double*) malloc(N * sizeof(double));
    for (unsigned int i = 0; i < N; ++i)
    {
      ub[i]     = HUGE_VAL  ; 
      lb[i]     = -HUGE_VAL ; 
    }
    double tol = 1e-6;
    optDataEik* data = new optDataEik ( N, Dx, Dy, Hess, rhs, cu_opt );  
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
                        unsigned int ne,
                        double* Dx, 
                        double* Dy,
                        double* Hess,
                        double* rhs,
                        double* Ve, 
                        double* cu_opt )
  {
    double *lob, *upb, *toleq; 
    lob   = (double*) malloc(N * sizeof(double));
    upb   = (double*) malloc(N * sizeof(double));
    toleq = (double*) malloc(ne * sizeof(double)); 
    double tol = 1e-14;
    for (unsigned int i = 0; i < N; ++i)
    {
      upb[i]     = 1000 ; 
      lob[i]     = -1000; 
    }
    for (unsigned int i = 0; i < ne; ++i) { toleq[i] = tol; }

    optDataEik* data = new optDataEik ( N, Dx, Dy, Hess, rhs, cu_opt, 
                                  ne, Ve); 
    nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, N); 
    //nlopt_opt local_opt = nlopt_create(NLOPT_LD_CCSAQ, N);
    nlopt_set_lower_bounds(opt, lob);
    nlopt_set_upper_bounds(opt, upb);
    nlopt_set_min_objective(opt, F, data);
    nlopt_add_equality_mconstraint(opt, ne, cl, data, toleq);
    nlopt_set_xtol_rel(opt, tol);
    nlopt_set_ftol_rel(opt, tol);
    //nlopt_set_xtol_rel(local_opt, tol);
    //nlopt_set_ftol_rel(local_opt, tol);
    //nlopt_set_local_optimizer(opt, local_opt);
    double minF;
    nlopt_result result = nlopt_optimize(opt, data->cu, &minF);
    if (result < 0 && result != -4) 
    {
      std::cerr << "NLOPT failed! \n"; 
      std::cerr << nlopt_result_to_string(result);
      std::cerr << "\n Exiting ..\n";
      exit(1);
    }
    else
    {
      if (result == -4) 
      { 
        std::cout << nlopt_result_to_string(result) << std::endl;
      }
      std::cout << "NLOPT converged with residual " << minF << std::endl;
    }
    count = 0;
    delete data; 
    nlopt_destroy(opt);
    //nlopt_destroy(local_opt);
  } 
}
