#ifndef _EIKONAL_H
#define _EIKONAL_H

#include<cblas.h>
#include<lapacke.h>
#include<nlopt.h>
#include<cstdlib>
#include<cmath>
#include<iostream>
#include<cstring>
#include<legQuad.hh>
#include<sFactors.hh>
#include<kMat.hh>
#include<dMat.hh>
#include<lMat.hh>
#include<random>

struct optDataEik
{

  unsigned int N;
  double *vecx;
  double *store;
  double *Dx, *Dy, *rhs;
  double *Ve;
  unsigned int ne;
  double *cu;
  double *Hess;
 
  optDataEik (  unsigned int _N,
                double* _Dx, double* _Dy,
                double* _Hess,
                double* _rhs, double* _cu,
                unsigned int _ne = 0,
                double* _Ve = nullptr )
    : N(_N), Dx(_Dx), Dy(_Dy), Hess(_Hess), rhs(_rhs), 
      Ve(_Ve), ne(_ne), cu(_cu) 
  {
    // Coeff derivs in x,y
    vecx = (double*) calloc(N, sizeof(double));
    // storage for contraint result
    store = (double*) calloc(ne, sizeof(double));
  }

  ~optDataEik ()
  {
    free(vecx); 
    free(store); 
  }

};

void transpose  ( const double* A, 
                  const unsigned int m, 
                  const unsigned int n,
                  double* AT )
{
  for (unsigned int j = 0; j < n; ++j)
  {
    for (unsigned int i = 0; i < m; ++i)
    {
      AT[j + n*i] = A[i + m*j];
    }
  }
}

struct eikonal
{
  unsigned int n, N, Ne, Np, nthreads;
  unsigned int nl, nb, nh;
  double a, b, c;

 
  jPoly<double>* polyabc    = 0;
  jPoly<double>* polya1b1c1 = 0;
  jPoly<double>* polyedge   = 0;
  legQuad<double>* edgel    = 0; 
  legQuad<double>* edgeb    = 0; 
  legQuad<double>* edgeh    = 0; 

  double *cu      = 0;
  double *cu_eq   = 0;
  double *crhs0   = 0, *crhs1   = 0;
  double *crhs2   = 0, *crhs    = 0;
  double *Xedge   = 0, *Yedge   = 0; 
  double *X       = 0, *Y       = 0;
  double *W       = 0;

  double *Habc    = 0, *Ha1bc   = 0;
  double *Ha1b1c  = 0, *Ha1b1c1 = 0;
  double *Ha1bc1  = 0, *Hab1c   = 0;
  double *Hab1c1  = 0;

  double* Kabc_a1bc     = 0;
  double* Ka1bc_a1b1c   = 0;
  double* Ka1b1c_a1b1c1 = 0;
  double* Ka1bc1_a1b1c1 = 0;
  double* Kab1c1_a1b1c1 = 0;

  double *Dx0 = 0, *Dy0 = 0;
  double *Dx  = 0, *Dy  = 0;
  double *Wx  = 0, *Wy  = 0;
  double *Lx  = 0, *Ly  = 0;
  
  double* G   = 0, *G0  = 0;
  double* G1  = 0;

  double *S     = 0, *UT      = 0;
  double *V     = 0, *superb  = 0; 
  double* nullP = 0, *AT      = 0;

  bool weighted = false;
  bool unconstrained = false;
   
  unsigned int Nopt;


  eikonal ( unsigned int _n,
            unsigned int _N,
            unsigned int _nl,
            unsigned int _nb,
            unsigned int _nh,
            double* frhs,
            double* fu,
            double* _X,
            double* _Y, 
            double* _W,
            unsigned int _nthreads,
            bool _weighted = false,
            bool _unconstrained = false)
    : n(_n), N(_N), Ne(_nl+_nb+_nh),
      nl(_nl), nb(_nb), nh(_nh),
      X(_X), Y(_Y), W(_W), 
      nthreads(_nthreads),
      weighted(_weighted),
      unconstrained(_unconstrained)
  {

    this->a             = 0.5; 
    this->b             = 0.5; 
    this->c             = 0.5;
    this->polyabc       = new jPoly<double>(N, n, a, b, c, nthreads); 
    this->Np            = polyabc->Np;
    if (Ne > Np)
    {
      std::cerr << "ERROR: Number of edge constraints must be less than\
                    number of polynomials" << std::endl;
      exit(1);
    }
    this->cu            = (double*) calloc(Np, sizeof(double));
    if (this->weighted)
    {
      this->polya1b1c1  = new jPoly<double>(N, n, a+1, b+1, c+1, nthreads);
      this->polya1b1c1->makeWeighted(); this->unconstrained = true;
      this->Lx          = (double*) calloc(Np*Np, sizeof(double));
      this->Ly          = (double*) calloc(Np*Np, sizeof(double));
      this->Wx          = (double*) calloc(Np*Np, sizeof(double));
      this->Wy          = (double*) calloc(Np*Np, sizeof(double));
      this->polya1b1c1->computeCoeffs(fu, X, Y, W, cu);

    }
    if (not this->weighted)
    {
      this->polyedge    = new jPoly<double>(Ne, n, a, b, c, nthreads);
      this->edgel       = new legQuad<double>(nl); this->edgel->shift();
      this->edgeb       = new legQuad<double>(nb); this->edgeb->shift();
      this->edgeh       = new legQuad<double>(nh); this->edgeh->shift();
    }
    this->crhs0         = (double*) calloc(Np, sizeof(double)); 
    this->crhs1         = (double*) calloc(Np, sizeof(double)); 
    this->crhs2         = (double*) calloc(Np, sizeof(double)); 
    this->crhs          = (double*) calloc(Np, sizeof(double)); 
    this->Xedge         = (double*) calloc(Ne, sizeof(double));   
    this->Yedge         = (double*) calloc(Ne, sizeof(double));   
    this->Habc          = (double*) calloc((n+2)*(n+2), sizeof(double));
    this->Ha1bc         = (double*) calloc((n+2)*(n+2), sizeof(double));
    this->Ha1b1c        = (double*) calloc((n+2)*(n+2), sizeof(double));
    this->Ha1b1c1       = (double*) calloc((n+2)*(n+2), sizeof(double));
    this->Ha1bc1        = (double*) calloc((n+2)*(n+2), sizeof(double));
    this->Hab1c         = (double*) calloc((n+2)*(n+2), sizeof(double));
    this->Hab1c1        = (double*) calloc((n+2)*(n+2), sizeof(double));

    if (not this->weighted)
    {
      this->Kabc_a1bc     = (double*) calloc(Np*Np, sizeof(double));    
      this->Ka1bc_a1b1c   = (double*) calloc(Np*Np, sizeof(double));    
      this->Ka1b1c_a1b1c1 = (double*) calloc(Np*Np, sizeof(double));    
      this->Ka1bc1_a1b1c1 = (double*) calloc(Np*Np, sizeof(double));    
      this->Kab1c1_a1b1c1 = (double*) calloc(Np*Np, sizeof(double));    
      this->Dx0           = (double*) calloc(Np*Np, sizeof(double));
      this->Dy0           = (double*) calloc(Np*Np, sizeof(double));
    }
    unsigned int dimS = Ne;
    this->polyabc->computeCoeffs(frhs, X, Y, W, crhs0);
    this->Dx            = (double*) calloc(Np*Np, sizeof(double));
    this->Dy            = (double*) calloc(Np*Np, sizeof(double));
    if (not this->unconstrained)
    {
      this->polyabc->computeCoeffs(fu, X, Y, W, cu);
    }


    sFactors(n+2, a, b, c, Habc);  
    sFactors(n+2, a+1, b, c, Ha1bc);  
    sFactors(n+2, a+1, b+1, c, Ha1b1c);  
    sFactors(n+2, a+1, b+1, c+1, Ha1b1c1);  
    sFactors(n+2, a+1, b, c+1, Ha1bc1);  
    sFactors(n+2, a, b+1, c, Hab1c);  
    sFactors(n+2, a, b+1, c+1, Hab1c1);  
   
    if (not this->weighted)
    { 
      kMat(a, b, c, Habc, Ha1bc, n, 0, Kabc_a1bc);
      kMat(a+1, b, c, Ha1bc, Ha1b1c, n, 1, Ka1bc_a1b1c);
      kMat(a+1, b+1, c, Ha1b1c, Ha1b1c1, n, 2, Ka1b1c_a1b1c1);
      kMat(a+1, b, c+1, Ha1bc1, Ha1b1c1, n, 1, Ka1bc1_a1b1c1);
      kMat(a, b+1, c+1, Hab1c1, Ha1b1c1, n, 0, Kab1c1_a1b1c1);
      dMat(a, b, c, Habc, Ha1bc1, n, 0, Dx0);
      dMat(a, b, c, Habc, Hab1c1, n, 1, Dy0); 
      cblas_dgemv ( CblasColMajor, CblasNoTrans,
                    Np, Np, 1.0, Kabc_a1bc, Np, crhs0, 1, 0.0, crhs1, 1);
      cblas_dgemv ( CblasColMajor, CblasNoTrans,
                    Np, Np, 1.0, Ka1bc_a1b1c, Np, crhs1, 1, 0.0, crhs2, 1);
      cblas_dgemv ( CblasColMajor, CblasNoTrans,
                    Np, Np, 1.0, Ka1b1c_a1b1c1, Np, crhs2, 1, 0.0, crhs, 1);
      
      cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans,
                    Np, Np, Np, 1.0, Ka1bc1_a1b1c1, Np, 
                    Dx0, Np, 0.0, Dx, Np);
      cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans,
                    Np, Np, Np, 1.0, Kab1c1_a1b1c1, Np, 
                    Dy0, Np, 0.0, Dy, Np);
      double *xl = edgel->x, *xb = edgeb->x, *xh = edgeh->x;
      for (unsigned int i = 0; i < nl; ++i)     { Yedge[i] = xl[i];           }
      for (unsigned int i = nl; i < nl+nb; ++i) { Xedge[i] = xb[i-nl];        } 
      for (unsigned int i = nl+nb; i < Ne; ++i) { Xedge[i] = xh[i-nl-nb]    ; }
      this->polyedge->computeV(Xedge, Yedge); 
      this->G0            = (double*) calloc(Np*Np, sizeof(double));
      cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans,
                    Np, Np, Np, 1.0, Dx, Np, 
                    Dx, Np, 0.0, G0, Np );
      cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans,
                    Np, Np, Np, 1.0, Dy, Np, 
                    Dy, Np, 1.0, G0, Np );
      if (this->unconstrained)
      {
        this->AT            = (double*) calloc(Np*Ne, sizeof(double));
        this->S             = (double*) calloc(dimS, sizeof(double));
        this->UT            = (double*) calloc(Ne*Ne, sizeof(double)); 
        this->V             = (double*) calloc(Np*Np, sizeof(double));
        this->superb        = (double*) calloc(dimS, sizeof(double));
        
        transpose(polyedge->V, Ne, Np, this->AT);
        LAPACKE_dgesvd( LAPACK_COL_MAJOR, 'A','A', 
                        Np, Ne, AT, Np, S, 
                        V, Np, UT, Ne, superb );
        
        this->Nopt          = Np - Ne;
        this->nullP         = &V[Np*Ne];
        this->G1            = (double*) calloc(Np*Nopt, sizeof(double));
        this->G             = (double*) calloc(Nopt*Nopt, sizeof(double));
        // nullP is Np x Nopt, G0 is Np x Np, G0*nullP is Np x Nopt
        cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans,
                      Np, Nopt, Np, 1.0, G0, Np, nullP, Np, 0.0, G1, Np );
        // nullPT is Nopt x Np, G1 is Np x Nopt, nullPT*G1 is Nopt x Nopt
        cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans,
                      Nopt, Nopt, Np, 1.0, nullP, Np, G1, Np, 0.0, G, Nopt );

        this->cu_eq = (double*) calloc(Nopt, sizeof(double));
        for (unsigned int i = 0; i < Nopt; ++i)
        {
          cu_eq[i] = std::pow(-1,i)*1.0;
        }
      }
      else
      {
        this->G             = G0;
        this->Nopt          = Np;
      }  
    }
    else if (this->weighted)
    {
      lMat(a, b+1, c, Hab1c, Habc, n, 1, Lx);
      lMat(a+1, b, c, Ha1bc, Habc, n, 0, Ly);
      dMat(a+1, b+1, c+1, Ha1b1c1, Hab1c, n, 0, Wx, weighted);
      dMat(a+1, b+1, c+1, Ha1b1c1, Ha1bc, n, 1, Wy, weighted);
      cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans,
                    Np, Np, Np, 1.0, Lx, Np, 
                    Wx, Np, 0.0, Dx, Np );
      cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans,
                    Np, Np, Np, 1.0, Ly, Np, 
                    Wy, Np, 0.0, Dy, Np );

      cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans,
                    Np, Np, Np, 1.0, Dx, Np, 
                    Dx, Np, 0.0, G, Np );
      cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans,
                    Np, Np, Np, 1.0, Dy, Np, 
                    Dy, Np, 1.0, G, Np );
    } 

 
  }


  void btransform(double* u)
  {
    if (not this->weighted)
    {
      if (this->unconstrained)
      {
        cblas_dgemv ( CblasColMajor, CblasNoTrans,
                      Np, Nopt, 1.0, nullP, Np, cu_eq, 1, 0.0, cu, 1);
      }
      cblas_dgemv ( CblasColMajor, CblasNoTrans,
                    N, Np, 1.0, polyabc->V, N, cu, 1, 0.0, u, 1);
    }
    else
    {
      cblas_dgemv ( CblasColMajor, CblasNoTrans,
                    N, Np, 1.0, polya1b1c1->V, N, cu, 1, 0.0, u, 1);
    }
  }
  
  ~eikonal()
  {
    if (polyabc)        { delete polyabc;       polyabc       = 0; }
    if (polya1b1c1)     { delete polya1b1c1;    polya1b1c1    = 0; }
    if (polyedge)       { delete polyedge;      polyedge      = 0; }
    if (edgel)          { delete edgel;         edgel         = 0; }
    if (edgeb)          { delete edgeb;         edgeb         = 0; }
    if (edgeh)          { delete edgeh;         edgeh         = 0; }
    if (cu)             { free(cu);             cu            = 0; }
    if (cu_eq)          { free(cu_eq);          cu_eq         = 0; }
    if (crhs)           { free(crhs);           crhs          = 0; }
    if (crhs0)          { free(crhs0);          crhs0         = 0; }
    if (crhs1)          { free(crhs1);          crhs1         = 0; }
    if (crhs2)          { free(crhs2);          crhs2         = 0; }
    if (Xedge)          { free(Xedge);          Xedge         = 0; }
    if (Yedge)          { free(Yedge);          Yedge         = 0; }
    if (Habc)           { free(Habc);           Habc          = 0; }
    if (Ha1bc)          { free(Ha1bc);          Ha1bc         = 0; }
    if (Ha1b1c1)        { free(Ha1b1c1);        Ha1b1c1       = 0; }
    if (Ha1bc1)         { free(Ha1bc1);         Ha1bc1        = 0; }
    if (Hab1c)          { free(Hab1c);          Hab1c         = 0; }
    if (Hab1c1)         { free(Hab1c1);         Hab1c1        = 0; }
    if (Ha1b1c)         { free(Ha1b1c);         Ha1b1c        = 0; }
    if (Kabc_a1bc)      { free(Kabc_a1bc);      Kabc_a1bc     = 0; }
    if (Ka1bc_a1b1c)    { free(Ka1bc_a1b1c);    Ka1bc_a1b1c   = 0; }
    if (Ka1b1c_a1b1c1)  { free(Ka1b1c_a1b1c1);  Ka1b1c_a1b1c1 = 0; }
    if (Ka1bc1_a1b1c1)  { free(Ka1bc1_a1b1c1);  Ka1bc1_a1b1c1 = 0; }
    if (Kab1c1_a1b1c1)  { free(Kab1c1_a1b1c1);  Kab1c1_a1b1c1 = 0; }
    if (Dx0)            { free(Dx0);            Dx0           = 0; }
    if (Dy0)            { free(Dx0);            Dy0           = 0; }
    if (Dx)             { free(Dx);             Dx            = 0; }
    if (Dy)             { free(Dy);             Dy            = 0; }
    if (Lx)             { free(Lx);             Lx            = 0; }
    if (Ly)             { free(Ly);             Ly            = 0; }
    if (Wx)             { free(Wx);             Wx            = 0; }
    if (Wy)             { free(Wy);             Wy            = 0; }
    if (G && 
        unconstrained)  { free(G);              G             = 0; }
    if (G0)             { free(G0);             G0            = 0; }
    if (G1)             { free(G1);             G1            = 0; }
    if (superb)         { free(superb);         superb        = 0; }
    if (UT)             { free(UT);             UT            = 0; }
    if (V)              { free(V);              V             = 0; }
    if (S)              { free(S);              S             = 0; }
    if (AT)             { free(AT);             AT            = 0; }
  } 
 
};



void pre(unsigned int n, const double* cu, const double* v, double* vpre, void* _data)
{
  optDataEik* data  = (optDataEik*) _data;
  unsigned int N    = data->N;
  double* Hess      = data->Hess;
  cblas_dgemv ( CblasColMajor, CblasNoTrans,
                N, N, 1.0, Hess, N, v, 1, 0.0, vpre, 1);

}

static unsigned long counteik = 0;
double Fhlp  ( const double* cu, void* _data )
{
  ++counteik;
  optDataEik* data  = (optDataEik*) _data;
  unsigned int N    = data->N;
  double* vecx      = data->vecx; 
  double* rhs       = data->rhs;
  double* Hess      = data->Hess;
 
  cblas_dgemv ( CblasColMajor, CblasNoTrans,
                N, N, 1.0, Hess, N, cu, 1, 0.0, vecx, 1);

  double eik    = cblas_ddot(N, cu, 1, vecx, 1);
  double rhssq  = cblas_ddot(N, rhs, 1, rhs, 1);
  double Fval   = eik - rhssq;
  Fval *= Fval;
  if ( !(counteik % 100) )
  { 
    std::cout << "Eval #" << counteik 
              << " : F = " << Fval
              << std::endl; 
  }
  
  return Fval;

}


                

double F  ( unsigned int n, const double* cu, 
            double* grad, void* _data  )
{

  if (grad)
  {
    optDataEik* data  = (optDataEik*) _data;
    double* Hess      = data->Hess;
    cblas_dgemv ( CblasColMajor,
                  CblasNoTrans,
                  n, n, 2.0, Hess, n,
                  cu, 1,
                  0, grad, 1 );
    
    double df = 2.0 * std::sqrt(Fhlp(cu, _data));
    for (unsigned int i = 0; i < n; ++i)
    {
      grad[i] *= df; 
    } 
  }

  return Fhlp (cu, _data);

}


void cl ( unsigned int m, double* result, unsigned int n, 
          const double* cu, double* grad, void* f_data)
{
  optDataEik* data   = (optDataEik*) f_data;
  double* Ve      = data->Ve;
  double* store   = data->store;
  
  if (grad)
  {
    for (unsigned int j = 0; j < n; ++j)
    {
      for (unsigned int i = 0; i < m; ++i)
      {
        grad[j + n*i] = Ve[i + m*j];
      }
    }
  }

  cblas_dgemv ( CblasColMajor, CblasNoTrans,  
                m, n, 1.0, Ve, m, cu, 1, 0.0, store, 1 );
  
  for (unsigned int i = 0; i < m; ++i)
  {
    result[i] = store[i];
  }
}


#endif
