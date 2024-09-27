#ifndef _JMAT_H
#define _JMAT_H

#include<iomanip>
#include<fstream>
#include<string>
#include<algorithm>
#include<sFactors.h>
#include<jWeight.h>
#include<jPoly.h>
#include<cblas.h>
#include<mapTensorQuad.h>

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

template <typename T> 
struct jMat
{
  unsigned int n, N, nlg; 
  unsigned int dim; 
  T *Jn1 = 0, *Jn2 = 0, *Jn3 = 0;
  T *H = 0, *A = 0, *B = 0, *C = 0;
  T *D = 0, *E = 0, *F = 0, *G = 0;
  T **A1 = 0, **A2 = 0, **B1 = 0, **B2 = 0;
  T **A3 = 0, **B3 = 0;
  T a, b, c, d, kap;
  T *x, *w; // 1d quad rule
  T *J, *avecON, *bvec;

  unsigned int blkind;

  jMat(unsigned int _n, T _a, T _b)
    : n(_n), a(_a), b(_b), dim(1)
  {
    J       = (T*) calloc(n*n, sizeof(T));
    bvec    = (T*) calloc(n, sizeof(T));
    avecON  = (T*) calloc(n, sizeof(T));
    T av, bv, asq = a * a, bsq = b * b;
    bvec[0]   = -(0.5 * (a - b)) / (0.5 * (a + b) + 1.0);
    avecON[0] = (2.0 / (a + b + 2.0)) * sqrt( (a + 1.0) * (b + 1.0) / (a + b + 3.0) );
    #pragma omp simd
    for (unsigned int i = 1; i < n; ++i)
    {
      av = 
        (2.0 * i + a + b + 1.0) * 
        (2.0 * i + a + b + 2.0) / 
        ( 2.0 * (i + 1.0) * (i + a + b + 1.0) );
      
      bv = 
        (asq - bsq) * (2.0 * i + a + b + 1.0) / 
        ( 2.0 * (i + 1.0) * (i + a + b + 1.0) * (2.0 * i + a + b) );
      
      bvec[i] = -bv / av;
      
      avecON[i] = 
        2 / (a + b + 2.0 * i + 2.0) * 
        sqrt(
              (a + i + 1.0) * (b + i + 1.0) * 
              (i + 1.0) * (a + b + i + 1.0) /
              ( (a + b + 2.0 * i + 1.0) * 
                (a + b + 2.0 * i + 3.0) )  
            );
    }

    #pragma omp simd
    for (unsigned int i = 0; i < n; ++i)
    {
      J[i + n*i] = bvec[i];

    }

    #pragma omp simd
    for (unsigned int i = 0; i < n-1; ++i)
    {
      J[i + n*(i+1)] = avecON[i];
      J[i+1 + n*i] = avecON[i];
    }
   
  } 
    
  jMat (  unsigned int _n, 
          T _a, T _b, T _c , T _d, 
          unsigned int _nlg,
          T* _x, T* _w,
          unsigned int nthreads )
    : n(_n), a(_a), b(_b), c(_c), d(_d), 
      nlg(_nlg), x(_x), w(_w), dim(3)
  {
    if( a != 0.5 || b != 0.5 || c != 0.5 || d != 0.5 )
    {
      std::cerr << "ERROR: only legendre analog polys are supported\n"
                << "       for numerical evaluation of Jacobi matrices with d>2";
      exit(1);
    }
    N = dimPI3(n-1);
    Jn1 = (T*) calloc(N*N, sizeof(T));
    Jn2 = (T*) calloc(N*N, sizeof(T));
    Jn3 = (T*) calloc(N*N, sizeof(T));
    /* block array of coefficient matrices in 3-term recurrence*/
    A1 = (T**) malloc((n+1) * sizeof(T*));
    A2 = (T**) malloc((n+1) * sizeof(T*));
    A3 = (T**) malloc((n+1) * sizeof(T*));
    B1 = (T**) malloc((n+1) * sizeof(T*));
    B2 = (T**) malloc((n+1) * sizeof(T*));
    B3 = (T**) malloc((n+1) * sizeof(T*));
    unsigned int M, NN;
    for (unsigned int nn = 0; nn <= (n-1); ++nn)
    {
      M = rn3(nn); NN = rn3(nn+1);
      A1[nn]    = (T*) calloc(M*NN, sizeof(T));
      A2[nn]    = (T*) calloc(M*NN, sizeof(T));
      A3[nn]    = (T*) calloc(M*NN, sizeof(T));
      B1[nn]    = (T*) calloc(M*M, sizeof(T));
      B2[nn]    = (T*) calloc(M*M, sizeof(T));
      B3[nn]    = (T*) calloc(M*M, sizeof(T));
    }

    // mapped quad rule on tet
    mapTensorQuad<T>* C2T = new mapTensorQuad<T>(nlg, x, w);
    // evaluate vandermonde on abscissa from mapped rule on tet
    jPoly<T>* Pn1 = new jPoly<T>(nlg*nlg*nlg, n+2, a, b, c, d, nthreads); 
    Pn1->computeV(C2T->X, C2T->Y, C2T->Z);
    // get weight function evaluator for tet
    T params[4] = {a, b, c, d}; 
    jWeight<T,3>* Jw = new jWeight<T,3>(params, x, w);
    // accumulate inner products into into blocks 
    blkind = 0;
    for (unsigned int nn = 0; nn <= (n-1); ++nn)
    {
      jBlock3(nn, C2T, Pn1, Jw, A1[nn], A2[nn], A3[nn], B1[nn], B2[nn], B3[nn]); 
      blkind += rn3(nn);
    }
    delete Jw;
    delete Pn1; 
    delete C2T;
    
    unsigned int i, inds1[2], inds[2] = {1, 1}, bsz = 1;
    for (unsigned int nn = 1; nn <= (n-1); ++nn)
    {
      i = 0;
      for (unsigned int col = inds[0]; col <= inds[1]; ++col)
      {
        for (unsigned int row = inds[0]; row <= inds[1]; ++row)
        {
          Jn1[(row-1) + N*(col-1)] = B1[nn-1][i]; 
          Jn2[(row-1) + N*(col-1)] = B2[nn-1][i];
          Jn3[(row-1) + N*(col-1)] = B3[nn-1][i];
          i += 1; 
        }
      }
      bsz += nn + 1;
      inds1[0] = inds[0] + bsz - (nn + 1);
      inds1[1] = inds[1] + bsz;      
      i = 0;
      for (unsigned int col = inds1[0]; col <= inds1[1]; ++col)
      {
        for (unsigned int row = inds[0]; row <= inds[1]; ++row)
        {
          Jn1[(row-1) + N*(col-1)] = A1[nn-1][i];
          Jn1[(col-1) + N*(row-1)] = A1[nn-1][i];
          Jn2[(row-1) + N*(col-1)] = A2[nn-1][i];
          Jn2[(col-1) + N*(row-1)] = A2[nn-1][i];
          Jn3[(row-1) + N*(col-1)] = A3[nn-1][i];
          Jn3[(col-1) + N*(row-1)] = A3[nn-1][i];
          i += 1;
        }
      }
      inds[0] = inds[0] + bsz - (nn + 1);
      inds[1] = inds[1] + bsz;
    }
    // get last block
    i = 0;
    for (unsigned int col = inds[0]; col <= inds[1]; ++col)
    {
      for (unsigned int row = inds[0]; row <= inds[1]; ++row)
      {
        Jn1[(row-1) + N*(col-1)] = B1[n-1][i]; 
        Jn2[(row-1) + N*(col-1)] = B2[n-1][i];
        Jn3[(row-1) + N*(col-1)] = B3[n-1][i];
        i += 1;
      }
    }
  }

  void jBlock3( unsigned int n,
                mapTensorQuad<T>* C2T,
                jPoly<T>* Pn1,
                jWeight<T,3>* Jw,
                T* Ax, T* Ay, T* Az,
                T* Bx, T* By, T* Bz )
  {
    T *W = C2T->W;
    const T *X = Pn1->X, *Y = Pn1->Y, *Z = Pn1->Z, *V = Pn1->V;
    T alphax, alphay, alphaz, wval, Xt[3];
    unsigned int M = rn3(n), NN = rn3(n+1);
    unsigned int coln = dimPI3(n), colnp1 = dimPI3(n+1);
    std::cout << coln << " " << colnp1 << std::endl;
    std::cout << M << " " << NN << std::endl;
    unsigned int nlg3 = nlg *nlg *nlg;
    unsigned int rn = rn3(n); 
    T* Pn   = (T*) calloc(M, sizeof(T));
    T* Pnp1 = (T*) calloc(NN, sizeof(T));


    for (unsigned int i = 0; i < nlg3; ++i)
    {
      for (unsigned int iblk = 0; iblk < M; ++iblk)
      { 
        Pn[iblk] = V[i + nlg3*(blkind+iblk)]; 
      }
      for (unsigned int iblk = 0; iblk < NN; ++iblk)
      { 
        Pnp1[iblk] = V[i + nlg3*(blkind+rn+iblk)]; 
      }
      Xt[0] = X[i]; Xt[1] = Y[i]; Xt[2] = Z[i];
      wval = Jw->w(Xt);
      alphax = wval * Xt[0] * W[i];
      alphay = wval * Xt[1] * W[i];
      alphaz = wval * Xt[2] * W[i];
      cblas_dger(CblasColMajor, M, NN, alphax, Pn, 1, Pnp1, 1, Ax, M); 
      cblas_dger(CblasColMajor, M, NN, alphay, Pn, 1, Pnp1, 1, Ay, M); 
      cblas_dger(CblasColMajor, M, NN, alphaz, Pn, 1, Pnp1, 1, Az, M); 
      cblas_dger(CblasColMajor, M, M, alphax, Pn, 1, Pn, 1, Bx, M); 
      cblas_dger(CblasColMajor, M, M, alphay, Pn, 1, Pn, 1, By, M); 
      cblas_dger(CblasColMajor, M, M, alphaz, Pn, 1, Pn, 1, Bz, M); 
    }
    free(Pn);
    free(Pnp1);
  }

  jMat (unsigned int _n, T _a, T _b, T _c, T _d, std::string dir)
    : n(_n), a(_a), b(_b), c(_c), d(_d), dim(3)
  {
    N = dimPI3(n-1);
    Jn1 = (T*) calloc(N*N, sizeof(T));
    Jn2 = (T*) calloc(N*N, sizeof(T));
    Jn3 = (T*) calloc(N*N, sizeof(T));
    std::stringstream ssx, ssy, ssz; 
    ssx << dir << "J" << n << "x_tet.txt";
    ssy << dir << "J" << n << "z_tet.txt";
    ssz << dir << "J" << n << "y_tet.txt";
    std::ifstream Jxfile(ssx.str());
    std::ifstream Jyfile(ssy.str());
    std::ifstream Jzfile(ssz.str());
    for (unsigned int j = 0; j < N*N; ++j)
    {
      Jxfile >> Jn1[j];
      Jyfile >> Jn2[j];
      Jzfile >> Jn3[j];
    }
    Jxfile.close(); Jyfile.close(); Jzfile.close();
    ssx.str() = ""; ssy.str() = ""; ssz.str() = "";
  }
    

  jMat ( unsigned int _n, T _a, T _b, T _c )
    : n(_n), a(_a), b(_b), c(_c), dim(2)
  {
    N = static_cast<unsigned int>(0.5 * n * (n + 1)); 
    kap = abs(a + b + c);

    // normalization constants
    H = (T*) calloc((n+3)*(n+3), sizeof(T));
    sFactors<T>(n+3, a, b, c, H);
    /* temporary upper triangular matrices for
       constructing coefficient matrices in recurrence */
    A = (T*) calloc((n+1)*(n+1), sizeof(T)); 
    B = (T*) calloc((n+1)*(n+1), sizeof(T)); 
    C = (T*) calloc((n+1)*(n+1), sizeof(T)); 
    D = (T*) calloc((n+1)*(n+1), sizeof(T)); 
    E = (T*) calloc((n+1)*(n+1), sizeof(T)); 
    F = (T*) calloc((n+1)*(n+1), sizeof(T)); 
    G = (T*) calloc((n+1)*(n+1), sizeof(T)); 
    /* block array of coefficient matrices in 3-term recurrence*/
    A1 = (T**) malloc((n+1) * sizeof(T*));
    A2 = (T**) malloc((n+1) * sizeof(T*));
    B1 = (T**) malloc((n+1) * sizeof(T*));
    B2 = (T**) malloc((n+1) * sizeof(T*));
    Jn1 = (T*) calloc(N*N, sizeof(T));
    Jn2 = (T*) calloc(N*N, sizeof(T));

    // populated tmp up tri mats
    for (unsigned int nn = 0; nn <= n; ++nn)
    {
      for (unsigned int kk = 0; kk <= nn; ++kk)
      {
        A[kk + nn*(n+1)] = 
          (H[kk + (nn+1)*(n+3)] / H[kk + nn*(n+3)]) * 
          (nn - kk + 1.0) * (nn + kk + kap + 0.5) /
          ( (2.0 * nn + kap + 0.5) * (2.0 * nn + kap + 1.5) );
         
        B[kk + nn*(n+1)] = 
          0.5 - ( pow(b + c + 2.0 * kk, 2.0) - pow(a - 0.5, 2.0) ) /
          ( 2.0 * (2.0 * nn + kap - 0.5) * (2.0 * nn + kap + 1.5) );
        
        if (nn > 0)
        {
          for (unsigned int kk1 = 1; kk1 <= nn; ++kk1)
          {
            C[kk1-1 + nn*(n+1)] = 
              ( (H[kk1-1 + (nn+1)*(n+3)] / H[kk1 + nn*(n+3)]) * (nn - kk1 + 1.0) * 
                (nn - kk1 + 2.0) * (kk1 + b - 0.5) * (kk1 + c - 0.5) ) /
              ( (2.0 * nn + kap + 0.5) * (2.0 * nn + kap + 1.5) * 
                (2.0 * kk1 + b + c) * (2.0 * kk1 + b + c - 1.0) );          
          }
        }
        
        E[kk + nn*(n+1)] = 
          -( 1.0 + ( pow(b - 0.5, 2.0) - pow(c - 0.5, 2.0) ) /
            ( (2.0 * kk + b + c + 1.0) * 
              (2.0 * kk + b + c - 1.0) ) ) * 
          A[kk + nn*(n+1)] / 2.0;

        F[kk + nn*(n+1)] = 
          ( 1.0 + ( pow(b - 0.5, 2.0) - pow(c - 0.5, 2.0) ) /
            ( (2.0 * kk + b + c + 1.0) * 
              (2.0 * kk + b + c - 1.0) ) ) * 
          (1.0 - B[kk + nn*(n+1)]) / 2.0;  
       
        if (a == 0.5 && b == 0.5 && c == 0.5)
        {
          E[nn*(n+1)] = -A[nn*(n+1)] / 2.0;
          F[nn*(n+1)] = (1.0 - B[nn*(n+1)]) / 2.0;
        } 
        
        D[kk + nn*(n+1)] = 
          (H[kk+1 + (nn+1)*(n+3)] / H[kk + nn*(n+3)]) * 
          (nn + kk + kap + 0.5) * (nn + kk + kap + 1.5) * 
          (kk + 1.0) * (kk + b + c) /
          ( (2.0 * nn + kap + 0.5) * 
            (2.0 * nn + kap + 1.5) * 
            (2.0 * kk + b + c) * 
            (2.0 * kk + b + c + 1.0) 
          );

        G[kk + nn*(n+1)] = 
          (-2.0 * H[kk+1 + nn*(n+3)] / (H[kk + nn*(n+3)]) ) *
          (nn - kk + a - 0.5) * (nn + kk + kap + 0.5) * 
          (kk + 1) * (kk + b + c) /
          ( (2.0 * nn + kap - 0.5) * 
            (2.0 * nn + kap + 1.5) * 
            (2.0 * kk + b + c) * 
            (2.0 * kk + b + c + 1.0)  
          );
      }
    }
    // prepoulate blocks for easy assignment into Jn_i
    // efficiency/redundancy doesn't matter 
    // - amortize out precompute, go for readablity instead
    for (unsigned int nn = 0; nn <= (n-1); ++nn)
    {
      A1[nn]    = (T*) calloc((nn+1)*(nn+2), sizeof(T));
      A2[nn]    = (T*) calloc((nn+1)*(nn+2), sizeof(T));
      B1[nn]    = (T*) calloc((nn+1)*(nn+1), sizeof(T));
      B2[nn]    = (T*) calloc((nn+1)*(nn+1), sizeof(T));
      for (unsigned int ii = 0; ii <= nn; ++ii) 
      { 
        A1[nn][ii + ii*(nn+1)]        = A[ii + nn*(n+1)];
        B1[nn][ii + ii*(nn+1)]        = B[ii + nn*(n+1)];
        A2[nn][ii + ii*(nn+1)]        = E[ii + nn*(n+1)];
        A2[nn][ii + (ii+1)*(nn+1)]    = D[ii + nn*(n+1)];
        B2[nn][ii + ii*(nn+1)]        = F[ii + nn*(n+1)];
        if (ii < nn)
        {
          A2[nn][ii+1 + ii*(nn+1)]      = C[ii + nn*(n+1)];
          B2[nn][ii + (ii+1)*(nn+1)]    = G[ii + nn*(n+1)];
          B2[nn][ii+1 + ii*(nn+1)]      = G[ii + nn*(n+1)];
        }
      }
    }
    // copy into Jn_i
    unsigned int inds[2] = {1,1};
    unsigned int inds1[2];
    unsigned int i;
    for (unsigned int nn = 1; nn <= (n-1); ++nn)
    {
      i = 0;
      for (unsigned int col = inds[0]; col <= inds[1]; ++col)
      {
        for (unsigned int row = inds[0]; row <= inds[1]; ++row)
        {
          Jn1[(row-1) + N*(col-1)] = B1[nn-1][i]; 
          Jn2[(row-1) + N*(col-1)] = B2[nn-1][i];
          i += 1; 
        }
      }
      inds1[0] = inds[0] + nn;
      inds1[1] = inds[1] + nn + 1;
      i = 0;
      for (unsigned int col = inds1[0]; col <= inds1[1]; ++col)
      {
        for (unsigned int row = inds[0]; row <= inds[1]; ++row)
        {
          Jn1[(row-1) + N*(col-1)] = A1[nn-1][i];
          Jn1[(col-1) + N*(row-1)] = A1[nn-1][i];
          Jn2[(row-1) + N*(col-1)] = A2[nn-1][i];
          Jn2[(col-1) + N*(row-1)] = A2[nn-1][i];
          
          i += 1;
        }
      }
      inds[0] += nn;
      inds[1] += nn + 1;
    }
    // get last block
    i = 0;
    for (unsigned int col = inds[0]; col <= inds[1]; ++col)
    {
      for (unsigned int row = inds[0]; row <= inds[1]; ++row)
      {
        Jn1[(row-1) + N*(col-1)] = B1[n-1][i]; 
        Jn2[(row-1) + N*(col-1)] = B2[n-1][i];
        i += 1;
      }
    }
  }
  




  ~jMat()
  {
    if (this->dim == 1)
    {
      if (J) { free(J); J = 0; }
      if (avecON) { free(avecON); avecON = 0; }
      if (bvec) { free(bvec); bvec = 0; }
    }
    if (this->dim == 2)
    { 
      for (unsigned int nn = 0; nn < n; ++nn)
      {
        free(A1[nn]); free(A2[nn]);
        free(B1[nn]); free(B2[nn]);
      }
      free(A1); free(B1); 
      free(A2); free(B2); 
      free(H); free(A); 
      free(B); free(C); 
      free(D); free(E); 
      free(F); free(G);
      free(Jn1); free(Jn2);
    }
    else if (this->dim == 3)
    {
      free(Jn1); free(Jn2); free(Jn3);
      for (unsigned int nn = 0; nn < n; ++nn)
      {
        free(A1[nn]); free(A2[nn]);
        free(B1[nn]); free(B2[nn]);
        free(A3[nn]); free(B3[nn]);
      }
      free(A1); free(A2); free(A3);
      free(B1); free(B2); free(B3);
    } 
  }
}; 

template struct jMat<double>;
   
#endif
