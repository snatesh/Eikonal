#ifndef _JMAT_H
#define _JMAT_H
#include<sFactors.h>
#include<iomanip>
#include<fstream>
#include<string>

template <typename T> 
struct jMat
{
  unsigned int n, N; 
  unsigned int dim; 
  T *Jn1 = 0, *Jn2 = 0, *Jn3 = 0;
  T *H = 0, *A = 0, *B = 0, *C = 0;
  T *D = 0, *E = 0, *F = 0, *G = 0;
  T **A1 = 0, **A2 = 0, **B1 = 0, **B2 = 0;
  T a, b, c, d, kap;

  T *J, *avecON, *bvec;

  jMat(unsigned int _n, T _a, T _b)
    : n(_n), a(_a), b(_b)
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
    A1 = (T**) malloc((n+1) * sizeof(T**));
    A2 = (T**) malloc((n+1) * sizeof(T**));
    B1 = (T**) malloc((n+1) * sizeof(T**));
    B2 = (T**) malloc((n+1) * sizeof(T**));
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
    } 
  }
}; 

template struct jMat<double>;
template struct jMat<float>;
   
#endif
