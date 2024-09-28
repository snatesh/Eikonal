#include<jevd.hh>
#include<jMat.hh>
#include<legQuad.hh>
#include<fstream>
#include<cstdio>
//void printMat(const double* A, const unsigned int m, const unsigned int n)
//{
//  for (unsigned int i = 0; i < m; ++i)
//  {
//    for (unsigned int j = 0; j < n; ++j)
//    {
//      std::cout << std::setw(10);
//      std::cout << A[i + m*j] << " ";
//    }
//    std::cout << std::endl;
//  }
//  std::cout << std::endl;
//}

int main(int argc, char* argv[])
{


  double a = 0.5, b = 0.5, c = 0.5, d = 0.5; 
  unsigned int norder = std::stoi(argv[1]);
  unsigned int nlg = std::stoi(argv[2]);
  
  unsigned int N = dimPI3(norder-1);

  double* J = (double*) calloc(N*N*3, sizeof(double));
  double* X0 = (double*) calloc(N, sizeof(double)); 
  double* Y0 = (double*) calloc(N, sizeof(double)); 
  double* Z0 = (double*) calloc(N, sizeof(double)); 

  legQuad<double>* legq = new legQuad<double>(nlg); 
  
  jMat<double>* jmat = 
    new jMat<double>(norder, a, b, c, d, nlg, legq->x, legq->w, 1);

  //printMat(jmat->Jn1, N, N); 
  //printMat(jmat->Jn2, N, N); 
  //printMat(jmat->Jn3, N, N); 

  mapTensorQuad<double>* C2T = new mapTensorQuad<double>(nlg, legq->x, legq->w);
  
  std::ofstream XXfile("Xtet.txt");
  std::ofstream YYfile("Ytet.txt");
  std::ofstream ZZfile("Ztet.txt");
  std::ofstream WWfile("Wtet.txt");
  for (unsigned int i = 0; i < nlg*nlg*nlg; ++i)
  {
    XXfile << C2T->X[i] << std::endl;
    YYfile << C2T->Y[i] << std::endl;
    ZZfile << C2T->Z[i] << std::endl;
    WWfile << C2T->W[i] << std::endl;
  }  
  
  XXfile.close();
  YYfile.close();
  ZZfile.close();
  

  for (unsigned int j = 0; j < N; ++j)
  {
    for (unsigned int i = 0; i < N; ++i)
    {
      J[i + N*j]          = jmat->Jn1[i + j*N];
      J[i + N*(j + N)]    = jmat->Jn2[i + j*N];
      J[i + N*(j + 2*N)]  = jmat->Jn3[i + j*N];
    }
  } 

  jointDiag<double>* jevd = new jointDiag(N, 3, 1e-10, J, 1); 
  std::ofstream Xfile("xtet_new.txt");
  std::ofstream Yfile("ytet_new.txt");
  std::ofstream Zfile("ztet_new.txt");

  for (unsigned int i = 0; i < N; ++i) 
  { 
    Xfile << J[i + N*i] << std::endl; 
    Yfile << J[i + N*(i + N)] << std::endl;
    Zfile << J[i + N*(i + 2*N)] << std::endl; 
  } 

  Xfile.close();
  Yfile.close();
  Zfile.close();
 
  free(J);
  free(X0);
  free(Y0);
  free(Z0);
}
