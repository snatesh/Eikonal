#include<legQuad.hh>

extern "C"
{
   
  legQuad<double>* legendre (unsigned int n) 
  {
    legQuad<double>* lgq = new legQuad<double>(n);
    return lgq;
  } 
  
  void copyQuad  ( legQuad<double>* lgq,
                   double* legx, 
                   double* legw )
  {
    for (unsigned int i = 0; i < lgq->nlg; ++i)
    {
      legx[i] = lgq->x[i];
      legw[i] = lgq->w[i];
    }
  }

  void deleteQuad(legQuad<double>* leg)
  {
    delete leg;
  }

}
