#ifndef _JCODEC_H
#define _JCODEC_H

#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkSmartPointer.h>
#include <vtkImageInterpolator.h>
#include <vtkIncrementalOctreePointLocator.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>
#include <fstream>
#include <iostream>
#include <string>

#include<cblas.h>
#include<lapacke.h>
#include<dMat.hh>
#include<jPoly.hh>


vtkSmartPointer<vtkDelaunay2D> triangulateUniform ( int* dims, 
                                                    double* origin,
                                                    unsigned int nSamp );


void readQuad ( const std::string& trix, 
                const std::string& triy, 
                const std::string& triw,
                unsigned int N,
                double* X, double* Y, 
                double*  W  );

struct Triangulator
{
  unsigned int N, m, M, nthreads;
  unsigned int nSamp, nRuns;
  double a, b, c;
  jPoly<double> *Pm = 0, *Pmx = 0, *Pmy = 0;
  double *Habc = 0, *Ha1bc1 = 0, *Hab1c1 = 0;
  double *Dx = 0, *Dy = 0; 
  double *X = 0, *Y = 0, *W = 0; 
  double *cimg = 0, *cdimg = 0, *dimgr = 0, *dimgt = 0;
  double *interpc = 0; 
  vtkSmartPointer<vtkImageData> imagedata;
  vtkSmartPointer<vtkImageInterpolator> interpolator; 
  vtkSmartPointer<vtkPolyData> polytri1, polytri2, polytri3; 
  
  std::string trix, triy, triw; 
  
  Triangulator  ( unsigned int _N, unsigned int _m,
                  double _a, double _b, double _c,
                  unsigned int _nSamp, unsigned int _nRuns,
                  const std::string& _trix, 
                  const std::string& _triy,
                  const std::string& _triw,
                  vtkSmartPointer<vtkImageData> _imagedata,
                  vtkSmartPointer<vtkImageInterpolator> _interpolator, 
                  unsigned int _nthreads  );
  
  vtkSmartPointer<vtkDelaunay2D> triangulateEntropy ( double* intgn, unsigned int channel  );
                                                
  void run();    
  
  ~Triangulator();
};

struct Compressor
{
  // we can get coeffs up to/including the 21st poly subspace
  // since integral(p21 * f) is exact for f of deg <= 21 with 
  // quad order m=42
  //unsigned int N = 325, mmax = 21, nthreads, Mmax;
  unsigned int N = 325, mmax = 10, nthreads, Mmax;
  std::string trix = "xtri_N325_n24_M946_m42.txt";
  std::string triy = "ytri_N325_n24_M946_m42.txt";
  std::string triw = "wtri_N325_n24_M946_m42.txt";
  double *R = 0, *S = 0, *W = 0, *cimg = 0;
  double *interpc = 0;
  jPoly<double>* Pm;
  double a, b, c;
  vtkSmartPointer<vtkImageInterpolator> interpolator; 
  vtkSmartPointer<vtkPolyData> polytri1, polytri2, polytri3; 
  vtkSmartPointer<vtkDoubleArray> coeffs1, coeffs2, coeffs3;
  vtkSmartPointer<vtkIntArray> offsets1, offsets2, offsets3;
  vtkSmartPointer<vtkUnsignedIntArray> orders1, orders2, orders3;
 
  Compressor  ( Triangulator* T );

  void compressChannel  ( unsigned int channel, double blknormtol );
 
  void run  ( double ctol );

  ~Compressor();

};

struct Decompressor
{
  vtkSmartPointer<vtkPolyData> polytri1, polytri2, polytri3; 
  vtkSmartPointer<vtkDoubleArray> coeffs1, coeffs2, coeffs3;
  vtkSmartPointer<vtkIntArray> offsets1, offsets2, offsets3;
  vtkSmartPointer<vtkUnsignedIntArray> orders1, orders2, orders3;
  vtkSmartPointer<vtkImageData> imagedata;
  vtkSmartPointer<vtkUnsignedShortArray> colors;
  vtkSmartPointer<vtkPoints> pixels;
  double a, b, c;
  //unsigned int mmax = 21, Mmax = 253, nthreads;
  //unsigned int mmax = 5, Mmax = 21, nthreads;
  //unsigned int mmax = 15, Mmax = 136, nthreads;
  unsigned int mmax = 10, Mmax = 66, nthreads;
  double bounds[6];


  Decompressor  ( const char* channel1, 
                  const char* channel2, 
                  const char* channel3 );

  void decompressChannel  ( unsigned int channel  );

  void writeImage (const char* fname);

  void run ( unsigned int nthreads );
  
  ~Decompressor(){}

};



#endif
