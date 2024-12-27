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
  unsigned int N = 325, mmax = 21, nthreads, Mmax;
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

/*  
     1) will need to create an empty imagedata with same dims
        as original - gives pixel points.
     2) create map between cell index (triangle) and 
        the indices of pixel points in that triangle
     3) use enconding we defined in vtp files to interpolate
        using jacobi interp mat evauated on pixel points
        in that triangle with coeffs in that tri 
        (this uses the coeffs field array and offsets cell
         array we stored in the vtp during compression)
     4) This is to be stored as a vtkInt32/16/8 Array in vtkPointData
        with each tuple having 3 components (1 for each channel)
     5) Finally, we write it back to some lossless format like png
        (NOTE: we don't want to write to a lossy compressed format) 
*/
struct Decompressor
{
  vtkSmartPointer<vtkPolyData> polytri1, polytri2, polytri3; 
  vtkSmartPointer<vtkDoubleArray> coeffs1, coeffs2, coeffs3;
  vtkSmartPointer<vtkIntArray> offsets1, offsets2, offsets3;
  vtkSmartPointer<vtkUnsignedIntArray> orders1, orders2, orders3;
  vtkSmartPointer<vtkImageData> imagedata;
  vtkSmartPointer<vtkUnsignedCharArray> colors;
  vtkSmartPointer<vtkPoints> pixels;
  double a, b, c;
  unsigned int mmax = 21, Mmax = 253;
  double bounds[6];


  Decompressor  ( const char* channel1, 
                  const char* channel2, 
                  const char* channel3 );

  void decompressChannel  ( unsigned int channel  );

  void writeImage (const char* fname);

  void run ();
  
  ~Decompressor(){}

};



#endif
