#ifndef _JCODEC_H
#define _JCODEC_H

#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkSmartPointer.h>
#include <vtkImageInterpolator.h>
#include <vtkCellTreeLocator.h>
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
#include <legQuad.hh>


vtkSmartPointer<vtkDelaunay2D> triangulateUniform ( int* dims, 
                                                    double* origin,
                                                    unsigned int nSamp,
                                                    vtkIdType& ll, int& lr,
                                                    vtkIdType& ur, int& ul );


void readQuad ( const std::string& trix, 
                const std::string& triy, 
                const std::string& triw,
                unsigned int N, double* X, 
                double* Y, double*  W  );

/* Triangulator for compression */
struct Triangulator
{
  unsigned int N, m, M;
  unsigned int Mtarget;
  unsigned int nSamp, nRuns;
  int dims[3], Npix;
  double a, b, c;
  jPoly<double> *Pm = 0, *Pmx = 0, *Pmy = 0;
  double *Habc = 0, *Ha1bc1 = 0, *Hab1c1 = 0;
  double *Dx = 0, *Dy = 0; 
  double *X = 0, *Y = 0, *W = 0; 
  double *cimg = 0, *cdimg = 0, *dimgr = 0, *dimgt = 0;
  double *interpc = 0, *interpc_jacobi = 0;
  bool useMultiChannel = false;
  double bpp_target;
  bool* torefine;
  vtkIdType ll, lr, ur, ul;
 
  vtkSmartPointer<vtkImageData> imagedata;
  vtkSmartPointer<vtkImageInterpolator> interpolator; 
  vtkSmartPointer<vtkPolyData> polytri, polytri1, polytri2, polytri3; 
  vtkSmartPointer<vtkPolyData> polyBbox;
  
  std::string trix, triy, triw; 
  
  Triangulator  ( unsigned int _N, unsigned int _m,
                  double _a, double _b, double _c,
                  unsigned int _nSamp, unsigned int _nRuns,
                  const std::string& _trix, 
                  const std::string& _triy,
                  const std::string& _triw,
                  vtkSmartPointer<vtkImageData> _imagedata,
                  vtkSmartPointer<vtkImageInterpolator> _interpolator,
                  bool useMultiChannel );
  
  Triangulator  ( unsigned int _N, unsigned int _m,
                  double _a, double _b, double _c,
                  unsigned int _nSamp, double bpp_target,
                  const std::string& _trix, 
                  const std::string& _triy,
                  const std::string& _triw,
                  vtkSmartPointer<vtkImageData> _imagedata,
                  vtkSmartPointer<vtkImageInterpolator> _interpolator );

 
  double getBPP ( ); 
  
  vtkSmartPointer<vtkDelaunay2D> triangulateEntropyGreedy ( double* interr );

  vtkSmartPointer<vtkDelaunay2D> triangulateEntropy ( double* intgn, unsigned int channel  );
  double triangulateEntropy_help  ( double* intgn, unsigned int channel,
                                    vtkSmartPointer<vtkPolyData> polytri, double& stdev );
  double triangulateEntropyNoGrad_help  ( double* intgn, unsigned int channel,
                                          vtkSmartPointer<vtkPolyData> polytri,
                                          double& stdev);
  void triangulateEntropyGreedyL1J_help  ( double* interr, vtkSmartPointer<vtkPolyData> polytri );
  void triangulateEntropyGreedyL2J_help  ( double* interr, vtkSmartPointer<vtkPolyData> polytri );
  void triangulateEntropyGreedyL1P_help  ( double* interr, vtkSmartPointer<vtkPolyData> polytri );
                                                
  void run  ( );    
  
  ~Triangulator();
};

/* Compression after triangulation */
struct Compressor
{
  // we can get coeffs up to/including the 21st poly subspace
  // since integral(p21 * f) is exact for f of deg <= 21 with 
  // quad order m=42
  unsigned int Mmax, morder;
  unsigned int N = 496;  
  std::string trix = "xtri_N496_n30_M1378_m51.txt";
  std::string triy = "ytri_N496_n30_M1378_m51.txt";
  std::string triw = "wtri_N496_n30_M1378_m51.txt";
  double *R = 0, *S = 0, *W = 0, *cimg = 0;
  double *interpc = 0;
  bool useMultiChannel;
  jPoly<double>* Pm;
  double a, b, c;
  double totalBytes;
  vtkSmartPointer<vtkImageInterpolator> interpolator; 
  vtkSmartPointer<vtkPolyData> polytri, polytri1, polytri2, polytri3; 
  vtkSmartPointer<vtkDoubleArray> coeffs, coeffs1, coeffs2, coeffs3;
  vtkSmartPointer<vtkIntArray> offsets, offsets1, offsets2, offsets3, celltypes;
  vtkSmartPointer<vtkUnsignedIntArray> orders, orders1, orders2, orders3;
  legQuad<double>* legq = 0; 
  double ave0, ave1, ave2;
  double stdev0, stdev1, stdev2;
 
  Compressor  ( Triangulator* T , unsigned int order );

  void compressChannel  ( unsigned int channel );
  void compressChannel_help ( unsigned int channel,
                              vtkSmartPointer<vtkPolyData> polytri,
                              vtkSmartPointer<vtkDoubleArray> coeffs,
                              vtkSmartPointer<vtkIntArray> offsets,
                              vtkSmartPointer<vtkUnsignedIntArray> orders,
                              double& ave, double& stdev);
  void smoothCoeffs ( );
  void pruneCoeffs ( );
  void smoothCoeffs_alt ( );
  double smoothCoeffs_help  ( unsigned int channel,
                              unsigned int icell,
                              unsigned int nleg,
                              jPoly<double>* lPm,
                              jPoly<double>* bPm,
                              jPoly<double>* hPm,
                              double* cimg, double* cimg1,
                              double* imgbnd, double* imgbnd1,
                              double* pcoords10, double* pcoords11,
                              int offset1, int subid, double* wts, double dist2,
                              vtkSmartPointer<vtkIdList> cellPtIds,
                              vtkSmartPointer<vtkIdList> neighborCellIds, 
                              std::map<int, std::vector<int>>& neighbors, 
                              unsigned int edgenum, bool check = false);
  void smoothCoeffs_alt_help  ( unsigned int channel,
                                unsigned int icell,
                                unsigned int nb,
                                unsigned int nh,
                                unsigned int nl,
                                jPoly<double>* bPm,
                                jPoly<double>* hPm,
                                jPoly<double>* lPm,
                                double* cimg, double* cimg1,
                                double* imgbndb, double* imgbndb1,
                                double* imgbndh, double* imgbndh1,
                                double* imgbndl, double* imgbndl1,
                                double* pcoords10, double* pcoords11,
                                int offset1, int subid, double* wts, double dist2,
                                vtkSmartPointer<vtkIdList> cellPtIds,
                                vtkSmartPointer<vtkIdList> neighborCellIds,
                                std::map<int, std::vector<int>>& neighbors,
                                unsigned int edgenum, 
                                double* cRHS );
  void run  ( );

  ~Compressor();

};

/* Decompression after compression */
struct Decompressor
{
  vtkSmartPointer<vtkPolyData> polytri, polytri1, polytri2, polytri3; 
  vtkSmartPointer<vtkDoubleArray> coeffs, coeffs1, coeffs2, coeffs3;
  vtkSmartPointer<vtkIntArray> offsets, offsets1, offsets2, offsets3;
  vtkSmartPointer<vtkUnsignedIntArray> orders, orders1, orders2, orders3;
  vtkSmartPointer<vtkImageData> imagedata;
  vtkSmartPointer<vtkUnsignedShortArray> colors;
  vtkSmartPointer<vtkPoints> pixels;
  vtkSmartPointer<vtkCellTreeLocator> triloc;
  double a, b, c;
  unsigned int mmax, Mmax;
  double bounds[6];
  bool useMultiChannel;

  Decompressor  ( bool useMultiChannel,
                  const char* channel1, 
                  const char* channel2 = 0, 
                  const char* channel3 = 0 );

  void decompressChannel  ( unsigned int channel  );
  void decompressChannel_help ( unsigned int channel,
                                vtkSmartPointer<vtkPolyData> polytri,
                                vtkSmartPointer<vtkDoubleArray> coeffs,
                                vtkSmartPointer<vtkIntArray> offsets,
                                vtkSmartPointer<vtkUnsignedIntArray> orders );

  void writeImage ( const std::string& pref,
                    const std::string& ext );

  void run ();
  
  ~Decompressor(){}

};

/* wrappers for compression pipeline */
void writeVTP(vtkSmartPointer<vtkPolyData> polytri, const char* ofname);
void writeVTKLegacy(vtkSmartPointer<vtkPolyData> polytri, const char* ofname);
void writeSTL(vtkSmartPointer<vtkPolyData> polytri, const char* ofname);

Triangulator* jcompress_triangulate ( const char* fname, 
                                      unsigned int nSamp,
                                      unsigned int nRuns,
                                      unsigned int mtarget,
                                      bool useMultiChannel,
                                      bool viz );

double jcompress  ( Triangulator* T,
                    const char* fname, 
                    unsigned int nSamp,
                    unsigned int nRuns,
                    unsigned int order,
                    bool useMultiChannel,
                    bool viz );

double jcompress  ( const char* fname, 
                    unsigned int nSamp,
                    unsigned int nRuns,
                    unsigned int order,
                    bool useMultiChannel,
                    bool viz );

void jdecompress  ( bool useMultiChannel,
                    const char* fmt, 
                    const char* channel1,
                    const char* channel2 = 0,
                    const char* channel3 = 0 );


/* benchmarking utilities (reading, computing SSIM) */
vtkSmartPointer<vtkImageData> readImage ( const std::string& pref,
                                          const std::string& ext );

double ssim ( vtkSmartPointer<vtkImageData> img1,
              vtkSmartPointer<vtkImageData> img2 ); 
double ssim ( const char* f1WithExt, const char* f2WithExt );




#endif
