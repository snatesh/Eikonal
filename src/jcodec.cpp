#include <jcodec.hh>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkPolyLine.h>
#include <vtkFeatureEdges.h>
#include <vtkPolyDataEdgeConnectivityFilter.h>
#include <vtkLine.h>
#include <vtkExtractEdges.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkInformation.h>
#include <vtkPNGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkPNMWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkImageWriter.h>
#include <vtkImageCast.h>
#include <vtkImageViewer2.h>
#include <vtkJPEGReader.h>
#include <vtkPNGReader.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAttributeSmoothingFilter.h>
#include <vtkRenderer.h>
#include <vtkImageSSIM.h>
#include <vtkImageReader2.h>
#include <vtkTIFFReader.h>
#include <vtkImageSSIM.h>
#include <vtkPointData.h>
#include <vtkPNMReader.h>
#include <vtkLZMADataCompressor.h>
#include <vtkLZ4DataCompressor.h>
#include <vtkZLibDataCompressor.h>
#include <vtkSTLWriter.h>
#include <vtkDataSetWriter.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageFFT.h>
#include <vtkImageMagnitude.h>
#include <vtkImageToStructuredGrid.h>
#include <vtkImageFourierCenter.h>
#include <vtkImageLogarithmicScale.h>
#include <vtkStructuredGridWriter.h>
#include <map>
#include <vector>
#include <timer.hh>
#include <filesystem>
#include <fstream>
#include <half.hpp>
#include <array>
#include <unordered_set>
#include <algorithm>

using half_float::half;

vtkSmartPointer<vtkDelaunay2D> triangulateUniform ( int* dims, 
                                                    double* origin,
                                                    unsigned int nSamp, 
                                                    vtkIdType& ll, vtkIdType& lr,
                                                    vtkIdType& ur, vtkIdType& ul )
{
  if (nSamp < 2)
  {
    std::cerr << "Initial grid subsampling must have at least 2 points per axis\n";
    exit(1);
  }
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  unsigned int numPoints = dims[0]*dims[1];
  points->SetNumberOfPoints(nSamp * nSamp);
  double xstride =  static_cast<double>((dims[0]-1-origin[0])) / 
                    static_cast<double>((nSamp - 1));
  double ystride =  static_cast<double>((dims[1]-1-origin[0])) / 
                    static_cast<double>((nSamp - 1));
   
  int ptId = 0; double pt[3]; pt[2] = 0;
  for (unsigned int iY = 0; iY < nSamp; ++iY)
  {
    for (unsigned int iX = 0; iX < nSamp; ++iX)
    {
      pt[0] = iX * xstride; pt[1] = iY * ystride;
      points->SetPoint(ptId, pt);
      if (iX == 0 && iY == 0) { ll = ptId; }
      if (iX == nSamp-1 && iY == 0) { lr = ptId; }
      if (iX == nSamp-1 && iY == nSamp-1) { ur = ptId; }
      if (iX == 0 && iY == nSamp-1) { ul = ptId; }
      ptId += 1;  
    }
  }
  
  std::cout << dims[0] << " " << dims[1] << std::endl;
 
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New(); 
  polydata->SetPoints(points);
  vtkSmartPointer<vtkDelaunay2D> triangulator = vtkSmartPointer<vtkDelaunay2D>::New();
  triangulator->SetInputData(polydata);
  triangulator->Update();
  return triangulator;
}

void readQuad ( const std::string& trix, 
                const std::string& triy, 
                const std::string& triw,
                unsigned int N,
                double* X, double* Y, 
                double*  W  )
{

  std::ifstream xfile(trix);
  std::ifstream yfile(triy);
  std::ifstream wfile(triw);
  for (unsigned int i = 0; i < N; ++i)
  {
    xfile >> X[i];
    yfile >> Y[i];
    wfile >> W[i];
  }

}

Triangulator::Triangulator  ( unsigned int _N, unsigned int _m,
                              double _a, double _b, double _c,
                              unsigned int _nSamp, double _bpp_target,
                              const std::string& _trix, 
                              const std::string& _triy,
                              const std::string& _triw,
                              vtkSmartPointer<vtkImageData> _imagedata,
                              vtkSmartPointer<vtkImageInterpolator> _interpolator )
  : N(_N), m(_m), a(_a), b(_b), c(_c), nSamp(_nSamp), 
    bpp_target(_bpp_target), trix(_trix), triy(_triy), triw(_triw)
  {
    unsigned int nthreads = 1;
    this->M = static_cast<unsigned int>(0.5 * (m + 1) * (m + 2));
    this->Pm = new jPoly<double>(N, m, a, b, c, nthreads); 
    this->Pmx = new jPoly<double>(N, m, a+1, b, c+1, nthreads); 
    this->Pmy = new jPoly<double>(N, m, a, b+1, c+1, nthreads); 
    this->Habc  = (double*) calloc((m+2)*(m+2), sizeof(double));
    this->Ha1bc1  = (double*) calloc((m+2)*(m+2), sizeof(double));
    this->Hab1c1 = (double*) calloc((m+2)*(m+2), sizeof(double));
    this->Dx = (double*) calloc(M*M, sizeof(double));
    this->Dy = (double*) calloc(M*M, sizeof(double));
    this->X = (double*) calloc(N, sizeof(double));
    this->Y = (double*) calloc(N, sizeof(double));
    this->W = (double*) calloc(N, sizeof(double));
    sFactors(m+2, a, b, c, Habc);  
    sFactors(m+2, a+1, b, c+1, Ha1bc1);  
    sFactors(m+2, a, b+1, c+1, Hab1c1);  
    dMat(a, b, c, Habc, Ha1bc1, m, 0, Dx);
    dMat(a, b, c, Habc, Hab1c1, m, 1, Dy); 
    readQuad(trix, triy, triw, N, X, Y, W);
    // compute Jacobi interpolation matrices
    this->Pm->computeV(X,Y);
    this->Pmx->computeV(X,Y);
    this->Pmy->computeV(X,Y);
    // storage for computing sizefield in adapative triangulation
    this->cimg = (double*) calloc(M*3, sizeof(double));
    
    this->cdimg = (double*) calloc(M*2, sizeof(double));
    this->dimgr = (double*) calloc(N*2, sizeof(double));
    this->dimgt = (double*) calloc(2*N, sizeof(double));
    // interpolated val storage for each channel
    this->interpc = (double*) calloc(N * 3, sizeof(double));
    this->interpc_jacobi = (double*) calloc(N * 3, sizeof(double));
    
    // read the image
    this->imagedata = _imagedata;
    this->pixels = imagedata->GetPoints();
    imagedata->GetDimensions(this->dims);
    this->Npix = dims[0]*dims[1];
    double origin[3]; imagedata->GetOrigin(origin);
    // get the interpolator
    this->interpolator = _interpolator;
  
    // triangulate uniform subsampled grid on image
    vtkSmartPointer<vtkDelaunay2D> triangulator = 
      triangulateUniform(dims, origin, nSamp, ll, lr, ur, ul);
    polytri = triangulator->GetOutput();
    std::cout << "Num cells on subsampled grid: " 
              << polytri->GetNumberOfCells() << std::endl;

    /* create constraint polydata to ensure the image boundary is preserved
       through each triangulation run
    */
    
    // define a counter-clockwise ordering of lines for bounding box
    vtkSmartPointer<vtkPolyLine> bboxLines = vtkSmartPointer<vtkPolyLine>::New();
    bboxLines->GetPointIds()->SetNumberOfIds(8);
    // line 1
    bboxLines->GetPointIds()->SetId(0, ll);
    bboxLines->GetPointIds()->SetId(1, lr);
    // line 2
    bboxLines->GetPointIds()->SetId(2, lr);
    bboxLines->GetPointIds()->SetId(3, ur);
    // line 3
    bboxLines->GetPointIds()->SetId(4, ur);
    bboxLines->GetPointIds()->SetId(5, ul);
    // line 4
    bboxLines->GetPointIds()->SetId(6, ul);
    bboxLines->GetPointIds()->SetId(7, ll);
    // load the bounding box ordered lines into a cell array
    vtkSmartPointer<vtkCellArray> bboxCells = vtkSmartPointer<vtkCellArray>::New();
    bboxCells->InsertNextCell(bboxLines);
    // create a polydata instance of the bounding box
    this->polyBbox = vtkSmartPointer<vtkPolyData>::New();
    polyBbox->SetLines(bboxCells);

  }

Triangulator::Triangulator  ( unsigned int _N, unsigned int _m,
                              double _a, double _b, double _c,
                              unsigned int _nSamp, unsigned int _nRuns,
                              const std::string& _trix, 
                              const std::string& _triy,
                              const std::string& _triw,
                              vtkSmartPointer<vtkImageData> _imagedata,
                              vtkSmartPointer<vtkImageInterpolator> _interpolator,
                              bool _useMultiChannel )
  : N(_N), m(_m), a(_a), b(_b), c(_c), nSamp(_nSamp), nRuns(_nRuns),
    trix(_trix), triy(_triy), triw(_triw), useMultiChannel(_useMultiChannel)
  {
    unsigned int nthreads = 1;
    this->M = static_cast<unsigned int>(0.5 * (m + 1) * (m + 2));
    this->Pm = new jPoly<double>(N, m, a, b, c, nthreads); 
    this->Pmx = new jPoly<double>(N, m, a+1, b, c+1, nthreads); 
    this->Pmy = new jPoly<double>(N, m, a, b+1, c+1, nthreads); 
    this->Habc  = (double*) calloc((m+2)*(m+2), sizeof(double));
    this->Ha1bc1  = (double*) calloc((m+2)*(m+2), sizeof(double));
    this->Hab1c1 = (double*) calloc((m+2)*(m+2), sizeof(double));
    this->Dx = (double*) calloc(M*M, sizeof(double));
    this->Dy = (double*) calloc(M*M, sizeof(double));
    this->X = (double*) calloc(N, sizeof(double));
    this->Y = (double*) calloc(N, sizeof(double));
    this->W = (double*) calloc(N, sizeof(double));
    sFactors(m+2, a, b, c, Habc);  
    sFactors(m+2, a+1, b, c+1, Ha1bc1);  
    sFactors(m+2, a, b+1, c+1, Hab1c1);  
    dMat(a, b, c, Habc, Ha1bc1, m, 0, Dx);
    dMat(a, b, c, Habc, Hab1c1, m, 1, Dy); 
    readQuad(trix, triy, triw, N, X, Y, W);
    // compute Jacobi interpolation matrices
    this->Pm->computeV(X,Y);
    this->Pmx->computeV(X,Y);
    this->Pmy->computeV(X,Y);
    // storage for computing sizefield in adapative triangulation
    this->cimg = (double*) calloc(M, sizeof(double));
    
    this->cdimg = (double*) calloc(M*2, sizeof(double));
    this->dimgr = (double*) calloc(N*2, sizeof(double));
    this->dimgt = (double*) calloc(2*N, sizeof(double));
    // interpolated val storage for each channel
    this->interpc = (double*) calloc(N, sizeof(double));
    this->interpc_jacobi = (double*) calloc(N, sizeof(double));
    
    // read the image
    //this->imagedata = _imagedata; //smoothImage(_imagedata);
    this->imagedata = smoothImage(_imagedata);
    this->pixels = imagedata->GetPoints();
    imagedata->GetDimensions(this->dims);
    this->Npix = dims[0]*dims[1];
    double origin[3]; imagedata->GetOrigin(origin);
    // get the interpolator
    this->interpolator = _interpolator;
  
  
    // triangulate uniform subsampled grid on image
    vtkSmartPointer<vtkDelaunay2D> triangulator = 
      triangulateUniform(dims, origin, nSamp, ll, lr, ur, ul);
    polytri = triangulator->GetOutput();
    std::cout << "Num cells on subsampled grid: " 
              << polytri->GetNumberOfCells() << std::endl;

    
    /* create constraint polydata to ensure the image boundary is preserved
       through each triangulation run
    */
    
    // define a counter-clockwise ordering of lines for bounding box
    vtkSmartPointer<vtkPolyLine> bboxLines = vtkSmartPointer<vtkPolyLine>::New();
    bboxLines->GetPointIds()->SetNumberOfIds(8);
    // line 1
    bboxLines->GetPointIds()->SetId(0, ll);
    bboxLines->GetPointIds()->SetId(1, lr);
    // line 2
    bboxLines->GetPointIds()->SetId(2, lr);
    bboxLines->GetPointIds()->SetId(3, ur);
    // line 3
    bboxLines->GetPointIds()->SetId(4, ur);
    bboxLines->GetPointIds()->SetId(5, ul);
    // line 4
    bboxLines->GetPointIds()->SetId(6, ul);
    bboxLines->GetPointIds()->SetId(7, ll);
    // load the bounding box ordered lines into a cell array
    vtkSmartPointer<vtkCellArray> bboxCells = vtkSmartPointer<vtkCellArray>::New();
    bboxCells->InsertNextCell(bboxLines);
    // create a polydata instance of the bounding box
    this->polyBbox = vtkSmartPointer<vtkPolyData>::New();
    polyBbox->SetLines(bboxCells);
    if (useMultiChannel)
    { 
      this->polytri1 = vtkSmartPointer<vtkPolyData>::New();
      this->polytri2 = vtkSmartPointer<vtkPolyData>::New();
      this->polytri3 = vtkSmartPointer<vtkPolyData>::New();
      this->polytri1->DeepCopy(polytri);
      this->polytri2->DeepCopy(polytri);
      this->polytri3->DeepCopy(polytri);
    }
    //vtkSmartPointer<vtkImageData> fftimg = imageFFT(_imagedata);
    //vtkSmartPointer<vtkImageToStructuredGrid> sgrid 
    //  = vtkSmartPointer<vtkImageToStructuredGrid>::New();
    //sgrid->SetInputData(fftimg);
    //sgrid->Update();
    //vtkSmartPointer<vtkStructuredGridWriter> gridwriter
    //  = vtkSmartPointer<vtkStructuredGridWriter>::New();
    //gridwriter->SetInputData(sgrid->GetOutput());
    //gridwriter->SetFileName("FFT.vtk");
    //gridwriter->SetFileTypeToBinary();
    //gridwriter->Write();
    
  }

Triangulator::~Triangulator()
{
  if (Pm) { delete Pm; Pm = 0; }
  if (Pmx) { delete Pmx; Pmx = 0; }
  if (Pmy) { delete Pmy; Pmy = 0; }
  if (Habc) { free(Habc); Habc = 0; }
  if (Ha1bc1) { free(Ha1bc1); Ha1bc1 = 0; }
  if (Hab1c1) { free(Hab1c1); Hab1c1 = 0; }
  if (Dx) { free(Dx); Dx = 0; }
  if (Dy) { free(Dy); Dy = 0; }
  if (X) { free(X); X = 0; } 
  if (Y) { free(Y); Y = 0; } 
  if (W) { free(W); W = 0; }
  if (cdimg) { free(cdimg); cdimg = 0; }
  if (dimgr) { free(dimgr); dimgr = 0; }
  if (dimgt) { free(dimgt); dimgt = 0; } 
  if (interpc) { free(interpc); interpc = 0; }
  if (interpc_jacobi) { free(interpc_jacobi); interpc_jacobi = 0; }
}

double Triangulator::getBPP( )
{
  double totalBytes =  polytri->GetNumberOfCells() * (Mtarget * 12 + 6) +
                       polytri->GetNumberOfPoints() * 8; 
  return 8 * totalBytes / this->Npix; // LZ-compressor usually deflates by 3-4x
}

void Triangulator::triangulateEntropyGreedyL2J_help  ( double* interr, vtkSmartPointer<vtkPolyData> polytri )
{
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  double v1t[3], v2t[3], v3t[3], tmpwts[3];
  double err, areaT;
  double interpc_cpy[3];
  vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
  // loop over each triangle 
  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    // get tri verts
    polytri->GetCellPoints(i, ptids);
    polytri->GetPoints()->GetPoint(ptids->GetId(0), v1t);
    polytri->GetPoints()->GetPoint(ptids->GetId(1), v2t);
    polytri->GetPoints()->GetPoint(ptids->GetId(2), v3t);
    areaT = vtkTriangle::TriangleArea(v1t, v2t, v3t);
    // interpolate pixel channel vals onto quadrature points
    for (unsigned int j = 0; j < N; ++j) 
    {
      // get quad point in ref
      xqtr[0] = X[j]; xqtr[1] = Y[j];
      // map ref pt xqtr to tri point xqt
      polytri->GetCell(i)->EvaluateLocation(i, xqtr, xqt, tmpwts);
      interpolator->Interpolate(xqt, interpc_cpy);
      interpc[j] = interpc_cpy[0];
      interpc[j + N] = interpc_cpy[1];
      interpc[j + 2*N] = interpc_cpy[2];
    }
    // coefficients of interpolated data in each channel
    Pm->computeCoeffs(interpc, X, Y, W, cimg, 1);
    Pm->computeCoeffs(interpc+N, X, Y, W, cimg+M, 1);
    Pm->computeCoeffs(interpc+2*N, X, Y, W, cimg+2*M, 1);
    // check error in coefficient expansion
    cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans,
                  N, 3, M, 1.0, Pm->V, N, cimg, M, 0.0, interpc_jacobi, N ); 
    err = 0;
    #pragma omp simd reduction(+:err)
    for (unsigned int j = 0; j < N; ++j)
    {
      err += W[j] * ( std::pow(interpc[j] - interpc_jacobi[j], 2.0)
                      + std::pow(interpc[j+N] - interpc_jacobi[j+N], 2.0)
                      + std::pow(interpc[j+2*N] - interpc_jacobi[j+2*N], 2.0) ) / 3.0;
    }
    interr[i] = areaT * err / 2.0;
  }
}
void Triangulator::triangulateEntropyGreedyL1J_help  ( double* interr, vtkSmartPointer<vtkPolyData> polytri )
{
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  double v1t[3], v2t[3], v3t[3], tmpwts[3];
  double err;
  double interpc_cpy[3], areaT;
  vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
  // loop over each triangle 
  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    // get tri verts
    polytri->GetCellPoints(i, ptids);
    polytri->GetPoints()->GetPoint(ptids->GetId(0), v1t);
    polytri->GetPoints()->GetPoint(ptids->GetId(1), v2t);
    polytri->GetPoints()->GetPoint(ptids->GetId(2), v3t);
    areaT = vtkTriangle::TriangleArea(v1t, v2t, v3t);
    // interpolate pixel channel vals onto quadrature points
    for (unsigned int j = 0; j < N; ++j) 
    {
      // get quad point in ref
      xqtr[0] = X[j]; xqtr[1] = Y[j];
      // map ref pt xqtr to tri point xqt
      polytri->GetCell(i)->EvaluateLocation(i, xqtr, xqt, tmpwts);
      interpolator->Interpolate(xqt, interpc_cpy);
      interpc[j] = interpc_cpy[0];
      interpc[j + N] = interpc_cpy[1];
      interpc[j + 2*N] = interpc_cpy[2];
    }
    // coefficients of interpolated data in each channel
    Pm->computeCoeffs(interpc, X, Y, W, cimg, 1);
    Pm->computeCoeffs(interpc+N, X, Y, W, cimg+M, 1);
    Pm->computeCoeffs(interpc+2*N, X, Y, W, cimg+2*M, 1);
    // check error in coefficient expansion
    cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans,
                  N, 3, M, 1.0, Pm->V, N, cimg, M, 0.0, interpc_jacobi, N ); 

    err = 0;
    #pragma omp simd reduction(+:err)
    for (unsigned int j = 0; j < N; ++j)
    {
      err += W[j] * ( std::abs(interpc[j] - interpc_jacobi[j]) +
                      std::abs(interpc[j+N] - interpc_jacobi[j+N]) +
                      std::abs(interpc[j+2*N] - interpc_jacobi[j+2*N]) ) / 3.0;
    }
    interr[i] = areaT * err / 2.0 / N;
  }
}
void Triangulator::triangulateEntropyGreedyL1P_help  ( double* interr, vtkSmartPointer<vtkPolyData> polytri )
{
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  double err;
  vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
  double interpc_cpy[3], imgref[3];
  vtkSmartPointer<vtkCellTreeLocator> triloc = vtkSmartPointer<vtkCellTreeLocator>::New();
  // build locator on image triangulation
  triloc->Initialize();
  triloc->SetDataSet(polytri);
  triloc->BuildLocator(); 

  // build pix-tri map (pixInTri[j] are the pixels indices in triangle j)
  std::map<int, std::vector<int>> pixInTri;
  int triInd;
  for (int ipt = 0; ipt < pixels->GetNumberOfPoints(); ++ipt)
  {
    triInd = triloc->FindCell(pixels->GetPoint(ipt));
    if (triInd == - 1) { std::cout << "not in cell " << ipt << std::endl; }
    pixInTri[triInd].push_back(ipt);
  } 

  // get max num pix in tri
  auto it = pixInTri.begin(); 
  unsigned int maxpix = 0;
  while (it != pixInTri.end())
  {
    if (it->second.size() > maxpix)
    {
      maxpix = it->second.size();
    }
    ++it;
  }

  // allocate buffer for pixels in tri to ref 
  double* rspix = (double*) calloc(maxpix * 2, sizeof(double));  
  double* img = (double*) calloc(maxpix * 3, sizeof(double));
  double* img_ref = (double*) calloc(maxpix * 3, sizeof(double));
  // large buffer jacobi interp op
  jPoly<double>* Pm = new jPoly<double>(maxpix, M, a, b, c, 1); 
  // iterate over triangles
  double v1t[3], v2t[3], v3t[3], tmpwts[3], r[3], dist2;
  int offset, npix, subid;
  unsigned int m, M;
  unsigned short color;
  for (it = pixInTri.begin(); it != pixInTri.end(); ++it)
  {
    // interpolate pixel channel vals onto quadrature points
    for (unsigned int j = 0; j < N; ++j) 
    {
      // get quad point in ref
      xqtr[0] = X[j]; xqtr[1] = Y[j];
      // map ref pt xqtr to tri point xqt
      polytri->GetCell(it->first)->EvaluateLocation(subid, xqtr, xqt, tmpwts);
      interpolator->Interpolate(xqt, interpc_cpy);
      interpc[j] = interpc_cpy[0];
      interpc[j + N] = interpc_cpy[1];
      interpc[j + 2*N] = interpc_cpy[2];
    }
    // coefficients of interpolated data in each channel
    this->Pm->computeCoeffs(interpc, X, Y, W, cimg, 1);
    this->Pm->computeCoeffs(interpc+N, X, Y, W, cimg+M, 1);
    this->Pm->computeCoeffs(interpc+2*N, X, Y, W, cimg+2*M, 1);
  
    // get pixels inside current triangle,
    // get ref img colors
    // and get position in reference
    npix = it->second.size();
    for (unsigned int i = 0; i < npix; ++i)
    { 
      polytri->GetCell(it->first)->
        EvaluatePosition(pixels->GetPoint(it->second[i]), 
                                          nullptr, 
                                          subid, r, 
                                          dist2, tmpwts);
      
      interpolator->Interpolate(pixels->GetPoint(it->second[i]), imgref); 
      img_ref[i] = imgref[0];
      img_ref[i + npix] = imgref[1];
      img_ref[i + 2*npix] = imgref[2];
      rspix[i] = r[0]; rspix[i+npix] = r[1];
    }

    // compute jpoly interp matrix
    Pm->Nx = npix;
    Pm->computeV(rspix, rspix+npix);
    // check error in coefficient expansion
    cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans,
                  npix, 3, M, 1.0, Pm->V, npix, cimg, M, 0.0, img, npix ); 
    err = 0;
    #pragma omp simd reduction(+:err)
    for (unsigned int j = 0; j < npix; ++j)
    {
      err += (  std::abs(img[j] - img_ref[j]) +
                std::abs(img[j+npix] - img_ref[j+npix]) +
                std::abs(img[j+2*npix] - img_ref[j + 2*npix]) ) / 3.0;
    }
    interr[it->first] = err / npix;
  }
  free(rspix);
  free(img);
  free(img_ref);
}

double Triangulator::triangulateEntropyNoGrad_help  ( double* intgn, unsigned int channel,
                                                      vtkSmartPointer<vtkPolyData> polytri,
                                                      double& stdev)
{
  double v1t[3], v2t[3], v3t[3], tmpwts[3];
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  double err, areaT;
  vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
  // loop over each triangle 
  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    // interpolate pixel channel vals onto quadrature points
    for (unsigned int j = 0; j < N; ++j) 
    {
      // get quad point in ref
      xqtr[0] = X[j]; xqtr[1] = Y[j];
      // map ref pt xqtr to tri point xqt
      polytri->GetCell(i)->EvaluateLocation(i, xqtr, xqt, tmpwts);
      interpc[j] = interpolator->Interpolate(xqt[0], xqt[1], xqt[2], channel);
    }
    // get tri verts
    polytri->GetCellPoints(i, ptids);
    polytri->GetPoints()->GetPoint(ptids->GetId(0), v1t);
    polytri->GetPoints()->GetPoint(ptids->GetId(1), v2t);
    polytri->GetPoints()->GetPoint(ptids->GetId(2), v3t);
    areaT = vtkTriangle::TriangleArea(v1t, v2t, v3t);
    // coefficients of interpolated data in channel
    Pm->computeCoeffs(interpc, X, Y, W, cimg, 1);
    // check error in coefficient expansion
    cblas_dgemv(CblasColMajor, CblasNoTrans, N, M, 1.0, Pm->V, N, cimg, 1, 0.0, interpc_jacobi, 1);

    err = 0;
    #pragma omp simd reduction(+:err)
    for (unsigned int j = 0; j < N; ++j)
    {
      err += W[j] * std::pow(interpc[j] - interpc_jacobi[j], 2.0) ; 
    }
    intgn[i] = std::sqrt(areaT * err / 2.0); 
  }
 
  double sum = 0;
  #pragma omp simd reduction(+:sum)
  for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    sum += intgn[i];
  }
  for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    intgn[i] /= sum;
  } 

  double ave = 0;
  #pragma omp simd reduction(+:ave)
  for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    ave += intgn[i];
  } 
  
  ave /= polytri->GetNumberOfCells();

  
  stdev = 0;
  for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    stdev += std::pow((ave - intgn[i]), 2);
  }
  stdev /= polytri->GetNumberOfCells();
  stdev = std::sqrt(stdev);

  return ave;
}

double Triangulator::triangulateEntropy_help  ( double* intgn, unsigned int channel,
                                                vtkSmartPointer<vtkPolyData> polytri,
                                                double& stdev)
{
  double v1t[3], v2t[3], v3t[3], tmpwts[3];
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  double invIxe[4], gradnorm, areaT;
  vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
  // loop over each triangle 
  double ave_err = 0;
  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    // interpolate pixel channel vals onto quadrature points
    for (unsigned int j = 0; j < N; ++j) 
    {
      // get quad point in ref
      xqtr[0] = X[j]; xqtr[1] = Y[j];
      // map ref pt xqtr to tri point xqt
      polytri->GetCell(i)->EvaluateLocation(i, xqtr, xqt, tmpwts);
      interpc[j] = interpolator->Interpolate(xqt[0], xqt[1], xqt[2], channel);
    }
    // get tri verts
    polytri->GetCellPoints(i, ptids);
    polytri->GetPoints()->GetPoint(ptids->GetId(0), v1t);
    polytri->GetPoints()->GetPoint(ptids->GetId(1), v2t);
    polytri->GetPoints()->GetPoint(ptids->GetId(2), v3t);
    areaT = vtkTriangle::TriangleArea(v1t, v2t, v3t);
    // inverse of deformation map 
    invIxe[0] = v3t[1] - v1t[1];
    invIxe[1] = -(v2t[1] - v1t[1]);
    invIxe[2] = -(v3t[0] - v1t[0]); 
    invIxe[3] = v2t[0] - v1t[0];
    // coefficients of interpolated data in channel
    Pm->computeCoeffs(interpc, X, Y, W, cimg, 1);
    // coefficients of xy derivatives of interpolated data in channel
    cblas_dgemv(CblasColMajor, CblasNoTrans, M, M, 1.0, Dx, M, cimg, 1, 0.0, cdimg, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, M, M, 1.0, Dy, M, cimg, 1, 0.0, cdimg+M, 1);
    // values of derivative at points in rs-triangle
    cblas_dgemv(CblasColMajor, CblasNoTrans, N, M, 1.0, Pmx->V, N, cdimg, 1, 0.0, dimgr, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, N, M, 1.0, Pmy->V, N, cdimg+M, 1, 0.0, dimgr+N, 1);
    // map derivative to xy-triangle
    cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasTrans,
                  2, N, 2, 1.0, invIxe, 2, dimgr, N, 0.0, dimgt, 2 );

    gradnorm = 0;
    #pragma omp simd reduction(+:gradnorm)  
    for (unsigned int j = 0; j < N; ++j)
    {
      gradnorm += std::sqrt(dimgt[j] * dimgt[j] + dimgt[j+N] * dimgt[j+N]) * W[j];
    }
    intgn[i] = areaT * gradnorm / 2.0;
  }
 
  double sum = 0;
  #pragma omp simd reduction(+:sum)
  for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    sum += intgn[i];
  }
  for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    intgn[i] /= sum;
  } 

  double ave = 0;
  #pragma omp simd reduction(+:ave)
  for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    ave += intgn[i];
  } 
  
  ave /= polytri->GetNumberOfCells();
  
  stdev = 0;
  for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    stdev += std::pow((ave - intgn[i]), 2);
  }
  stdev /= polytri->GetNumberOfCells();
  stdev = std::sqrt(stdev);

  return ave;
}

vtkSmartPointer<vtkDelaunay2D> Triangulator::triangulateEntropyGreedy ( double* interr )
{
  if (useMultiChannel)
  {
    std::cerr << "Not supported for different triangulations per channel\n";
    exit(1);
  } 
 

  vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
  cellPtIds->SetNumberOfIds(3);
  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
  idList->SetNumberOfIds(2);
  vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();


  //triangulateEntropyGreedy_help  ( interr, polytri );
  //triangulateEntropyGreedyL2J_help  ( interr, polytri );
  //triangulateEntropyGreedyL1P_help  ( interr, polytri );
  triangulateEntropyGreedyL1J_help  ( interr, polytri );

  vtkSmartPointer<vtkPointSet> pointSet = vtkSmartPointer<vtkPointSet>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  // find cell with lowest error
  double min = HUGE_VAL;
  unsigned int mincell = 0;
  for (unsigned int icell = 0; icell < polytri->GetNumberOfCells(); ++icell)
  {
    int nneighbors = 0;
    polytri->GetCellPoints(icell, cellPtIds);
    // find neighbors of first edge
    idList->SetId(0, cellPtIds->GetId(0));
    idList->SetId(1, cellPtIds->GetId(1));
    polytri->GetCellNeighbors(icell, idList, neighborCellIds);
    nneighbors += neighborCellIds->GetNumberOfIds();
    // find neighbors of second  edge
    idList->SetId(0, cellPtIds->GetId(1));
    idList->SetId(1, cellPtIds->GetId(2));
    polytri->GetCellNeighbors(icell, idList, neighborCellIds);
    nneighbors += neighborCellIds->GetNumberOfIds();
    // find neighbors of this edge 
    idList->SetId(0, cellPtIds->GetId(2));
    idList->SetId(1, cellPtIds->GetId(0));
    polytri->GetCellNeighbors(icell, idList, neighborCellIds);
    nneighbors += neighborCellIds->GetNumberOfIds();
    if (nneighbors == 3)
    {
      if (interr[icell] <= min)
      {
        min = interr[icell]; mincell = icell; 
      }
    }
  }

  polytri->GetCellPoints(mincell, cellPtIds);
  int v0 = cellPtIds->GetId(0);
  int v1 = cellPtIds->GetId(1);
  int v2 = cellPtIds->GetId(2);
  double pt[3];
  vtkIdType ptid; int subid;
  vtkIdType ll1, lr1, ur1, ul1;
  for (unsigned int i = 0; i < polytri->GetNumberOfPoints(); ++i)
  {
    // insert any points which aren't part of mincell
    if (i != v0 && i != v1 && i != v2)
    {
      if (i == ll)
      {
        ll1 = points->InsertNextPoint(polytri->GetPoints()->GetPoint(ll));
      }
      else if (i == lr)
      {
        lr1 = points->InsertNextPoint(polytri->GetPoints()->GetPoint(lr));
      }
      else if (i == ur)
      {
        ur1 = points->InsertNextPoint(polytri->GetPoints()->GetPoint(ur));
      }
      else if (i == ul)
      {
        ul1 = points->InsertNextPoint(polytri->GetPoints()->GetPoint(ul));
      }
      else
      {
        points->InsertNextPoint(polytri->GetPoints()->GetPoint(i));
      }
    }
    else
    {
      if (i == ll)
      {
        ll1 = points->InsertNextPoint(polytri->GetPoints()->GetPoint(ll));
      }
      else if (i == lr)
      {
        lr1 = points->InsertNextPoint(polytri->GetPoints()->GetPoint(lr));
      }
      else if (i == ur)
      {
        ur1 = points->InsertNextPoint(polytri->GetPoints()->GetPoint(ur));
      }
      else if (i == ul)
      {
        ul1 = points->InsertNextPoint(polytri->GetPoints()->GetPoint(ul));
      }
    }
  }

  // redefine topology of image bounding box w.r.t new points 

  pointSet->SetPoints(points);
  ll = ll1; lr = lr1; ur = ur1; ul = ul1;

  vtkSmartPointer<vtkPolyLine> bboxLines = vtkSmartPointer<vtkPolyLine>::New();
  bboxLines->GetPointIds()->SetNumberOfIds(8);
  // line 1
  bboxLines->GetPointIds()->SetId(0, ll);
  bboxLines->GetPointIds()->SetId(1, lr);
  // line 2
  bboxLines->GetPointIds()->SetId(2, lr);
  bboxLines->GetPointIds()->SetId(3, ur);
  // line 3
  bboxLines->GetPointIds()->SetId(4, ur);
  bboxLines->GetPointIds()->SetId(5, ul);
  // line 4
  bboxLines->GetPointIds()->SetId(6, ul);
  bboxLines->GetPointIds()->SetId(7, ll);
  // load the bounding box ordered lines into a cell array
  vtkSmartPointer<vtkCellArray> bboxCells = vtkSmartPointer<vtkCellArray>::New();
  bboxCells->InsertNextCell(bboxLines);
  // create a polydata instance of the bounding box
  polyBbox->Initialize();
  polyBbox->SetLines(bboxCells);
  
  vtkSmartPointer<vtkDelaunay2D> triangulator = vtkSmartPointer<vtkDelaunay2D>::New();
  triangulator->SetInputData(pointSet);
  triangulator->SetSourceData(polyBbox);
  triangulator->Update();
  return triangulator;  
}


vtkSmartPointer<vtkDelaunay2D> Triangulator::triangulateEntropy ( double* intgn, unsigned int channel  )
{
  vtkSmartPointer<vtkPolyData> polytri;
  double ave, ave_gn, ave_err;
  double stdev_gn, stdev1_gn, stdev2_gn;
  double stdev_err, stdev1_err, stdev2_err;
  double* interr;
  if (useMultiChannel)
  {
    switch (channel)
    { 
      case 0:
        polytri = polytri1;
        break;
      case 1:
        polytri = polytri2;
        break;
      case 2:
        polytri = polytri3;
        break;
      default:
      {
        std::cerr << "channel must be 0,1,2\n"; 
        exit(1);
      }
    }
    ave = triangulateEntropy_help ( intgn, channel, polytri, stdev_gn );
  }
  else
  {
    polytri = this->polytri;
    double* intgn1 = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
    double* intgn2 = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));

    interr = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
    double* interr1 = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
    double* interr2 = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));

    ave_gn = ( triangulateEntropy_help ( intgn,  0, polytri, stdev_gn ) +
               triangulateEntropy_help ( intgn1, 1, polytri, stdev1_gn ) +
               triangulateEntropy_help ( intgn2, 2, polytri, stdev2_gn )  ) / 3.0;
    ave_err = ( triangulateEntropyNoGrad_help ( interr,  0, polytri, stdev_err ) +
            triangulateEntropyNoGrad_help ( interr1, 1, polytri, stdev1_err ) +
            triangulateEntropyNoGrad_help ( interr2, 2, polytri, stdev2_err )  ) / 3.0;
    for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
    {
      intgn[i] += intgn1[i] + intgn2[i];
      interr[i] += interr1[i] + interr2[i]; 
      intgn[i] /= 3.0;
      interr[i] /= 3.0;
    }
    free(intgn1);
    free(intgn2);
    free(interr1);
    free(interr2);
    stdev_gn = (stdev_gn + stdev1_gn + stdev2_gn) / 3.0;
    stdev_err = (stdev_err + stdev1_err + stdev2_err) / 3.0;
  }

  vtkSmartPointer<vtkIncrementalOctreePointLocator> locator = 
    vtkSmartPointer<vtkIncrementalOctreePointLocator>::New();
  
  locator->SetDataSet(polytri);
  locator->BuildLocator();
  double v1r[3] = {0.5, 0, 0};
  double v2r[3] = {0.5, 0.5, 0};
  double v3r[3] = {0, 0.5, 0};
  double v1t[3], v2t[3], v3t[3], tmpwts[3];
  double cent[3] = {1./3., 1./3., 0};
  vtkIdType ptid; int subid;
  for (int icell = 0; icell < polytri->GetNumberOfCells(); ++icell)
  {
    if (interr[icell] >= ave_err - 0.5*stdev_err)
    {
      polytri->GetCell(icell)->EvaluateLocation(subid, v1r, v1t, tmpwts);
      polytri->GetCell(icell)->EvaluateLocation(subid, v2r, v2t, tmpwts);
      polytri->GetCell(icell)->EvaluateLocation(subid, v3r, v3t, tmpwts);
      locator->InsertUniquePoint(v1t, ptid);  
      locator->InsertUniquePoint(v2t, ptid);  
      locator->InsertUniquePoint(v3t, ptid);
    }
    //else if (intgn[icell] >= ave_gn)
    //{
    //  polytri->GetCell(icell)->EvaluateLocation(subid, v1r, v1t, tmpwts);
    //  polytri->GetCell(icell)->EvaluateLocation(subid, v2r, v2t, tmpwts);
    //  polytri->GetCell(icell)->EvaluateLocation(subid, v3r, v3t, tmpwts);
    //  locator->InsertUniquePoint(v1t, ptid);  
    //  locator->InsertUniquePoint(v2t, ptid);  
    //  locator->InsertUniquePoint(v3t, ptid);
    //}
  } 

  vtkSmartPointer<vtkPointSet> pointSet0 = vtkSmartPointer<vtkPointSet>::New();
  pointSet0->SetPoints(locator->GetLocatorPoints());

  vtkSmartPointer<vtkDelaunay2D> triangulator0 = vtkSmartPointer<vtkDelaunay2D>::New();
  triangulator0->SetInputData(pointSet0);
  triangulator0->Update();

  vtkSmartPointer<vtkFeatureEdges> featureEdges 
    = vtkSmartPointer<vtkFeatureEdges>::New();
  featureEdges->SetInputData(triangulator0->GetOutput());
  featureEdges->BoundaryEdgesOn();
  featureEdges->FeatureEdgesOff();
  featureEdges->ManifoldEdgesOff();
  featureEdges->NonManifoldEdgesOff();
  featureEdges->Update();
  writeVTP(featureEdges->GetOutput(), "bndry0.vtp");
  vtkSmartPointer<vtkPolyData> polybbox = featureEdges->GetOutput();
  vtkSmartPointer<vtkPoints> linepts;
  vtkIdType linept0, linept1;
  double pt0[3], pt1[3], vpt[2];
  double xstart, ystart, xend, yend; 
  xstart = 0; ystart = 0;
  xend = dims[0]-1; yend = dims[1]-1;
  vtkSmartPointer<vtkCellArray> bboxCells = vtkSmartPointer<vtkCellArray>::New();

  for (unsigned int iline = 0; iline < polybbox->GetNumberOfCells(); ++iline)
  {
    vtkSmartPointer<vtkLine> bboxLines = vtkSmartPointer<vtkLine>::New();
    bboxLines->GetPointIds()->SetNumberOfIds(2);
    linepts = polybbox->GetCell(iline)->GetPoints();
    linepts->GetPoint(0, pt0);
    linepts->GetPoint(1, pt1);
    // constrain lines to boundary
    vpt[0] = pt1[0] - pt0[0];
    vpt[1] = pt1[1] - pt0[1];
    if (vpt[0] != 0 && vpt[1] != 0) // line not parallel to boundary
    {
      // one point is on top/bot boundary
      if (pt0[1] == yend)
      {
        pt1[1] = yend; 
      }
      else if (pt0[1] == ystart)
      {
        pt1[1] = ystart;
      }
      if (pt1[1] == yend)
      {
        pt0[1] = yend;
      }
      else if (pt1[1] == ystart)
      {
        pt0[1] = ystart;
      }
      // both pts off boundary (but not par to bndry)
      else if (pt0[1] != yend && pt1[1] != yend && pt0[0] != pt1[0])
      {
        // closer to/on bottom
        if (pt0[1]-ystart < yend-pt0[1]) 
        { 
          pt0[1] = ystart; 
          pt1[1] = ystart;
        }
        // closer to/on top
        else 
        { 
          pt0[1] = yend; 
          pt1[1] = yend;
        }
      }
      // one point is on left/right boundary
      if (pt0[0] == xend)
      {
        pt1[0] = xend; 
      }
      else if (pt0[0] == xstart)
      {
        pt1[0] = xstart;
      }
      if (pt1[0] == xend)
      {
        pt0[0] = xend;
      }
      else if (pt1[0] == xstart)
      {
        pt0[0] = xstart;
      }
      // both pts off boundary (but not par to bndry)
      else if (pt0[0] != xend && pt1[0] != xend && pt1[1] != pt0[1])
      {
        // closer to/on bottom
        if (pt0[0]-xstart < xend-pt0[0]) 
        { 
          pt0[0] = xstart; 
          pt1[0] = xstart;
        }
        // closer to/on top
        else 
        { 
          pt0[0] = xend; 
          pt1[0] = xend;
        }
      }
    }

 
    locator->InsertUniquePoint(pt0, linept0); 
    locator->InsertUniquePoint(pt1, linept1); 
 
    bboxLines->GetPointIds()->SetId(0,linept0);
    bboxLines->GetPointIds()->SetId(1,linept1);
    bboxCells->InsertNextCell(bboxLines);
  }


  polyBbox->Initialize();
  polyBbox->SetLines(bboxCells);
  vtkSmartPointer<vtkPointSet> pointSet = vtkSmartPointer<vtkPointSet>::New();
  pointSet->SetPoints(locator->GetLocatorPoints());


  vtkSmartPointer<vtkDelaunay2D> triangulator = vtkSmartPointer<vtkDelaunay2D>::New();
  triangulator->SetInputData(pointSet);
  triangulator->SetSourceData(polyBbox);
  triangulator->Update();
  vtkSmartPointer<vtkFeatureEdges> featureEdges1 
    = vtkSmartPointer<vtkFeatureEdges>::New();
  featureEdges1->SetInputData(triangulator->GetOutput());
  featureEdges1->BoundaryEdgesOn();
  featureEdges1->FeatureEdgesOff();
  featureEdges1->ManifoldEdgesOff();
  featureEdges1->NonManifoldEdgesOff();
  featureEdges1->Update();
  writeVTP(featureEdges1->GetOutput(), "bndry.vtp");
  free(interr);
  return triangulator;
}


void Triangulator::run()
{
  if (useMultiChannel)
  {
    for (unsigned int iRun = 0; iRun < nRuns; ++iRun)
    {
      double* intgn1 = (double*) calloc(polytri1->GetNumberOfCells(), sizeof(double));
      double* intgn2 = (double*) calloc(polytri2->GetNumberOfCells(), sizeof(double));
      double* intgn3 = (double*) calloc(polytri3->GetNumberOfCells(), sizeof(double));
      vtkSmartPointer<vtkDelaunay2D> triangulator1 = triangulateEntropy ( intgn1, 0 );
      vtkSmartPointer<vtkDelaunay2D> triangulator2 = triangulateEntropy ( intgn2, 1 );
      vtkSmartPointer<vtkDelaunay2D> triangulator3 = triangulateEntropy ( intgn3, 2 );

      polytri1 = triangulator1->GetOutput();
      polytri2 = triangulator2->GetOutput();
      polytri3 = triangulator3->GetOutput();
      free(intgn1);
      free(intgn2);
      free(intgn3);
      std::cout << "Refinement step: " << iRun + 1 << std::endl;
      std::cout << "Num cells on each grid: "
                << polytri1->GetNumberOfCells() << ", "
                << polytri2->GetNumberOfCells() << ", "
                << polytri3->GetNumberOfCells() << "\n\n";
    }
  }
  else
  {
    double xend = dims[0]-1;
    double yend = dims[1]-1;
    double ptll[3] = {0.0, 0.0, 0.0};    
    double ptlr[3] = {xend, 0.0, 0.0};
    double ptur[3] = {xend, yend, 0.0};
    double ptul[3] = {0.0, yend, 0.0};

    unsigned int iRun = 0;
    double bpp = getBPP();
    bpp_target = 1.25;
    double tmul = 1.5;
    //while (bpp > bpp_target)
    //{
    //  double* interr = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
    //  vtkSmartPointer<vtkDelaunay2D> triangulator = triangulateEntropyGreedy(interr);
    //  free(interr);
    //  polytri = triangulator->GetOutput();
    //  bpp = getBPP();
    //  std::cout << "Unrefinement step: " << iRun + 1 << std::endl;
    //  std::cout << "Num cells and points on grid: "
    //            << polytri->GetNumberOfCells() << " " 
    //            << polytri->GetNumberOfPoints() << "\n";
    //  std::cout << "BPP = " << bpp << std::endl;
    //  iRun += 1;
    //}
    
    //while (bpp < bpp_target * tmul)
    for (unsigned int iRun = 0; iRun < nRuns; ++iRun)
    {
      double* intgn = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
      vtkSmartPointer<vtkDelaunay2D> triangulator = triangulateEntropy ( intgn, 0 );
      vtkSmartPointer<vtkPolyData> polytri_old = polytri;
      polytri = triangulator->GetOutput();
      //bpp = getBPP();
      //if (bpp > bpp_target * tmul)
      //{
      //  polytri = polytri_old;
      //  break;
      //}
      free(intgn);
      std::cout << "Refinement step: " << iRun + 1 << std::endl;
      std::cout << "Num cells on grid: "
                << polytri->GetNumberOfCells() << "\n";
      //std::cout << "BPP = " << bpp << std::endl; 
      //iRun += 1;
    }
  }
}

Compressor::Compressor  ( Triangulator* T, unsigned int order)
{
  unsigned int nthreads = 1;
  this->useMultiChannel = T->useMultiChannel;
  this->morder = order;
  if (useMultiChannel)
  { 
    this->polytri1 = T->polytri1; 
    this->polytri2 = T->polytri2; 
    this->polytri3 = T->polytri3;
    this->coeffs1 = vtkSmartPointer<vtkDoubleArray>::New();
    this->coeffs2 = vtkSmartPointer<vtkDoubleArray>::New();
    this->coeffs3 = vtkSmartPointer<vtkDoubleArray>::New();
    this->offsets1 = vtkSmartPointer<vtkIntArray>::New();
    this->offsets2 = vtkSmartPointer<vtkIntArray>::New();
    this->offsets3 = vtkSmartPointer<vtkIntArray>::New();
    this->orders1 = vtkSmartPointer<vtkUnsignedIntArray>::New();
    this->orders2 = vtkSmartPointer<vtkUnsignedIntArray>::New();
    this->orders3 = vtkSmartPointer<vtkUnsignedIntArray>::New();
  }
  else
  {
    this->polytri = vtkSmartPointer<vtkPolyData>::New();
    this->polytri->DeepCopy(T->polytri);
    this->coeffs = vtkSmartPointer<vtkDoubleArray>::New();
    this->offsets = vtkSmartPointer<vtkIntArray>::New();
    this->orders = vtkSmartPointer<vtkUnsignedIntArray>::New();
    this->triloc = vtkSmartPointer<vtkCellTreeLocator>::New();
    this->triloc->SetDataSet(polytri);
    this->triloc->BuildLocator(); 
  }
  this->imagedata = T->imagedata;
  this->pixels = T->pixels;
  this->interpolator = T->interpolator;
  this->R = (double*) calloc(this->N, sizeof(double));
  this->S = (double*) calloc(this->N, sizeof(double));
  this->W = (double*) calloc(this->N, sizeof(double));
  this->interpc = (double*) calloc(this->N*3, sizeof(double));
  this->a = T->a; this->b = T->b; this->c = T->c;
  this->Pm = new jPoly<double>(N, morder, a, b, c, nthreads); 
  this->Mmax = Pm->Np;
  // storage buffer for img coefs
  this->cimg = (double*) calloc(3 * Mmax, sizeof(double)); 
  readQuad(trix, triy, triw, N, R, S, W);
}

Compressor::~Compressor()
{
  if (Pm) { delete Pm; Pm = 0; }
  if (legq) { delete legq; legq = 0; }
  if (R) { free(R); R = 0; }
  if (S) { free(S); S = 0; }
  if (W) { free(W); W = 0; }
  if (interpc) { free(interpc); interpc = 0; }
  if (cimg) { free(cimg); cimg = 0; }
}

void Compressor::compressChannel_alt1_help  ( vtkSmartPointer<vtkPolyData> polytri,
                                              vtkSmartPointer<vtkDoubleArray> coeffs,
                                              vtkSmartPointer<vtkIntArray> offsets,
                                              vtkSmartPointer<vtkUnsignedIntArray> orders,
                                              double* ave, double* stdev  )
{
  if (useMultiChannel)
  {
    std::cerr << "Not supported for separate triangulations per channel\n";
    exit(1);
  }

  unsigned int N = this->N, buff = 50;
  double* X = (double*) calloc(N+buff, sizeof(double)); 
  double* Y = (double*) calloc(N+buff, sizeof(double)); 
  unsigned int nx = static_cast<unsigned int>(std::floor(std::sqrt(2*N)));
  double* x = (double*) calloc(nx, sizeof(double));
  double xstride = 1.0 / (nx - 1);
  for (unsigned int j = 0; j < nx; ++j) { x[j] = j * xstride; }
  // bottom
  for (unsigned int j = 0; j < nx; ++j) { X[j] = x[j]; }
  // left
  for (unsigned int j = nx+1; j < 2*nx; ++j) { Y[j-1] = x[j-nx]; }
  // hyp
  for (unsigned int j = 2*nx; j < (3*nx-2); ++j)
  {
    X[j-1] = x[j-2*nx+1];
    Y[j-1] = 1 - X[j-1];
  }
  unsigned int Nbndry = nx + (nx-1) + (nx-2);
  unsigned int Nremain = N-Nbndry;
  unsigned int ipt = 1;
  double xx, yy;
  while (ipt < Nremain)
  {
    for (unsigned int ii = 1; ii < nx-1; ++ii)
    {
      for (unsigned int jj = 1; jj < nx-1; ++jj)
      {
        xx = x[ii]; yy = x[jj];
        if (yy < 1-xx)
        {
          X[ipt+Nbndry-1] = xx;
          Y[ipt+Nbndry-1] = yy;
          ipt += 1;
        }
        if (ipt > Nremain) { break; }
      }
      if (ipt > Nremain) { break; }
    }
  }
 
 
  unsigned int Ntot = Nbndry + ipt - 1;

  std::cout << Ntot << std::endl;
  FILE* file = fopen("test.txt", "w");
  for (unsigned int i = 0; i < Ntot; ++i)
  {
    fprintf(file, "%.17g %17g\n", X[i], Y[i]);
  }
  fclose(file);
  
  // allocate buffer for pixels in tri to ref 
  double* interpall = (double*) calloc(3*Ntot, sizeof(double));
  // large buffer jacobi interp op
  jPoly<double>* Pm_buf = new jPoly<double>(Ntot, morder, a, b, c, 1); 
  // storage for coeffs
  unsigned int Mtot = polytri->GetNumberOfCells() * Mmax;
  coeffs->SetNumberOfTuples(Mtot); 
  double* cimg = (double*) calloc(3*Mmax, sizeof(double)); 
 
  // iterate over triangles
  double tmpwts[3], dist2;
  double color[3], coeff[3], offtup[3];
  double mdub = static_cast<double>(morder);
  double ordtup[3] = {mdub, mdub, mdub};
  int npix, subid;
  unsigned int m, M;
  ave[0] = ave[1] = ave[2] = 0; 
  stdev[0] = stdev[1] = stdev[2] = 0;
  int offset = 0;
  bool goodSol;
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  double res[3] = {-1, -1, -1};
  for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    // interpolate pixel channel vals onto uniform sample points
    for (unsigned int j = 0; j < Ntot; ++j) 
    {
      // get quad point in ref
      xqtr[0] = X[j]; xqtr[1] = Y[j];
      // map ref pt xqtr to tri point xqt
      polytri->GetCell(i)->EvaluateLocation(subid, xqtr, xqt, tmpwts);
      interpolator->Interpolate(xqt, color);
      interpall[j] = color[0];
      interpall[j + Ntot] = color[1];
      interpall[j + 2*Ntot] = color[2]; 
    }

    goodSol = Pm_buf->computeInterpCoeffs(interpall, X, Y, cimg, 3, res);
    if (not goodSol)
    {
      std::cout << "bad sol" << std::endl;
      // interpolate pixel channel vals onto quadrature points
      for (unsigned int j = 0; j < N; ++j) 
      {
        // get quad point in ref
        xqtr[0] = R[j]; xqtr[1] = S[j];
        // map ref pt xqtr to tri point xqt
        polytri->GetCell(i)->EvaluateLocation(subid, xqtr, xqt, tmpwts);
        interpolator->Interpolate(xqt, color);
        interpall[j] = color[0];
        interpall[j + N] = color[1];
        interpall[j + 2*N] = color[2]; 
      }
      // if res = -1 (rank deficient), or all res for OD prob are too large
      if (res[0] > 1e-7 && res[1] > 1e-7 && res[2] > 1e-7)
      {
        this->Pm->computeCoeffs(interpall, R, S, W, cimg, 3);
      }
      // last channel sol was fine
      else if (res[0] > 1e-7 && res[1] > 1e-7)
      {
        this->Pm->computeCoeffs(interpall, R, S, W, cimg, 2);
      }
      // if middle channel sol was fine
      else if (res[0] > 1e-7 && res[2] > 1e-7)
      {
        this->Pm->computeCoeffs(interpall, R, S, W, cimg, 1);
        this->Pm->computeCoeffs(interpall+2*N, R, S, W, cimg+2*Mmax, 1);
      }
      // if last two channel sols are fine
      else if (res[0] > 1e-7)
      {
        this->Pm->computeCoeffs(interpall, R, S, W, cimg, 1);
      }
      // if first and last channel sols are fine
      else if (res[1] > 1e-7)
      {
        this->Pm->computeCoeffs(interpall+N, R, S, W, cimg+Mmax, 1);
      }
      // if first and second channel sols are fine
      else if (res[2] > 1e-7)
      {
        this->Pm->computeCoeffs(interpall+2*N, R, S, W, cimg+2*Mmax, 1); 
      }
    }

  
    for (unsigned int j = 0; j < Mmax; ++j)
    {
      coeff[0] = cimg[j]; 
      coeff[1] = cimg[j+Mmax]; 
      coeff[2] = cimg[j+2*Mmax];
      coeffs->SetTuple(offset+j, coeff);
      ave[0] += std::abs(coeff[0]);
      ave[1] += std::abs(coeff[1]);
      ave[2] += std::abs(coeff[2]);
    }   
    offtup[0] = offtup[1] = offtup[2] = offset; 
    offsets->SetTuple(i, offtup);
    orders->SetTuple(i, ordtup);
    offset += Mmax;
  }
  ave[0] /= Mtot;
  ave[1] /= Mtot;
  ave[2] /= Mtot;
  offset = 0;
  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    for (unsigned int j = 0; j < Mmax; ++j)
    {
      coeffs->GetTuple(offset+j, coeff);
      stdev[0] += std::pow(std::abs(coeff[0])-ave[0], 2.0);
      stdev[1] += std::pow(std::abs(coeff[1])-ave[1], 2.0);
      stdev[2] += std::pow(std::abs(coeff[2])-ave[2], 2.0);
    }
    offset += Mmax;
  }

  stdev[0] = std::sqrt(stdev[0] / Mtot); 
  stdev[1] = std::sqrt(stdev[1] / Mtot); 
  stdev[2] = std::sqrt(stdev[2] / Mtot); 

  std::cout << ave[0] << " " << ave[1] << " " << ave[2] << std::endl; 
  std::cout << stdev[0] << " " << stdev[1] << " " << stdev[2] << std::endl;
 
  free(cimg);
  free(interpall);
  free(X);
  free(Y);
  free(x);
  delete Pm_buf;
}

void Compressor::compressChannel_alt_help ( vtkSmartPointer<vtkPolyData> polytri,
                                            vtkSmartPointer<vtkDoubleArray> coeffs,
                                            vtkSmartPointer<vtkIntArray> offsets,
                                            vtkSmartPointer<vtkUnsignedIntArray> orders,
                                            double* ave, double* stdev  )
{
  if (useMultiChannel)
  {
    std::cerr << "Not supported for separate triangulations per channel\n";
    exit(1);
  }


  
  // build pix-tri map (pixInTri[j] are the pixels indices in triangle j)
  std::map<int, std::vector<int>> pixInTri;
  int triInd;
  for (int ipt = 0; ipt < pixels->GetNumberOfPoints(); ++ipt)
  {
    triInd = triloc->FindCell(pixels->GetPoint(ipt));
    if (triInd == - 1) { std::cout << "not in cell " << ipt << std::endl; }
    pixInTri[triInd].push_back(ipt);
  } 

  // get max num pix in tri
  auto it = pixInTri.begin(); 
  unsigned int maxpix = 0;
  while (it != pixInTri.end())
  {
    if (it->second.size() > maxpix)
    {
      maxpix = it->second.size();
    }
    ++it;
  }

  // allocate buffer for pixels in tri to ref 
  double* rspix = (double*) calloc(maxpix * 2, sizeof(double));  
  double* img = (double*) calloc(3*(maxpix > Mmax ? maxpix : Mmax), sizeof(double));
  double* interpall = (double*) calloc(3*this->N, sizeof(double));
  // large buffer jacobi interp op
  jPoly<double>* Pm_buf = new jPoly<double>(maxpix, morder, a, b, c, 1); 
  // storage for coeffs
  unsigned int Mtot = polytri->GetNumberOfCells() * Mmax;
  coeffs->SetNumberOfTuples(Mtot); 
  double* cimg = (double*) calloc(3*Mmax, sizeof(double)); 
 
  // iterate over triangles
  double v1t[3], v2t[3], v3t[3], tmpwts[3], r[3], dist2;
  double color[3], coeff[3], offtup[3];
  double mdub = static_cast<double>(morder);
  double ordtup[3] = {mdub, mdub, mdub};
  int npix, subid;
  unsigned int m, M;
  ave[0] = ave[1] = ave[2] = 0; 
  stdev[0] = stdev[1] = stdev[2] = 0;
  int offset = 0;
  bool goodSol;
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  for (it = pixInTri.begin(); it != pixInTri.end(); ++it)
  {
    // get tri verts
    vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
    polytri->GetCellPoints(it->first, ptids);
    polytri->GetPoints()->GetPoint(ptids->GetId(0), v1t);
    polytri->GetPoints()->GetPoint(ptids->GetId(1), v2t);
    polytri->GetPoints()->GetPoint(ptids->GetId(2), v3t);
    
    // get pixels inside current triangle,
    // get colors at the pixel,
    // and get position in reference
    npix = it->second.size();
    for (unsigned int i = 0; i < npix; ++i)
    {
      polytri->GetCell(it->first)->
        EvaluatePosition(pixels->GetPoint(it->second[i]), 
                         nullptr, subid, r, 
                         dist2, tmpwts);
      rspix[i] = r[0]; rspix[i+npix] = r[1];
      imagedata->GetPointData()->GetArray(0)->GetTuple(it->second[i], color);
      img[i] = color[0];
      img[i + npix] = color[1];
      img[i + 2*npix] = color[2];
    }
    Pm_buf->Nx = npix;
    goodSol = Pm_buf->computeInterpCoeffs(img, rspix, rspix+npix, cimg, 3);
    /*  if we solved an underdetermined problem, a rank-deficient problem,
        or if the residual norm of the overdetermined full-rank problem is too large, 
           use quadrature to evaluate coeffs instead of lsq sol */
    if (not goodSol)
    {
      // interpolate pixel channel vals onto quadrature points
      for (unsigned int j = 0; j < N; ++j) 
      {
        // get quad point in ref
        xqtr[0] = R[j]; xqtr[1] = S[j];
        // map ref pt xqtr to tri point xqt
        polytri->GetCell(it->first)->EvaluateLocation(subid, xqtr, xqt, tmpwts);
        interpolator->Interpolate(xqt, color);
        interpall[j] = color[0];
        interpall[j + N] = color[1];
        interpall[j + 2*N] = color[2]; 
      }
      this->Pm->computeCoeffs(interpall, R, S, W, cimg, 3);
    }
    for (unsigned int j = 0; j < Mmax; ++j)
    {
      coeff[0] = cimg[j]; 
      coeff[1] = cimg[j+Mmax]; 
      coeff[2] = cimg[j+2*Mmax];
      coeffs->SetTuple(offset+j, coeff);
      ave[0] += std::abs(coeff[0]);
      ave[1] += std::abs(coeff[1]);
      ave[2] += std::abs(coeff[2]);
    }   
    offtup[0] = offtup[1] = offtup[2] = offset; 
    offsets->SetTuple(it->first, offtup);
    orders->SetTuple(it->first, ordtup);
    offset += Mmax;
  }
  ave[0] /= Mtot;
  ave[1] /= Mtot;
  ave[2] /= Mtot;
  //offset = 0;
  //for (it = pixInTri.begin(); it != pixInTri.end(); ++it)
  //{
  //  for (unsigned int i = 0; i < Mmax; ++i)
  //  {
  //    coeffs->GetTuple(offset+i, coeff);
  //    stdev[0] += std::pow(std::abs(coeff[0])-ave[0], 2.0);
  //    stdev[1] += std::pow(std::abs(coeff[1])-ave[1], 2.0);
  //    stdev[2] += std::pow(std::abs(coeff[2])-ave[2], 2.0);
  //  }
  //  offset += Mmax;
  //}

  //stdev[0] = std::sqrt(stdev[0] / Mtot); 
  //stdev[1] = std::sqrt(stdev[1] / Mtot); 
  //stdev[2] = std::sqrt(stdev[2] / Mtot); 

  std::cout << ave[0] << " " << ave[1] << " " << ave[2] << std::endl; 
 
  free(rspix);
  free(img);
  free(cimg);
  free(interpall);
  delete Pm_buf;
}


void Compressor::compressChannel_help ( vtkSmartPointer<vtkPolyData> polytri,
                                        vtkSmartPointer<vtkDoubleArray> coeffs,
                                        vtkSmartPointer<vtkIntArray> offsets,
                                        vtkSmartPointer<vtkUnsignedIntArray> orders,
                                        double* ave, double* stdev)
{
  double tmpwts[3];
  double color[3], coeff[3], offtup[3];
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  
  unsigned int Mstart = 7; 
  unsigned int Mtot = polytri->GetNumberOfCells() * (Mmax-Mstart);
  unsigned int MMtot = polytri->GetNumberOfCells() * Mmax;
  coeffs->SetNumberOfTuples(MMtot); 
  int offset = 0;
  ave[0] = ave[1] = ave[2] = 0; 
  stdev[0] = stdev[1] = stdev[2] = 0;
  double mdub = static_cast<double>(morder);
  double ordtup[3] = {mdub, mdub, mdub};
  // loop over each triangle 
  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    // interpolate pixel channel vals onto quadrature points
    for (unsigned int j = 0; j < N; ++j) 
    {
      // get quad point in ref
      xqtr[0] = R[j]; xqtr[1] = S[j];
      // map ref pt xqtr to tri point xqt
      polytri->GetCell(i)->EvaluateLocation(i, xqtr, xqt, tmpwts);
      interpolator->Interpolate(xqt, color);
      interpc[j] = color[0];
      interpc[j + N] = color[1];
      interpc[j + 2*N] = color[2]; 
    }
    Pm->computeCoeffs(interpc, R, S, W, cimg, 3);
  
    for (unsigned int j = 0; j < Mmax; ++j)
    {
      coeff[0] = cimg[j]; 
      coeff[1] = cimg[j+Mmax]; 
      coeff[2] = cimg[j+2*Mmax];
      coeffs->SetTuple(offset+j, coeff);
      if (j >= Mstart)
      {
        ave[0] += std::abs(coeff[0]);
        ave[1] += std::abs(coeff[1]);
        ave[2] += std::abs(coeff[2]);
      }
    }   
    offtup[0] = offtup[1] = offtup[2] = offset; 
    offsets->SetTuple(i, offtup);
    orders->SetTuple(i, ordtup);
    offset += Mmax;
  }
  ave[0] /= Mtot;
  ave[1] /= Mtot;
  ave[2] /= Mtot;
  offset = 0;
  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    for (unsigned int j = Mstart; j < Mmax; ++j)
    {
      coeffs->GetTuple(offset+j, coeff);
      stdev[0] += std::pow(std::abs(coeff[0])-ave[0], 2.0);
      stdev[1] += std::pow(std::abs(coeff[1])-ave[1], 2.0);
      stdev[2] += std::pow(std::abs(coeff[2])-ave[2], 2.0);
    }
    offset += Mmax;
  }

  stdev[0] = std::sqrt(stdev[0] / Mtot); 
  stdev[1] = std::sqrt(stdev[1] / Mtot); 
  stdev[2] = std::sqrt(stdev[2] / Mtot); 

  std::cout << ave[0] << " " << ave[1] << " " << ave[2] << std::endl; 
  std::cout << stdev[0] << " " << stdev[1] << " " << stdev[2] << std::endl;
  
}


void Compressor::compressChannel_help ( unsigned int channel,
                                        vtkSmartPointer<vtkPolyData> polytri,
                                        vtkSmartPointer<vtkDoubleArray> coeffs,
                                        vtkSmartPointer<vtkIntArray> offsets,
                                        vtkSmartPointer<vtkUnsignedIntArray> orders,
                                        double& ave, double& stdev)
{
  double tmpwts[3];
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  double invIxe[4], gradnorm;
  
  unsigned int m = this->morder; 
  unsigned int M = static_cast<unsigned int>(0.5 * (m + 1) * (m + 2));
  int offset = 0; ave = 0; stdev = 0;
  // loop over each triangle 
  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    // interpolate pixel channel vals onto quadrature points
    for (unsigned int j = 0; j < N; ++j) 
    {
      // get quad point in ref
      xqtr[0] = R[j]; xqtr[1] = S[j];
      // map ref pt xqtr to tri point xqt
      polytri->GetCell(i)->EvaluateLocation(i, xqtr, xqt, tmpwts);
      interpc[j] = interpolator->Interpolate(xqt[0], xqt[1], xqt[2], channel);
    }
    Pm->computeCoeffs(interpc, R, S, W, cimg, 1);
    //Pm->computeInterpCoeffs(interpc, R, S, cimg, 1);
    for (unsigned int j = 0; j < M; ++j)
    {
      ave += std::abs(cimg[j]);
    }

    if (useMultiChannel)
    {
      for (unsigned int j = 0; j < M; ++j)
      {

        coeffs->InsertNextValue(cimg[j]);
      }
      offsets->SetValue(i, offset);
      orders->SetValue(i, m);
    }
    else
    {
      //for (unsigned int iblk = 3; iblk  <= m; ++iblk)
      //{
      //  int Mstart = static_cast<int>(0.5 * iblk * (iblk+1));
      //  int Mend = static_cast<int>(0.5 * (iblk+1) * (iblk+2));
      //  double aveblk = 0;
      //  for (int iblkj = Mstart; iblkj < Mend; ++iblkj)
      //  {
      //    aveblk += std::abs(cimg[iblkj]);
      //  }     
      //  aveblk /= (Mend-Mstart);
      //  double stdevblk = 0;
      //  for (int iblkj = Mstart; iblkj < Mend; ++iblkj)
      //  {
      //    stdevblk += std::pow(std::abs(cimg[iblkj]) - aveblk, 2.0);
      //  }     
      //  stdevblk = std::sqrt(stdevblk / (Mend-Mstart));
      //
      //  if (stdevblk / aveblk < 1)       
      //  {
      //    for (int iblkj = Mstart; iblkj < Mend; ++iblkj)
      //    {
      //      if (std::abs(cimg[iblkj]) < aveblk/stdevblk) { cimg[iblkj] = 0; }
      //    }
      //  }
      //}
      for (unsigned int j = 0; j < M; ++j)
      {
        //if (std::abs(cimg[j]) < ave / stdev) { cimg[j] = 0; }
        coeffs->InsertComponent(offset+j, channel, cimg[j]);
      }
      offsets->SetComponent(i, channel, offset);
      orders->SetComponent(i, channel, m);
    }
    offset += M;
  }
  unsigned int Mtot = M*polytri->GetNumberOfCells();
  ave /= Mtot;
}

void Compressor::pruneCoeffs ( )
{
  unsigned int offset = 0;
  double coefftup[3];
  unsigned int m = this->morder; 
  unsigned int M = static_cast<unsigned int>(0.5 * (m + 1) * (m + 2));
  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    for (unsigned int j = 0; j < M; ++j)
    {
      coeffs->GetTuple(offset+j, coefftup);
      for (unsigned int k = 0; k < 3; ++k)
      {
        if (std::abs(coefftup[k]) < aves[k])
        {
          coeffs->SetComponent(offset+j, k, 0);
        }
      }     
    } 
    offset += M; 
  }
}
void Compressor::compressChannel  ( unsigned int channel )
{
  vtkSmartPointer<vtkPolyData> polytri;
  vtkSmartPointer<vtkDoubleArray> coeffs;
  vtkSmartPointer<vtkIntArray> offsets;
  vtkSmartPointer<vtkUnsignedIntArray> orders;
  if (useMultiChannel)
  {
    switch (channel)
    { 
      case 0:
      {
        polytri = polytri1;
        coeffs = coeffs1;
        offsets = offsets1;
        orders = orders1;
        break;
      }
      case 1:
      {
        polytri = polytri2;
        coeffs = coeffs2;
        offsets = offsets2;
        orders = orders2;
        break;
      }
      case 2:
      {
        polytri = polytri3;
        coeffs = coeffs3;
        offsets = offsets3;
        orders = orders3;
        break;
      }
      default:
      {
        std::cerr << "channel must be 0,1,2\n"; 
        exit(1);
      }
    }
    offsets->SetName("offsets");
    offsets->SetNumberOfComponents(1);
    offsets->SetNumberOfTuples(polytri->GetNumberOfCells());
    orders->SetName("orders");
    orders->SetNumberOfComponents(1);
    orders->SetNumberOfTuples(polytri->GetNumberOfCells());
    coeffs->SetName("coeffs");
    coeffs->SetNumberOfComponents(1);
    compressChannel_help ( channel, polytri, coeffs, offsets, orders, ave0, stdev0 );
  }
  else
  {
    polytri = this->polytri;
    coeffs = this->coeffs;
    offsets = this->offsets;
    orders = this->orders;
    offsets->SetName("offsets");
    offsets->SetNumberOfComponents(3);
    offsets->SetNumberOfTuples(polytri->GetNumberOfCells());
    orders->SetName("orders");
    orders->SetNumberOfComponents(3);
    orders->SetNumberOfTuples(polytri->GetNumberOfCells());
    coeffs->SetName("coeffs");
    coeffs->SetNumberOfComponents(3);
    compressChannel_help ( polytri, coeffs, offsets, orders, aves, stdevs );
    //compressChannel_alt1_help ( polytri, coeffs, offsets, orders, aves, stdevs );
    //compressChannel_alt_help ( polytri, coeffs, offsets, orders, aves, stdevs );
    pruneCoeffs();
  }
  
  polytri->GetFieldData()->AddArray(coeffs);
  polytri->GetCellData()->AddArray(offsets); 
  polytri->GetCellData()->AddArray(orders); 


}

double Compressor::smoothCoeffs_help  ( unsigned int channel,
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
                                        unsigned int edgenum, bool check )
{

  double sum1 = 0, sum2 = 0;
  if (neighborCellIds->GetNumberOfIds() > 0)
  {
    neighbors[icell].push_back(neighborCellIds->GetId(0));
    // get coeffs from neighbor cell
    offset1 = offsets->GetComponent(neighborCellIds->GetId(0), channel);
    for (unsigned int i = 0; i < Mmax; ++i)
    {
      cimg1[i] = coeffs->GetComponent(offset1+i, channel);
    }         
    if (edgenum == 0)
    {
      // eval pcoord of endpts of edge in neighbor cell
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(0)),
                          nullptr, subid, pcoords10, dist2, wts );
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(1)),
                          nullptr, subid, pcoords11, dist2, wts );
      // evaluate channel color over bottom edge in icell
      cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                    nleg, Mmax, 1.0, bPm->V, nleg, 
                    cimg, 1, 0.0, imgbnd, 1 );
      // if on bottom edge of neighbor cell
      if ((pcoords10[0] == 0 && pcoords10[1] == 0 &&
           pcoords11[0] == 1 && pcoords11[1] == 0) ||
          (pcoords10[0] == 1 && pcoords10[1] == 0 &&
           pcoords11[0] == 0 && pcoords10[1] == 0))
      {
        // evaluate channel color over bottom edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nleg, Mmax, 1.0, bPm->V, nleg, 
                      cimg1, 1, 0.0, imgbnd1, 1 );
        sum1 = 0, sum2 = 0;
        for (unsigned int k = 0; k < nleg; ++k)
        {
          sum1 += legq->w[k] * imgbnd[k];
          sum2 += legq->w[k] * imgbnd1[k];
        } 
      }
      // if on left edge of neighbor cell
      else if ((pcoords10[0] == 0 && pcoords10[1] == 1 && 
                pcoords11[0] == 0 && pcoords11[1] == 0) ||
               (pcoords10[0] == 0 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1))
      {
        // evaluate channel color over left edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nleg, Mmax, 1.0, lPm->V, nleg, 
                      cimg1, 1, 0.0, imgbnd1, 1 );
        sum1 = 0, sum2 = 0;
        for (unsigned int k = 0; k < nleg; ++k)
        {
          sum1 += legq->w[k] * imgbnd[k];
          sum2 += legq->w[k] * imgbnd1[k];
        } 
      }
      // if on hyp edge of neighbor cell
      else if ((pcoords10[0] == 1 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1) ||
               (pcoords10[0] == 0 && pcoords10[1] == 1 &&
                pcoords11[0] == 1 && pcoords11[1] == 0))
      {
        // evaluate channel color over hyp edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nleg, Mmax, 1.0, hPm->V, nleg, 
                      cimg1, 1, 0.0, imgbnd1, 1 );
        sum1 = 0, sum2 = 0;
        for (unsigned int k = 0; k < nleg; ++k)
        {
          sum1 += legq->w[k] * imgbnd[k];
          sum2 += legq->w[k] * imgbnd1[k];
        } 
      }
    }
    else if (edgenum == 1)
    {
      // eval pcoord of endpts of edge in neighbor cell
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(1)),
                          nullptr, subid, pcoords10, dist2, wts );
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(2)),
                          nullptr, subid, pcoords11, dist2, wts );
      // evaluate channel color over hyp edge in icell
      cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                    nleg, Mmax, 1.0, hPm->V, nleg, 
                    cimg, 1, 0.0, imgbnd, 1 );
      // if on bottom edge of neighbor cell
      if ((pcoords10[0] == 0 && pcoords10[1] == 0 &&
           pcoords11[0] == 1 && pcoords11[1] == 0) ||
          (pcoords10[0] == 1 && pcoords10[1] == 0 &&
           pcoords11[0] == 0 && pcoords10[1] == 0))
      {
        // evaluate channel color over bottom edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nleg, Mmax, 1.0, bPm->V, nleg, 
                      cimg1, 1, 0.0, imgbnd1, 1 );
        sum1 = 0, sum2 = 0;
        for (unsigned int k = 0; k < nleg; ++k)
        {
          sum1 += legq->w[k] * imgbnd[k];
          sum2 += legq->w[k] * imgbnd1[k];
        } 
      }
      // if on left edge of neighbor cell
      else if ((pcoords10[0] == 0 && pcoords10[1] == 1 && 
                pcoords11[0] == 0 && pcoords11[1] == 0) ||
               (pcoords10[0] == 0 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1))
      {
        // evaluate channel color over left edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nleg, Mmax, 1.0, lPm->V, nleg, 
                      cimg1, 1, 0.0, imgbnd1, 1 );
        sum1 = 0, sum2 = 0;
        for (unsigned int k = 0; k < nleg; ++k)
        {
          sum1 += legq->w[k] * imgbnd[k];
          sum2 += legq->w[k] * imgbnd1[k];
        } 
      }
      // if on hyp edge of neighbor cell
      else if ((pcoords10[0] == 1 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1) ||
               (pcoords10[0] == 0 && pcoords10[1] == 1 &&
                pcoords11[0] == 1 && pcoords11[1] == 0))
      {
        // evaluate channel color over hyp edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nleg, Mmax, 1.0, hPm->V, nleg, 
                      cimg1, 1, 0.0, imgbnd1, 1 );
        sum1 = 0, sum2 = 0;
        for (unsigned int k = 0; k < nleg; ++k)
        {
          sum1 += legq->w[k] * imgbnd[k];
          sum2 += legq->w[k] * imgbnd1[k];
        } 
      }

    }
    else if (edgenum == 2)
    {
      // eval pcoord of endpts of edge in neighbor cell
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(2)),
                          nullptr, subid, pcoords10, dist2, wts );
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(0)),
                          nullptr, subid, pcoords11, dist2, wts );
      // evaluate channel color over left edge in icell
      cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                    nleg, Mmax, 1.0, lPm->V, nleg, 
                    cimg, 1, 0.0, imgbnd, 1 );
      // if on bottom edge of neighbor cell
      if ((pcoords10[0] == 0 && pcoords10[1] == 0 &&
           pcoords11[0] == 1 && pcoords11[1] == 0) ||
          (pcoords10[0] == 1 && pcoords10[1] == 0 &&
           pcoords11[0] == 0 && pcoords10[1] == 0))
      {
        // evaluate channel color over bottom edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nleg, Mmax, 1.0, bPm->V, nleg, 
                      cimg1, 1, 0.0, imgbnd1, 1 );
        sum1 = 0, sum2 = 0;
        for (unsigned int k = 0; k < nleg; ++k)
        {
          sum1 += legq->w[k] * imgbnd[k];
          sum2 += legq->w[k] * imgbnd1[k];
        } 
      }
      // if on left edge of neighbor cell
      else if ((pcoords10[0] == 0 && pcoords10[1] == 1 && 
                pcoords11[0] == 0 && pcoords11[1] == 0) ||
               (pcoords10[0] == 0 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1))
      {
        // evaluate channel color over left edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nleg, Mmax, 1.0, lPm->V, nleg, 
                      cimg1, 1, 0.0, imgbnd1, 1 );
        sum1 = 0, sum2 = 0;
        for (unsigned int k = 0; k < nleg; ++k)
        {
          sum1 += legq->w[k] * imgbnd[k];
          sum2 += legq->w[k] * imgbnd1[k];
        } 
      }
      // if on hyp edge of neighbor cell
      else if ((pcoords10[0] == 1 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1) ||
               (pcoords10[0] == 0 && pcoords10[1] == 1 &&
                pcoords11[0] == 1 && pcoords11[1] == 0))
      {
        // evaluate channel color over hyp edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nleg, Mmax, 1.0, hPm->V, nleg, 
                      cimg1, 1, 0.0, imgbnd1, 1 );
        sum1 = 0, sum2 = 0;
        for (unsigned int k = 0; k < nleg; ++k)
        {
          sum1 += legq->w[k] * imgbnd[k];
          sum2 += legq->w[k] * imgbnd1[k];
        } 
      }
    }
  }
  if (check)
  {
    std::cout << sum1 << " " << sum2 << std::endl;
  }
  return sum2;
}

void Compressor::smoothCoeffs_alt_help  ( unsigned int channel,
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
                                          double* cRHS )
{
  /* edges are default numbered as
      - 0 for bottom
      - 1 for hyp
      - 2 for left 
     in parameteric space
  */
  if (neighborCellIds->GetNumberOfIds() > 0)
  {
    neighbors[icell].push_back(neighborCellIds->GetId(0));
    // get coeffs from neighbor cell
    offset1 = offsets->GetComponent(neighborCellIds->GetId(0), channel);
    for (unsigned int i = 0; i < Mmax; ++i)
    {
      cimg1[i] = coeffs->GetComponent(offset1+i, channel);
    }         
    if (edgenum == 0)
    {
      // eval pcoord of endpts of edge in neighbor cell
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(0)),
                          nullptr, subid, pcoords10, dist2, wts );
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(1)),
                          nullptr, subid, pcoords11, dist2, wts );
      // evaluate channel color over bottom edge in icell
      cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                    nb, Mmax, 1.0, bPm->V, nb, 
                    cimg, 1, 0.0, imgbndb, 1 );
      // if on bottom edge of neighbor cell
      if ((pcoords10[0] == 0 && pcoords10[1] == 0 &&
           pcoords11[0] == 1 && pcoords11[1] == 0) ||
          (pcoords10[0] == 1 && pcoords10[1] == 0 &&
           pcoords11[0] == 0 && pcoords10[1] == 0))
      {
        // evaluate channel color over bottom edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nb, Mmax, 1.0, bPm->V, nb, 
                      cimg1, 1, 0.0, imgbndb1, 1 );

        for (unsigned int k = 0; k < nb; ++k)
        {
          cRHS[1+k + channel*Mmax] = imgbndb1[k];
        }                 
      }
      // if on left edge of neighbor cell
      else if ((pcoords10[0] == 0 && pcoords10[1] == 1 && 
                pcoords11[0] == 0 && pcoords11[1] == 0) ||
               (pcoords10[0] == 0 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1))
      {
        // evaluate channel color over left edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nl, Mmax, 1.0, lPm->V, nl, 
                      cimg1, 1, 0.0, imgbndl1, 1 );
        for (unsigned int k = 0; k < nl; ++k)
        {
          cRHS[1+nb+nh+k + channel*Mmax] = imgbndl1[k];
        }
      }
      // if on hyp edge of neighbor cell
      else if ((pcoords10[0] == 1 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1) ||
               (pcoords10[0] == 0 && pcoords10[1] == 1 &&
                pcoords11[0] == 1 && pcoords11[1] == 0))
      {
        // evaluate channel color over hyp edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nh, Mmax, 1.0, hPm->V, nh, 
                      cimg1, 1, 0.0, imgbndh1, 1 );
        for (unsigned int k = 0; k < nh; ++k)
        {
          cRHS[1+nb+k + channel*Mmax] = imgbndh1[k];
        }
      }
    }
    else if (edgenum == 1)
    {
      // eval pcoord of endpts of edge in neighbor cell
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(1)),
                          nullptr, subid, pcoords10, dist2, wts );
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(2)),
                          nullptr, subid, pcoords11, dist2, wts );
      // evaluate channel color over hyp edge in icell
      cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                    nh, Mmax, 1.0, hPm->V, nh, 
                    cimg, 1, 0.0, imgbndh, 1 );
      // if on bottom edge of neighbor cell
      if ((pcoords10[0] == 0 && pcoords10[1] == 0 &&
           pcoords11[0] == 1 && pcoords11[1] == 0) ||
          (pcoords10[0] == 1 && pcoords10[1] == 0 &&
           pcoords11[0] == 0 && pcoords10[1] == 0))
      {
        // evaluate channel color over bottom edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nb, Mmax, 1.0, bPm->V, nb, 
                      cimg1, 1, 0.0, imgbndb1, 1 );
        for (unsigned int k = 0; k < nb; ++k)
        {
          cRHS[1+k + channel*Mmax] = imgbndb1[k];
        }                 
      }
      // if on left edge of neighbor cell
      else if ((pcoords10[0] == 0 && pcoords10[1] == 1 && 
                pcoords11[0] == 0 && pcoords11[1] == 0) ||
               (pcoords10[0] == 0 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1))
      {
        // evaluate channel color over left edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nl, Mmax, 1.0, lPm->V, nl, 
                      cimg1, 1, 0.0, imgbndl1, 1 );
        for (unsigned int k = 0; k < nl; ++k)
        {
          cRHS[1+nb+nh+k + channel*Mmax] = imgbndl1[k];
        }
      }
      // if on hyp edge of neighbor cell
      else if ((pcoords10[0] == 1 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1) ||
               (pcoords10[0] == 0 && pcoords10[1] == 1 &&
                pcoords11[0] == 1 && pcoords11[1] == 0))
      {
        // evaluate channel color over hyp edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nh, Mmax, 1.0, hPm->V, nh, 
                      cimg1, 1, 0.0, imgbndh1, 1 );
        for (unsigned int k = 0; k < nh; ++k)
        {
          cRHS[1+nb+k + channel*Mmax] = imgbndh1[k];
        }
      }
    }
    else if (edgenum == 2)
    {
      // eval pcoord of endpts of edge in neighbor cell
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(2)),
                          nullptr, subid, pcoords10, dist2, wts );
      polytri->GetCell(neighborCellIds->GetId(0))->
        EvaluatePosition( polytri->GetPoints()->GetPoint(cellPtIds->GetId(0)),
                          nullptr, subid, pcoords11, dist2, wts );
      // evaluate channel color over left edge in icell
      cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                    nl, Mmax, 1.0, lPm->V, nl, 
                    cimg, 1, 0.0, imgbndl, 1 );
      // if on bottom edge of neighbor cell
      if ((pcoords10[0] == 0 && pcoords10[1] == 0 &&
           pcoords11[0] == 1 && pcoords11[1] == 0) ||
          (pcoords10[0] == 1 && pcoords10[1] == 0 &&
           pcoords11[0] == 0 && pcoords10[1] == 0))
      {
        // evaluate channel color over bottom edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nb, Mmax, 1.0, bPm->V, nb, 
                      cimg1, 1, 0.0, imgbndb1, 1 );
        for (unsigned int k = 0; k < nb; ++k)
        {
          cRHS[1+k + channel*Mmax] = imgbndb1[k];
        }                 
      }
      // if on left edge of neighbor cell
      else if ((pcoords10[0] == 0 && pcoords10[1] == 1 && 
                pcoords11[0] == 0 && pcoords11[1] == 0) ||
               (pcoords10[0] == 0 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1))
      {
        // evaluate channel color over left edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nl, Mmax, 1.0, lPm->V, nl, 
                      cimg1, 1, 0.0, imgbndl1, 1 );
        for (unsigned int k = 0; k < nl; ++k)
        {
          cRHS[1+nb+nh+k + channel*Mmax] = imgbndl1[k];
        }

      }
      // if on hyp edge of neighbor cell
      else if ((pcoords10[0] == 1 && pcoords10[1] == 0 &&
                pcoords11[0] == 0 && pcoords11[1] == 1) ||
               (pcoords10[0] == 0 && pcoords10[1] == 1 &&
                pcoords11[0] == 1 && pcoords11[1] == 0))
      {
        // evaluate channel color over hyp edge of neighbor cell
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      nh, Mmax, 1.0, hPm->V, nh, 
                      cimg1, 1, 0.0, imgbndh1, 1 );
        for (unsigned int k = 0; k < nh; ++k)
        {
          cRHS[1+nb+k + channel*Mmax] = imgbndh1[k];
        }
      }
    }
  }
}

void Compressor::smoothCoeffs_alt ( )
{
  if (not useMultiChannel)
  {
    vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
    cellPtIds->SetNumberOfIds(3);
    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    idList->SetNumberOfIds(2);
    vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();

    std::map<int, std::vector<int>> neighbors;
    // create edge sampling
    unsigned int n = static_cast<unsigned int>(std::floor( (Mmax - 1.0) / 3.0 ));
    unsigned int nl = n;
    unsigned int nb = n;
    unsigned int nh = (Mmax - 1) - nl - nb;
    std::cout << nl << " " << nb << " " << nh << " " << nh+nl+nb << " " << Mmax << std::endl;
    // construct cart-prod for bottom, left and hyp edges of tref
    double* lX = (double*) calloc(nl, sizeof(double));
    double* lY = (double*) calloc(nl, sizeof(double));
    double* bX = (double*) calloc(nb, sizeof(double));
    double* bY = (double*) calloc(nb, sizeof(double));
    double* hX = (double*) calloc(nh, sizeof(double));
    double* hY = (double*) calloc(nh, sizeof(double));

    double lstride = static_cast<double>(1.0 / (nl - 1));
    double bstride = static_cast<double>(1.0 / (nb - 1));
    double hstride = static_cast<double>(1.0 / (nh - 1));

    for (unsigned int il = 0; il < nl; ++il) { lY[il] = il * lstride; }
    for (unsigned int ib = 0; ib < nb; ++ib) { bX[ib] = ib * bstride; }
    for (unsigned int ih = 0; ih < nh; ++ih)
    {
      hX[ih] = ih * hstride;
      hY[ih] = 1.0 - hX[ih];
    }

    // create polynomial evaluators for the boundary
    jPoly<double>* lPm = new jPoly<double>(lX, lY, nl, Mmax, a, b, c, 1);
    jPoly<double>* bPm = new jPoly<double>(bX, bY, nb, Mmax, a, b, c, 1);
    jPoly<double>* hPm = new jPoly<double>(hX, hY, nh, Mmax, a, b, c, 1);

    // storage for channel coeffs on tri
    double* cimg1 = (double*) calloc(Mmax, sizeof(double));
    double* cimg2 = (double*) calloc(Mmax, sizeof(double));
    double* cimg3 = (double*) calloc(Mmax, sizeof(double));
    double* cimg11 = (double*) calloc(Mmax, sizeof(double));

    // storage for color on bndry
    double* imgbndl = (double*) calloc(nl, sizeof(double));
    double* imgbndl1 = (double*) calloc(nl, sizeof(double));
    double* imgbndb = (double*) calloc(nb, sizeof(double));
    double* imgbndb1 = (double*) calloc(nb, sizeof(double));
    double* imgbndh = (double*) calloc(nh, sizeof(double));
    double* imgbndh1 = (double*) calloc(nh, sizeof(double));
  
    // continuity matrices

    // fully interior
    double* cMat = (double*) calloc(Mmax*Mmax, sizeof(double));
    double* cMat_introw = (double*) calloc(Mmax, sizeof(double));
    double* invcMat = (double*) calloc(Mmax*Mmax, sizeof(double));
    // singular value storage for gesdd
    double* S = (double*) calloc(Mmax, sizeof(double));
    double* invS = (double*) calloc(Mmax*Mmax, sizeof(double));
    double* U = (double*) calloc(Mmax*Mmax, sizeof(double));
    double* VT = (double*) calloc(Mmax*Mmax, sizeof(double));
    double* invSTUT = (double*) calloc(Mmax*Mmax, sizeof(double));

    // continuity right-hand-sides

    // fully interior
    double* cRHS = (double*) calloc(Mmax*3, sizeof(double));
    double* xRHS = (double*) calloc(Mmax*3, sizeof(double));

    /* populate cMat as
      | w^T P(int)  |
      |  P(e01)     | 
      |  P(e02)     | 
      |  P(e03)     | 
    */
   
    // compute and fill row 0 of cMat with \int_T P(x) dx
    cblas_dgemv ( CblasColMajor, CblasTrans,
                  Pm->Nx, Mmax, 1.0, Pm->V, Pm->Nx,
                  W, 1, 0.0, cMat_introw, 1 ); 

    for (unsigned int i = 0; i < Mmax; ++i)
    {
      cMat[i*Mmax] = cMat_introw[i];
    }    
    // fill rows 1 to nb with P(eb)
    for (unsigned int i = 0; i < nb; ++i)
    {
      for (unsigned int j = 0; j < Mmax; ++j)
      {
        cMat[i+1 + Mmax*j] = bPm->V[i + nb*j];
      }
    }
    // fill rows nb+1 to nb+nh with P(eh)
    for (unsigned int i = 0; i < nh; ++i)
    {
      for (unsigned int j = 0; j < Mmax; ++j)
      {
        cMat[i+nb+1 + Mmax*j] = hPm->V[i + nh*j];
      }
    } 
    // fill rows nb+nh+1 to nb+nh+nl with P(el)
    for (unsigned int i = 0; i < nl; ++i)
    {
      for (unsigned int j = 0; j < Mmax; ++j)
      {
        cMat[i+nb+nh+1 + Mmax*j] = lPm->V[i + nl*j];
      }
    }

    // compute SVD of cMat
    if (  LAPACKE_dgesdd  ( LAPACK_COL_MAJOR, 'A',
                            Mmax, Mmax, cMat, Mmax,
                            S, U, Mmax, VT, Mmax ) )
    {
      std::cerr << "ERROR: Lapack dgesdd - Could not compute SVD of cMat\n";
    }
    else
    {
      // compute S^{-1}
      for (unsigned int i = 0; i < Mmax; ++i)
      {
        invS[i + Mmax*i] = 1.0 / S[i];
      } 
      // compute S^{-T}UT
      cblas_dgemm ( CblasColMajor, CblasTrans, CblasTrans,
                    Mmax, Mmax, Mmax, 
                    1, invS, Mmax, U, Mmax, 
                    0, invSTUT, Mmax  );
      // compute cMat^{-1} = VS^{-1}UT
      cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans,
                    Mmax, Mmax, Mmax,
                    1, VT, Mmax, invSTUT, Mmax,
                    0, invcMat, Mmax );      
 
    }

    int subid, offset1, offset2, offset3, offset11;
    double pcoords10[3], pcoords11[3];
    double wts[3], dist2;
    lapack_int rank[1]; 
    // loop over triangles to classify cells
    for (int icell = 0; icell < polytri->GetNumberOfCells(); ++icell)
    {
      // first copy channel coeffs into stor
      offset1 = offsets->GetComponent(icell, 0);
      offset2 = offsets->GetComponent(icell, 1);
      offset3 = offsets->GetComponent(icell, 2);
      // copy  channel coeffs
      for (unsigned int i = 0; i < Mmax; ++i)
      {
        cimg1[i] = coeffs->GetComponent(offset1+i, 0);
        cimg2[i] = coeffs->GetComponent(offset2+i, 1);
        cimg3[i] = coeffs->GetComponent(offset3+i, 2);
      }

      // compute first row of RHS
      cRHS[0]       = cblas_ddot(Mmax, cMat_introw, 1, cimg1, 1);
      cRHS[Mmax]    = cblas_ddot(Mmax, cMat_introw, 1, cimg2, 1);
      cRHS[2*Mmax]  = cblas_ddot(Mmax, cMat_introw, 1, cimg3, 1); 

      // get ptids defining cell
      polytri->GetCellPoints(icell, cellPtIds);

      /* edges are default numbered as
          - 0 for bottom
          - 1 for hyp
          - 2 for left 
         in parameteric space
      */
      
      // find neighbors of first edge
      idList->SetId(0, cellPtIds->GetId(0));
      idList->SetId(1, cellPtIds->GetId(1));
      polytri->GetCellNeighbors(icell, idList, neighborCellIds);

      smoothCoeffs_alt_help ( 0, icell, nb, nh, nl, bPm, hPm, lPm, cimg1, cimg11, 
                              imgbndb, imgbndb1,
                              imgbndh, imgbndh1,
                              imgbndl, imgbndl1,
                              pcoords10, pcoords11,
                              offset11, subid, wts, dist2,
                              cellPtIds, neighborCellIds,
                              neighbors, 0, cRHS  );             
      
      smoothCoeffs_alt_help ( 1, icell, nb, nh, nl, bPm, hPm, lPm,
                              cimg2, cimg11, 
                              imgbndb, imgbndb1,
                              imgbndh, imgbndh1,
                              imgbndl, imgbndl1,
                              pcoords10, pcoords11,
                              offset11, subid, wts, dist2,
                              cellPtIds, neighborCellIds,
                              neighbors, 0, cRHS  );             
      
      smoothCoeffs_alt_help ( 2, icell, nb, nh, nl, bPm, hPm, lPm,
                              cimg3, cimg11, 
                              imgbndb, imgbndb1,
                              imgbndh, imgbndh1,
                              imgbndl, imgbndl1,
                              pcoords10, pcoords11,
                              offset11, subid, wts, dist2,
                              cellPtIds, neighborCellIds,
                              neighbors, 0, cRHS  );             

      // find neighbors of second edge
      idList->SetId(0, cellPtIds->GetId(1));
      idList->SetId(1, cellPtIds->GetId(2));
      polytri->GetCellNeighbors(icell, idList, neighborCellIds);
      
      smoothCoeffs_alt_help ( 0, icell, nb, nh, nl, bPm, hPm, lPm,
                              cimg1, cimg11, 
                              imgbndb, imgbndb1,
                              imgbndh, imgbndh1,
                              imgbndl, imgbndl1,
                              pcoords10, pcoords11,
                              offset11, subid, wts, dist2,
                              cellPtIds, neighborCellIds,
                              neighbors, 1, cRHS  );             
      
      smoothCoeffs_alt_help ( 1, icell, nb, nh, nl, bPm, hPm, lPm,
                              cimg2, cimg11, 
                              imgbndb, imgbndb1,
                              imgbndh, imgbndh1,
                              imgbndl, imgbndl1,
                              pcoords10, pcoords11,
                              offset11, subid, wts, dist2,
                              cellPtIds, neighborCellIds,
                              neighbors, 1, cRHS  );             
      
      smoothCoeffs_alt_help ( 2, icell, nb, nh, nl, bPm, hPm, lPm,
                              cimg3, cimg11, 
                              imgbndb, imgbndb1,
                              imgbndh, imgbndh1,
                              imgbndl, imgbndl1,
                              pcoords10, pcoords11,
                              offset11, subid, wts, dist2,
                              cellPtIds, neighborCellIds,
                              neighbors, 1, cRHS  );             
     
      // find neighbors of third edge
      idList->SetId(0, cellPtIds->GetId(2));
      idList->SetId(1, cellPtIds->GetId(0));
      polytri->GetCellNeighbors(icell, idList, neighborCellIds);
      
      smoothCoeffs_alt_help ( 0, icell, nb, nh, nl, bPm, hPm, lPm,
                              cimg1, cimg11, 
                              imgbndb, imgbndb1,
                              imgbndh, imgbndh1,
                              imgbndl, imgbndl1,
                              pcoords10, pcoords11,
                              offset11, subid, wts, dist2,
                              cellPtIds, neighborCellIds,
                              neighbors, 2, cRHS  );             
      
      smoothCoeffs_alt_help ( 1, icell, nb, nh, nl, bPm, hPm, lPm,
                              cimg2, cimg11, 
                              imgbndb, imgbndb1,
                              imgbndh, imgbndh1,
                              imgbndl, imgbndl1,
                              pcoords10, pcoords11,
                              offset11, subid, wts, dist2,
                              cellPtIds, neighborCellIds,
                              neighbors, 2, cRHS  );             
      
      smoothCoeffs_alt_help ( 2, icell, nb, nh, nl, bPm, hPm, lPm,
                              cimg3, cimg11, 
                              imgbndb, imgbndb1,
                              imgbndh, imgbndh1,
                              imgbndl, imgbndl1,
                              pcoords10, pcoords11,
                              offset11, subid, wts, dist2,
                              cellPtIds, neighborCellIds,
                              neighbors, 2, cRHS  );             

      if (neighbors[icell].size() == 9)
      {
        
        // solve for new coeffs
        cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans,
                      Mmax, 3, Mmax, 1, invcMat, Mmax, 
                      cRHS, Mmax, 0.0,
                      xRHS, Mmax );
        for (unsigned int i = 0; i < Mmax; ++i)
        {
          coeffs->SetComponent(offset1+i, 0, xRHS[i]);
          coeffs->SetComponent(offset2+i, 1, xRHS[i + Mmax]);
          coeffs->SetComponent(offset3+i, 2, xRHS[i + 2*Mmax]);
        }
        if (!(icell % 1000)) 
        { 
          std::cout << "Smoothing Progress : " << std::setprecision(2) 
                    << ((double)icell / (double) polytri->GetNumberOfCells()) * 100 
                    << "%\n";
        }
      } 
    }

    free(lX); free(lY);
    free(bX); free(bY);
    free(hX); free(hY);
    delete lPm;
    delete bPm;
    delete hPm;
    free(cimg1);
    free(cimg2);
    free(cimg3);
    free(cimg11);
    free(imgbndl);
    free(imgbndl1);
    free(imgbndb);
    free(imgbndb1);
    free(imgbndh);
    free(imgbndh1);
    free(cMat);
    free(cMat_introw);
    free(cRHS);
    free(S);
    free(U);  
    free(VT);
    free(invcMat);
    free(invS);
    free(invSTUT);
    free(xRHS);
  }
  else
  {
    std::cerr << "not supported for separate triangulations per channel\n";
    exit(1);
  }
  

}

void Compressor::smoothCoeffs ( )
{
  if (not useMultiChannel)
  {
    celltypes = vtkSmartPointer<vtkIntArray>::New();
    celltypes->SetNumberOfComponents(1);
    celltypes->SetNumberOfTuples(polytri->GetNumberOfCells());
    celltypes->SetName("celltype");

    vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
    cellPtIds->SetNumberOfIds(3);
    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    idList->SetNumberOfIds(2);
    vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();

    std::map<int, std::vector<int>> neighbors;
    // create legendre quad
    unsigned int nleg = 30;
    this->legq = new legQuad<double>(nleg); 
    legq->shift();
    // construct cart-prod for bottom, left and hyp edges of tref
    double* lX = (double*) calloc(nleg, sizeof(double));
    double* lY = (double*) calloc(nleg, sizeof(double));
    double* bX = (double*) calloc(nleg, sizeof(double));
    double* bY = (double*) calloc(nleg, sizeof(double));
    double* hX = (double*) calloc(nleg, sizeof(double));
    double* hY = (double*) calloc(nleg, sizeof(double));

    for (unsigned int i = 0; i < nleg; ++i)
    { 
      lY[i] = legq->x[i];
      bX[i] = legq->x[i];
      hX[i] = legq->x[i];
      hY[i] = 1.0 - hX[i];
    }
    // create polynomial evaluators for the boundary
    jPoly<double>* lPm = new jPoly<double>(lX, lY, nleg, Mmax, a, b, c, 1);
    jPoly<double>* bPm = new jPoly<double>(bX, bY, nleg, Mmax, a, b, c, 1);
    jPoly<double>* hPm = new jPoly<double>(hX, hY, nleg, Mmax, a, b, c, 1);
    

    // storage for channel coeffs on tri
    double* cimg1 = (double*) calloc(Mmax, sizeof(double));
    double* cimg2 = (double*) calloc(Mmax, sizeof(double));
    double* cimg3 = (double*) calloc(Mmax, sizeof(double));
    double* cimg11 = (double*) calloc(Mmax, sizeof(double));

    // storage for color on bndry
    double* imgbnd = (double*) calloc(nleg, sizeof(double));
    double* imgbnd1 = (double*) calloc(nleg, sizeof(double));
  
    // continuity matrices

    // fully interior
    double* cMat = (double*) calloc((N+3)*Mmax, sizeof(double));
    double* invcMat = (double*) calloc(Mmax*(N+3), sizeof(double));
    // on boundary with 2 neighbors
    double* cMat_bnd12 = (double*) calloc((N+2)*Mmax, sizeof(double));
    double* cMat_bnd13 = (double*) calloc((N+2)*Mmax, sizeof(double));
    double* cMat_bnd23 = (double*) calloc((N+2)*Mmax, sizeof(double));
    // on boundary with 1 neighbor
    double* cMat_bnd1 = (double*) calloc((N+1)*Mmax, sizeof(double));
    double* cMat_bnd2 = (double*) calloc((N+1)*Mmax, sizeof(double));
    double* cMat_bnd3 = (double*) calloc((N+1)*Mmax, sizeof(double));
    // singular value storage for gesdd
    double* S = (double*) calloc(N+3, sizeof(double));
    double* invS = (double*) calloc((N+3)*Mmax, sizeof(double));
    double* U = (double*) calloc((N+3)*(N+3), sizeof(double));
    double* VT = (double*) calloc(Mmax*Mmax, sizeof(double));
    double* invSTUT = (double*) calloc(Mmax*(N+3), sizeof(double));

    // continuity right-hand-sides

    // fully interior
    double* cRHS = (double*) calloc((N+3)*3, sizeof(double));
    double* xRHS = (double*) calloc(Mmax*3, sizeof(double));
    // on boundary with 1 neighbor
    double* cRHS_bnd1 = (double*) calloc(N+1, sizeof(double));
    // on boundary with 2 neighbors
    double* cRHS_bnd2 = (double*) calloc(N+2, sizeof(double));
    

    /* populate cMat as
      | P(int)    |
      | w^T P(e01)| 
      | w^T P(e02)| 
      | w^T P(e03)| 
    */
    
    double pmv;
    for (unsigned int j = 0; j < Mmax; ++j)
    {
      for (unsigned int i = 0; i < N; ++i)
      {
        pmv = Pm->V[i + N*j]; 
        cMat[i + (N+3)*j]       = pmv;
        cMat_bnd12[i + (N+2)*j] = pmv; 
        cMat_bnd13[i + (N+2)*j] = pmv; 
        cMat_bnd23[i + (N+2)*j] = pmv;
        cMat_bnd1[i + (N+1)*j]  = pmv;
        cMat_bnd2[i + (N+1)*j]  = pmv;
        cMat_bnd3[i + (N+1)*j]  = pmv;
      } 
    }


    // fully interior 
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, bPm->V, nleg,
                  legq->w, 1, 0.0, &cMat[N], N+3 ); 
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, hPm->V, nleg,
                  legq->w, 1, 0.0, &cMat[N+1], N+3 ); 
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, lPm->V, nleg,
                  legq->w, 1, 0.0, &cMat[N+2], N+3 );
    // compute SVD of cMat
    if (  LAPACKE_dgesdd  ( LAPACK_COL_MAJOR, 'A',
                            N+3, Mmax, cMat, N+3,
                            S, U, N+3, VT, Mmax ) )
    {
      std::cerr << "ERROR: Lapack dgesdd - Could not compute SVD of cMat\n";
    }
    else
    {
      // compute S^{-1}
      for (unsigned int i = 0; i < Mmax; ++i)
      {
        invS[i + (N+3)*i] = 1.0 / S[i];
      } 
      // compute S^{-T}UT
      cblas_dgemm ( CblasColMajor, CblasTrans, CblasTrans,
                    Mmax, N+3, N+3, 
                    1, invS, N+3, U, N+3, 
                    0, invSTUT, Mmax  );
      // compute cMat^{-1} = VS^{-1}UT
      cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans,
                    Mmax, N+3, Mmax,
                    1, VT, Mmax, invSTUT, Mmax,
                    0, invcMat, Mmax );      
 
    }
 

    // two neighbors - bottom and hyp
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, bPm->V, nleg,
                  legq->w, 1, 0.0, &cMat_bnd12[N], N+2 ); 
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, hPm->V, nleg,
                  legq->w, 1, 0.0, &cMat_bnd12[N+1], N+2 ); 
    // two neighbors - bottom and left  
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, bPm->V, nleg,
                  legq->w, 1, 0.0, &cMat_bnd13[N], N+2 ); 
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, lPm->V, nleg,
                  legq->w, 1, 0.0, &cMat_bnd13[N+1], N+2 ); 
    // two neighbors - hyp and left 
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, hPm->V, nleg,
                  legq->w, 1, 0.0, &cMat_bnd23[N], N+2 ); 
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, lPm->V, nleg,
                  legq->w, 1, 0.0, &cMat_bnd23[N+1], N+2 ); 

    // one neighbor - bottom
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, bPm->V, nleg,
                  legq->w, 1, 0.0, &cMat_bnd1[N], N+1 ); 
    // one neighbor - hyp
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, hPm->V, nleg,
                  legq->w, 1, 0.0, &cMat_bnd2[N], N+1 ); 
    // one neighbor - left
    cblas_dgemv ( CblasColMajor, CblasTrans, 
                  nleg, Mmax, 1.0, lPm->V, nleg,
                  legq->w, 1, 0.0, &cMat_bnd3[N], N+1 ); 
    

    int subid, offset1, offset2, offset3, offset11;
    double pcoords10[3], pcoords11[3];
    double wts[3], dist2;
    double sum1, sum2, sum3;
    double sum11, sum21, sum31;
    double sum12, sum22, sum32;
    lapack_int rank[1]; 
    // loop over triangles to classify cells
    for (int icell = 0; icell < polytri->GetNumberOfCells(); ++icell)
    {
      // first copy channel coeffs into stor
      offset1 = offsets->GetComponent(icell, 0);
      offset2 = offsets->GetComponent(icell, 1);
      offset3 = offsets->GetComponent(icell, 2);
      // copy  channel coeffs
      for (unsigned int i = 0; i < Mmax; ++i)
      {
        cimg1[i] = coeffs->GetComponent(offset1+i, 0);
        cimg2[i] = coeffs->GetComponent(offset2+i, 1);
        cimg3[i] = coeffs->GetComponent(offset3+i, 2);
      }

      // get ptids defining cell
      polytri->GetCellPoints(icell, cellPtIds);


      /* edges are default numbered as
          - 0 for bottom
          - 1 for hyp
          - 2 for left 
         in parameteric space
      */
      
      // find neighbors of first edge
      idList->SetId(0, cellPtIds->GetId(0));
      idList->SetId(1, cellPtIds->GetId(1));
      polytri->GetCellNeighbors(icell, idList, neighborCellIds);

      cblas_dcopy ( N, cRHS, 1, cRHS_bnd1, 1 );
      cblas_dcopy ( N, cRHS, 1, cRHS_bnd2, 1 );

      sum1 = smoothCoeffs_help  ( 0, icell, nleg, lPm, bPm, hPm,
                                  cimg1, cimg11, imgbnd, imgbnd1,
                                  pcoords10, pcoords11,
                                  offset11, subid, wts, dist2,
                                  cellPtIds, neighborCellIds,
                                  neighbors, 0 );              
      sum11 = smoothCoeffs_help ( 1, icell, nleg, lPm, bPm, hPm,
                                  cimg2, cimg11, imgbnd, imgbnd1,
                                  pcoords10, pcoords11,
                                  offset11, subid, wts, dist2,
                                  cellPtIds, neighborCellIds,
                                  neighbors, 0 );              
      sum12 = smoothCoeffs_help ( 2, icell, nleg, lPm, bPm, hPm,
                                  cimg3, cimg11, imgbnd, imgbnd1,
                                  pcoords10, pcoords11,
                                  offset11, subid, wts, dist2,
                                  cellPtIds, neighborCellIds,
                                  neighbors, 0 );              

      // find neighbors of second edge
      idList->SetId(0, cellPtIds->GetId(1));
      idList->SetId(1, cellPtIds->GetId(2));
      polytri->GetCellNeighbors(icell, idList, neighborCellIds);
      sum2 = smoothCoeffs_help  ( 0, icell, nleg, lPm, bPm, hPm,
                                  cimg1, cimg11, imgbnd, imgbnd1,
                                  pcoords10, pcoords11,
                                  offset11, subid, wts, dist2,
                                  cellPtIds, neighborCellIds,
                                  neighbors, 1 );              
      sum21 = smoothCoeffs_help ( 1, icell, nleg, lPm, bPm, hPm,
                                  cimg2, cimg11, imgbnd, imgbnd1,
                                  pcoords10, pcoords11,
                                  offset11, subid, wts, dist2,
                                  cellPtIds, neighborCellIds,
                                  neighbors, 1 );              
      sum22 = smoothCoeffs_help ( 2, icell, nleg, lPm, bPm, hPm,
                                  cimg3, cimg11, imgbnd, imgbnd1,
                                  pcoords10, pcoords11,
                                  offset11, subid, wts, dist2,
                                  cellPtIds, neighborCellIds,
                                  neighbors, 1 );              

      // find neighbors of third edge
      idList->SetId(0, cellPtIds->GetId(2));
      idList->SetId(1, cellPtIds->GetId(0));
      polytri->GetCellNeighbors(icell, idList, neighborCellIds);
      sum3 = smoothCoeffs_help  ( 0, icell, nleg, lPm, bPm, hPm,
                                  cimg1, cimg11, imgbnd, imgbnd1,
                                  pcoords10, pcoords11,
                                  offset11, subid, wts, dist2,
                                  cellPtIds, neighborCellIds,
                                  neighbors, 2 );             
      sum31 = smoothCoeffs_help ( 1, icell, nleg, lPm, bPm, hPm,
                                  cimg2, cimg11, imgbnd, imgbnd1,
                                  pcoords10, pcoords11,
                                  offset11, subid, wts, dist2,
                                  cellPtIds, neighborCellIds,
                                  neighbors, 2 );             
      sum32 = smoothCoeffs_help ( 2, icell, nleg, lPm, bPm, hPm,
                                  cimg3, cimg11, imgbnd, imgbnd1,
                                  pcoords10, pcoords11,
                                  offset11, subid, wts, dist2,
                                  cellPtIds, neighborCellIds,
                                  neighbors, 2 );             
      if (neighbors[icell].size() == 9)
      {
        // populate cRHS for channel 0 
        cblas_dgemv ( CblasColMajor, CblasNoTrans,
                      N, Mmax, 1.0, Pm->V, N,
                      cimg1, 1, 0.0, cRHS, 1 );      
        cRHS[N] = sum1;
        cRHS[N+1] = sum2;
        cRHS[N+2] = sum3;
        // populate cRHS for channel 1
        cblas_dgemv ( CblasColMajor, CblasNoTrans,
                      N, Mmax, 1.0, Pm->V, N,
                      cimg2, 1, 0.0, &cRHS[N+3], 1 );      
        cRHS[N + (N+3)] = sum11;
        cRHS[N+1 + (N+3)] = sum21;
        cRHS[N+2 + (N+3)] = sum31;
        // populate cRHS for channel 2
        cblas_dgemv ( CblasColMajor, CblasNoTrans,
                      N, Mmax, 1.0, Pm->V, N,
                      cimg3, 1, 0.0, &cRHS[2*(N+3)], 1 );      
        cRHS[N + 2*(N+3)] = sum12;
        cRHS[N+1 + 2*(N+3)] = sum22;
        cRHS[N+2 + 2*(N+3)] = sum32;
        
        // solve for new coeffs
        cblas_dgemm ( CblasColMajor, CblasNoTrans, CblasNoTrans,
                      Mmax, 3, N+3, 1, invcMat, Mmax, 
                      cRHS, N+3, 0.0,
                      xRHS, Mmax );
        for (unsigned int i = 0; i < Mmax; ++i)
        {
          coeffs->SetComponent(offset1+i, 0, xRHS[i]);
          coeffs->SetComponent(offset2+i, 1, xRHS[i + Mmax]);
          coeffs->SetComponent(offset3+i, 2, xRHS[i + 2*Mmax]);
        }
        if (!(icell % 1000)) 
        { 
          std::cout << "Smoothing Progress : " << std::setprecision(2) 
                    << ((double)icell / (double) polytri->GetNumberOfCells()) * 100 
                    << "%\n";
        }
      } 
    }
 

    //for (auto it = neighbors.begin(); it != neighbors.end(); ++it)
    //{
    //  for (unsigned int i = 0; i < it->second.size(); ++i)
    //  {
    //    std::cout << it->second[i] << " ";
    //  }
    //  std::cout << std::endl;
    //  celltypes->SetComponent(it->first, 0, it->second.size());
    //} 
    // 
    //polytri->GetCellData()->AddArray(celltypes); 
    //std::cout << N << " " << Mmax << std::endl;


    free(lX); free(lY);
    free(bX); free(bY);
    free(hX); free(hY);
    delete lPm;
    delete bPm;
    delete hPm;
    free(cimg1);
    free(cimg2);
    free(cimg3);
    free(cimg11);
    free(imgbnd);
    free(imgbnd1);
    free(cMat);
    free(cMat_bnd12);
    free(cMat_bnd13);
    free(cMat_bnd23);
    free(cMat_bnd1);
    free(cMat_bnd2);
    free(cMat_bnd3);
    free(cRHS);
    free(cRHS_bnd1);
    free(cRHS_bnd2);
    free(S);
    free(U);  
    free(VT);
    free(invcMat);
    free(invS);
    free(invSTUT);
    free(xRHS);

    /* TODO: Outline for coefficient smoothing 

    1) neighbors[i] gives the triangles in
       the element patch of triangle i
    2) for each of these adjacent triangles
        - get the shared edge and evaluate
          the end points parametrically
          - if p1[0] == 0 and p1[1] == 0 and 
               p2[0] == 1 and p2[1] == 0,
            we are on the bottom edge of Tref
          - if p1[0] == 1 and p1[1] == 0 and 
               p2[0] == 0 and p2[1] == 1,
            we are on the hypotenuse of Tref
          - if p1[0] == 0 and p1[1] == 1 and
               p2[0] == 0 and p2[1] == 0,
            we are on the left edge of Tref
          - NOTE: must check that this ordering is always the case
    3) based on the edge location, create an array of (x,y)
       tuples where 
          - if on bottom (x,y) = (legx, 0)
          - if on hyp    (x,y) = (legx, 1-legx)
          - if on left   (x,y) = (0, legx)
          - NOTE: these can be precomputed once and re-used
          - do this for each shared edge
    4) Assemble the matrix (eij is edge between triangle i an j)
        | P(interior) |
        | w'*P(ei1)   |
        | w'*P(ei2)   |
        | w'*P(ei3)   |
      - NOTE: for boundary triangles i, there are either 1 or 2
              adjacent triangles, so the number of weighted
              rows is reduced by 1 or 2.
      - NOTE: this can be precomputed and re-used 
    5) Assemble the right hand side vector 
        (ci_k are coeffs on triangle i for channel j)
        (adjacent triangles are numbered 1,2,3)
        | P(interior) * [ci_1 ci_2 ci_3] |
        | w'*P(ei1)   * [c1_1 c1_2 c1_3] |
        | w'*P(ei2)   * [c2_1 c2_2 c2_3] |
        | w'*P(ei3)   * [c3_1 c3_2 c3_3] |
      - NOTE: for boundary triangles i, there are either 1 or 2
              adjacent triangles, so the number of weighted
              rows is reduced by 1 or 2.
      - NOTE: Should use precomputed matrix to evaluate this RHS 
   6) Solve in the least-squares sense for coeffs 
      which enforce these continuity conditions
   7) Update the coeffs array (indexed by offset) in 
      the polytri instance 
  */


  }
  else
  {
    std::cerr << "not supported for separate triangulations per channel\n";
    exit(1);
  }
}

void Compressor::run  ( )
{
  if (useMultiChannel)
  { 
    compressChannel(0);
    std::cout << "compressed channel 1" << std::endl;
    compressChannel(1);
    std::cout << "compressed channel 2" << std::endl;
    compressChannel(2);
    std::cout << "compressed channel 3" << std::endl;
  }
  else
  {
    std::cout << "Number of triangles: " 
              << polytri->GetNumberOfCells() << std::endl;
    std::cout << "Order per triangle: " << morder << std::endl; 
    compressChannel(0);
    //smoothCoeffs_alt();
    //smoothCoeffs();
    unsigned int M = static_cast<unsigned int>(0.5 * (morder+1)*(morder+2));
    /* (ntri*(mcoeff floats/channel)*(3 channels)*(4ytes/float) 
      + ntri*(3 short indices)*(2 bytes/short) 
      + npts * (2 floats for x,y)*(4 bytes/float) 
    */
    totalBytes =  polytri->GetNumberOfCells() * (M * 12 + 6) +
                  polytri->GetNumberOfPoints() * 8; 
    //std::cout << "compressed all channels into "
    //          << totalBytes / 1e6 << " MB\n";
  }
}

Decompressor::Decompressor  ( bool _useMultiChannel,
                              const char* channel1, 
                              const char* channel2, 
                              const char* channel3  )
{
  if (_useMultiChannel && (!channel2 || !channel3))
  {
    std::cerr << "missing channel files" << std::endl;
  }
  this->useMultiChannel = _useMultiChannel;
  this->triloc = vtkSmartPointer<vtkCellTreeLocator>::New();

  if (useMultiChannel)
  {
    // read compressed channel files
    vtkSmartPointer<vtkXMLPolyDataReader> reader1 =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
    vtkSmartPointer<vtkXMLPolyDataReader> reader2 =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
    vtkSmartPointer<vtkXMLPolyDataReader> reader3 =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();

    reader1->SetFileName(channel1); reader1->Update();
    reader2->SetFileName(channel2); reader2->Update();
    reader3->SetFileName(channel3); reader3->Update();

    this->polytri1 = reader1->GetOutput();
    this->polytri2 = reader2->GetOutput();
    this->polytri3 = reader3->GetOutput();
    // read coefficient and offset arrays
    this->offsets1 = vtkIntArray::SafeDownCast(polytri1->GetCellData()->GetAbstractArray(0));
    this->offsets2 = vtkIntArray::SafeDownCast(polytri2->GetCellData()->GetAbstractArray(0));
    this->offsets3 = vtkIntArray::SafeDownCast(polytri3->GetCellData()->GetAbstractArray(0));
    this->orders1 = vtkUnsignedIntArray::SafeDownCast(polytri1->GetCellData()->GetAbstractArray(1));
    this->orders2 = vtkUnsignedIntArray::SafeDownCast(polytri2->GetCellData()->GetAbstractArray(1));
    this->orders3 = vtkUnsignedIntArray::SafeDownCast(polytri3->GetCellData()->GetAbstractArray(1));
    this->coeffs1 = vtkDoubleArray::SafeDownCast(polytri1->GetFieldData()->GetAbstractArray(0));
    this->coeffs2 = vtkDoubleArray::SafeDownCast(polytri2->GetFieldData()->GetAbstractArray(0));
    this->coeffs3 = vtkDoubleArray::SafeDownCast(polytri3->GetFieldData()->GetAbstractArray(0));
    // initialize decompressed image data
    polytri1->GetBounds(this->bounds);
  }
  else
  {
    // read compressed channel file
    vtkSmartPointer<vtkXMLPolyDataReader> reader =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(channel1); reader->Update();
    this->polytri = reader->GetOutput();
    
    // read coefficient and offset arrays
    this->offsets = vtkIntArray::SafeDownCast(polytri->GetCellData()->GetAbstractArray(0));
    this->orders = vtkUnsignedIntArray::SafeDownCast(polytri->GetCellData()->GetAbstractArray(1));
    this->coeffs = vtkDoubleArray::SafeDownCast(polytri->GetFieldData()->GetAbstractArray(0));
    polytri->GetBounds(this->bounds);
    
    // set decompression order from compressed file
    double tup[3];
    orders->GetTuple(0, tup);
    this->mmax = static_cast<unsigned int>(tup[0]); 
    this->Mmax = static_cast<unsigned int>(0.5*(mmax+1)*(mmax+2));
    // build locator
    triloc->SetDataSet(polytri);
    triloc->BuildLocator(); 
  }
  
  // initialize decompressed image data
  this->imagedata = vtkSmartPointer<vtkImageData>::New();
  int dims[3]; 
  dims[0] = bounds[1]+1; dims[1] = bounds[3]+1; dims[2] = 1;
  double origin[3] = {0, 0, 0};
  imagedata->SetDimensions(dims);
  imagedata->SetOrigin(origin);
  imagedata->SetExtent  ( (int) bounds[0], 
                          (int) bounds[1], 
                          (int) bounds[2], 
                          (int) bounds[3],
                          (int) bounds[4],
                          (int) bounds[5] );

  this->pixels = imagedata->GetPoints();
  this->a = 0.5; this->b = 0.5; this->c = 0.5;

  this->colors = vtkSmartPointer<vtkUnsignedShortArray>::New();
  this->colors->SetName("colors");
  this->colors->SetNumberOfComponents(3);
  this->colors->SetNumberOfTuples(imagedata->GetNumberOfPoints());
}

void Decompressor::writeImage ( const std::string& pref,
                                const std::string& ext )
{
  imagedata->GetPointData()->AddArray(colors);
  imagedata->GetPointData()->SetActiveScalars("colors"); 
  std::string fname = pref + ext;

  vtkSmartPointer<vtkImageWriter> writer;
  if (ext == ".jpg" || ext == ".jpeg")
  {
    vtkSmartPointer<vtkJPEGWriter> jpgwriter = 
      vtkSmartPointer<vtkJPEGWriter>::New();
    writer = jpgwriter; 
  }
  else if (ext == ".png")
  {
    vtkSmartPointer<vtkPNGWriter> pngwriter = 
      vtkSmartPointer<vtkPNGWriter>::New();
    writer = pngwriter;
  }
  else if (ext == ".tiff")
  {
    vtkSmartPointer<vtkTIFFWriter> tiffwriter = 
      vtkSmartPointer<vtkTIFFWriter>::New();
    writer = tiffwriter;
  }
  else if (ext == ".ppm")
  {
    vtkSmartPointer<vtkPNMWriter> ppmwriter = 
      vtkSmartPointer<vtkPNMWriter>::New();
    writer = ppmwriter; 
  }

  //vtkSmartPointer<vtkAttributeSmoothingFilter> smoother 
  //  = vtkSmartPointer<vtkAttributeSmoothingFilter>::New();
  //smoother->SetInputData(imagedata);
  //smoother->SetSmoothingStrategyToAllButBoundary();
  //smoother->SetNumberOfIterations(5);
  //smoother->SetRelaxationFactor(0.1); 
  //smoother->Update();
 
  vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
  cast->SetInputData(imagedata);
  //cast->SetInputData(smoother->GetOutput());
  cast->SetOutputScalarTypeToUnsignedChar();
  cast->Update(); 
  
  vtkSmartPointer<vtkImageData> imagecast = cast->GetOutput();
  writer->SetFileName(fname.c_str());
  writer->SetInputData(imagecast);
  writer->Write(); 
}

void Decompressor::decompressChannel_help ( unsigned int channel,
                                            vtkSmartPointer<vtkPolyData> polytri,
                                            vtkSmartPointer<vtkDoubleArray> coeffs,
                                            vtkSmartPointer<vtkIntArray> offsets,
                                            vtkSmartPointer<vtkUnsignedIntArray> orders )
{
  double* cimg;
  if (not useMultiChannel)
  {
    cimg = (double*) calloc(Mmax, sizeof(double)); 
  }
  else
  {
    // build locator on image triangulation
    triloc->Initialize();
    triloc->SetDataSet(polytri);
    triloc->BuildLocator(); 
  }
  
  // build pix-tri map (pixInTri[j] are the pixels indices in triangle j)
  std::map<int, std::vector<int>> pixInTri;
  int triInd;
  for (int ipt = 0; ipt < pixels->GetNumberOfPoints(); ++ipt)
  {
    triInd = triloc->FindCell(pixels->GetPoint(ipt));
    if (triInd == - 1) { std::cout << "not in cell " << ipt << std::endl; }
    pixInTri[triInd].push_back(ipt);
  } 

  // get max num pix in tri
  auto it = pixInTri.begin(); 
  unsigned int maxpix = 0;
  while (it != pixInTri.end())
  {
    if (it->second.size() > maxpix)
    {
      maxpix = it->second.size();
    }
    ++it;
  }

  // allocate buffer for pixels in tri to ref 
  double* rspix = (double*) calloc(maxpix * 2, sizeof(double));  
  double* img = (double*) calloc(maxpix, sizeof(double));
  // large buffer jacobi interp op
  jPoly<double>* Pm = new jPoly<double>(maxpix, mmax, a, b, c, 1); 
  
  // iterate over triangles
  double v1t[3], v2t[3], v3t[3], tmpwts[3], r[3], dist2;
  int offset, npix, subid;
  unsigned int m, M;
  unsigned short color;
  double cimg_dub;
  for (it = pixInTri.begin(); it != pixInTri.end(); ++it)
  {
    // get tri verts
    vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
    polytri->GetCellPoints(it->first, ptids);
    polytri->GetPoints()->GetPoint(ptids->GetId(0), v1t);
    polytri->GetPoints()->GetPoint(ptids->GetId(1), v2t);
    polytri->GetPoints()->GetPoint(ptids->GetId(2), v3t);
    
    // get pixels inside current triangle
    // and get position in reference
    npix = it->second.size();
    for (unsigned int i = 0; i < npix; ++i)
    { 
      polytri->GetCell(it->first)->
        EvaluatePosition(pixels->GetPoint(it->second[i]), 
                                          nullptr, 
                                          subid, r, 
                                          dist2, tmpwts);
      rspix[i] = r[0]; rspix[i+npix] = r[1];

    }

    // compute jpoly interp matrix
    Pm->Nx = npix;
    Pm->computeV(rspix, rspix+npix);
    
    if (useMultiChannel)
    {
      // get order and offset
      offset = offsets->GetTuple1(it->first);
      m = orders->GetTuple1(it->first); 
      M = static_cast<unsigned int>(0.5 * (m + 1) * (m + 2));
      // offset to correct channel coeffs
      cimg = coeffs->GetPointer(offset);      
      // evaluate image over ref tri;
      cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                    npix, M, 1.0, Pm->V, npix, 
                    cimg, 1, 0.0, img, 1);
      
      // save the image data in channel
      for (unsigned int i = 0; i < npix; ++i)
      {
        // clip out of bounds colors
        if (img[i] < 0) { img[i] = 0; }
        if (img[i] > 255) { img[i] = 255; }
        color = static_cast<unsigned short>(std::round(img[i])); 
        colors->SetComponent(it->second[i], channel, color);
      }
    }
    else
    {
      for (unsigned int ichannel = 0; ichannel < 3; ++ichannel)
      {
        // get order and offset
        offset = offsets->GetComponent(it->first, ichannel);
        m = mmax;
        M = Mmax; 
        // copy correct channel coeffs
        for (unsigned int i = 0; i < M; ++i)
        {
          cimg_dub = coeffs->GetComponent(offset+i, ichannel);
          cimg[i] = static_cast<double>(half(cimg_dub));
        }
        // evaluate image over ref tri;
        cblas_dgemv ( CblasColMajor, CblasNoTrans, 
                      npix, M, 1.0, Pm->V, npix, 
                      cimg, 1, 0.0, img, 1);
        
        // save the image data in channel
        for (unsigned int i = 0; i < npix; ++i)
        {
          if (img[i] < 0) { img[i] = 0; }
          if (img[i] > 255) { img[i] = 255; }
          color = static_cast<unsigned short>(std::round(img[i])); 
          colors->SetComponent(it->second[i], ichannel, color);
        }
      }
    }
    
    if (!(it->first % 1000)) 
    { 
      std::cout << "Progress : " << std::setprecision(2) 
                << ((double)it->first / (double) pixInTri.size()) * 100 
                << "%\n";
    } 
  }

  free(rspix);
  free(img);
  delete Pm;
  if (not useMultiChannel)
  {
    free(cimg);
  }

}

void Decompressor::decompressChannel  ( unsigned int channel )
{
  vtkSmartPointer<vtkPolyData> polytri;
  vtkSmartPointer<vtkDoubleArray> coeffs;
  vtkSmartPointer<vtkIntArray> offsets;
  vtkSmartPointer<vtkUnsignedIntArray> orders;

  if (useMultiChannel)
  {
    switch (channel)
    { 
      case 0:
      {
        polytri = polytri1;
        coeffs = coeffs1;
        offsets = offsets1;
        orders = orders1;
        break;
      }
      case 1:
      {
        polytri = polytri2;
        coeffs = coeffs2;
        offsets = offsets2;
        orders = orders2;
        break;
      }
      case 2:
      {
        polytri = polytri3;
        coeffs = coeffs3;
        offsets = offsets3;
        orders = orders3;
        break;
      }
      default:
      {
        std::cerr << "channel must be 0,1,2\n"; 
        exit(1);
      }
    }

  }
  else
  {
    polytri = this->polytri;
    coeffs = this->coeffs;
    offsets = this->offsets;
    orders = this->orders;
  }
  decompressChannel_help  ( channel, polytri, coeffs, offsets, orders );
}

void Decompressor::run()
{
  if (useMultiChannel)
  {
    decompressChannel(0);
    std::cout << "decompressed channel 1" << std::endl;
    decompressChannel(1);
    std::cout << "decompressed channel 2" << std::endl;
    decompressChannel(2); 
    std::cout << "decompressed channel 3" << std::endl;
  }
  else
  {
    decompressChannel(0);
    std::cout << "decompressed all channels" << std::endl;
  }
}


void writeVTKLegacy(vtkSmartPointer<vtkPolyData> polytri, const char* ofname)
{
  // write mesh
  vtkSmartPointer<vtkDataSetWriter> writer
    = vtkSmartPointer<vtkDataSetWriter>::New();
  writer->SetFileName(ofname);
  writer->SetInputData(polytri);  
  writer->Write();
}

void writeTri(vtkSmartPointer<vtkPolyData> polytri, const char* pref)
{

  std::string prefname = pref;
  std::string ptsname = prefname + ".pts";
  std::string cellname = prefname + ".cells";

  vtkSmartPointer<vtkPoints> points = polytri->GetPoints();
  int npts = points->GetNumberOfPoints();
  int ncells = polytri->GetNumberOfCells();
  int szpts = npts*2;
  int szcells = ncells*3;
  float* ptsbuffer = (float*) calloc(szpts, sizeof(float));
  half* ptsbuffer_half = (half*) calloc(szpts, sizeof(half));
  unsigned int* cellbuffer 
    = (unsigned int*) calloc(szcells, sizeof(unsigned int)); 
  double pt[3];
  for (int ipt = 0; ipt < npts; ++ipt)
  {
    points->GetPoint(ipt, pt);
    for (int j = 0; j < 2; ++j)
    {
      ptsbuffer[ipt + npts*j] = static_cast<float>(pt[j]); 
      ptsbuffer_half[ipt + npts*j] = half(pt[j]); 
    }
  }

  float diff = 0, fltptsnrm = 0;
  for (unsigned int ipt = 0; ipt < npts*2; ++ipt)
  {
    diff += std::pow(ptsbuffer[ipt] - static_cast<float>(ptsbuffer_half[ipt]), 2.0);
    fltptsnrm += std::pow(ptsbuffer[ipt], 2.0);
  }
  std::cout << "Absolute error: pts float vs half " << std::sqrt(diff/fltptsnrm) << std::endl;  

  std::cout << "pts dim: " << npts << " x " << 3 << std::endl;

  vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
  cellPtIds->SetNumberOfIds(3);
  for (int icell = 0; icell < ncells; ++icell)
  {
    polytri->GetCellPoints(icell, cellPtIds);
    for (int j = 0; j < 3; ++j)
    {
      cellbuffer[icell + ncells*j] = static_cast<unsigned int>(cellPtIds->GetId(j));
    }
  }
  
  std::cout << "cell dim: " << ncells << " x " << 3 << std::endl;
 
  FILE *pfile, *cfile;
  pfile = fopen(ptsname.c_str(), "wb");
  cfile = fopen(cellname.c_str(), "wb");

  size_t nelmpts = fwrite(ptsbuffer_half, sizeof(half), szpts, pfile);  
  size_t nelmcells = fwrite(cellbuffer, sizeof(unsigned int), szcells, cfile);  
  if (nelmpts != szpts)
  {
    perror("Error writing to file\n");
    fclose(pfile);
    exit(1);
  }
  if (nelmcells != szcells)
  {
    perror("Error writing to file\n");
    fclose(cfile);
    exit(1);
  }
  fclose(pfile); fclose(cfile);
  free(ptsbuffer);
  free(ptsbuffer_half);
  free(cellbuffer);
}

void writeVTP(vtkSmartPointer<vtkPolyData> polytri, const char* ofname)
{
  // write mesh
  vtkSmartPointer<vtkXMLPolyDataWriter> writer
    = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(ofname);
  writer->SetInputData(polytri);  
  writer->SetDataModeToBinary();
  writer->Write();

}

void writeSTL(vtkSmartPointer<vtkPolyData> polytri, const char* ofname)
{
  vtkSmartPointer<vtkSTLWriter> stlwriter
    = vtkSmartPointer<vtkSTLWriter>::New();
  stlwriter->SetFileName(ofname);
  stlwriter->SetInputData(polytri);  
  stlwriter->Write();
}

void writeCoeffs(vtkSmartPointer<vtkPolyData> polytri, const char* ofname)
{
  vtkSmartPointer<vtkIntArray> offsets 
    = vtkIntArray::SafeDownCast(polytri->GetCellData()->GetAbstractArray(0));
  vtkSmartPointer<vtkUnsignedIntArray> orders 
    = vtkUnsignedIntArray::SafeDownCast(polytri->GetCellData()->GetAbstractArray(1));
  vtkSmartPointer<vtkDoubleArray> coeffs 
    = vtkDoubleArray::SafeDownCast(polytri->GetFieldData()->GetAbstractArray(0));
  
  int ncoeffs = coeffs->GetNumberOfTuples() * 3;
  // ncells x M x 3
  half* coeffsbuffer = (half*) calloc(ncoeffs, sizeof(half));
  double ms[3]; orders->GetTuple(0, ms); 
  int M = static_cast<int>(0.5 * (ms[0] + 1) * (ms[0] + 2));

  // load coeffs 
  int offset1, offset2, offset3;
  float c1, c2, c3;
  for (int icell = 0; icell < polytri->GetNumberOfCells(); ++icell)
  {
    offset1 = offsets->GetComponent(icell, 0);
    offset2 = offsets->GetComponent(icell, 1);
    offset3 = offsets->GetComponent(icell, 2);
    for (unsigned int i = 0; i < M; ++i)
    {
      c1 = coeffs->GetComponent(offset1+i, 0); 
      c2 = coeffs->GetComponent(offset2+i, 1);
      c3 = coeffs->GetComponent(offset3+i, 2);
      coeffsbuffer[3*(i + M*icell)] = half(c1);
      coeffsbuffer[1+3*(i + M*icell)] = half(c2);
      coeffsbuffer[2+3*(i + M*icell)] = half(c3);
    }
  }

  std::cout << "coeffs dim: " << polytri->GetNumberOfCells() 
            << " x " << M << " x " << 3 << std::endl;

  // write coeffs
  FILE* ofile = fopen(ofname,"wb");
  size_t nelmcoeffs = fwrite(coeffsbuffer, sizeof(half), ncoeffs, ofile);  
  if (nelmcoeffs!= ncoeffs)
  {
    perror("Error writing to file\n");
    fclose(ofile);
    exit(1);
  }
  fclose(ofile);
  free(coeffsbuffer);
}

vtkSmartPointer<vtkImageData> imageFFT  ( vtkSmartPointer<vtkImageData> imagedata)
{
  vtkSmartPointer<vtkImageFFT> fft = vtkSmartPointer<vtkImageFFT>::New();
  fft->SetInputData(imagedata);
  fft->Update();
  
  //vtkSmartPointer<vtkImageMagnitude> magnitude
  //  = vtkSmartPointer<vtkImageMagnitude>::New();
  //magnitude->SetInputData(fft->GetOutput());
  //magnitude->Update();
  
  //vtkSmartPointer<vtkImageFourierCenter> center 
  //  = vtkSmartPointer<vtkImageFourierCenter>::New();
  //center->SetInputData(magnitude->GetOutput());
  //center->Update();

  //vtkSmartPointer<vtkImageLogarithmicScale> lgsc
  //  = vtkSmartPointer<vtkImageLogarithmicScale>::New();
  ////lgsc->SetInputData(center->GetOutput());
  ////lgsc->SetInputData(magnitude->GetOutput());
  //lgsc->SetInputData(fft->GetOutput());
  //lgsc->SetConstant(15);
  //lgsc->Update();

  return fft->GetOutput(); 
}

vtkSmartPointer<vtkImageData> smoothImage ( vtkSmartPointer<vtkImageData> imagedata)
{
  vtkSmartPointer<vtkImageGaussianSmooth> smoother
    = vtkSmartPointer<vtkImageGaussianSmooth>::New();
  smoother->SetDimensionality(2);
  smoother->SetInputData(imagedata);
  smoother->SetStandardDeviations(3.0, 3.0);
  smoother->SetRadiusFactors(3, 3);
  smoother->Update();
  return smoother->GetOutput();
}

void writeSmoothImage ( vtkSmartPointer<vtkImageData> imagedata, const std::string& fname )
{
  vtkSmartPointer<vtkImageGaussianSmooth> smoother
    = vtkSmartPointer<vtkImageGaussianSmooth>::New();
  smoother->SetInputData(imagedata);
  smoother->Update();

  
  vtkSmartPointer<vtkPNMWriter> ppmwriter = 
    vtkSmartPointer<vtkPNMWriter>::New();
  ppmwriter->SetFileName(fname.c_str());
  ppmwriter->SetInputData(smoother->GetOutput());
  ppmwriter->Write(); 
  
}

vtkSmartPointer<vtkImageData> readImage ( const std::string& pref,
                                          const std::string& ext )
{
  std::string fname = pref + ext;
  vtkSmartPointer<vtkImageReader2> reader;
  if (ext == ".jpg" || ext == ".jpeg")
  {
    vtkSmartPointer<vtkJPEGReader> jpgreader = 
      vtkSmartPointer<vtkJPEGReader>::New();
    reader = jpgreader; 
  }
  else if (ext == ".png")
  {
    vtkSmartPointer<vtkPNGReader> pngreader = 
      vtkSmartPointer<vtkPNGReader>::New();
    reader = pngreader;
  }
  else if (ext == ".tiff")
  {
    vtkSmartPointer<vtkTIFFReader> tiffreader = 
      vtkSmartPointer<vtkTIFFReader>::New();
    reader = tiffreader;
  }
  else if (ext == ".ppm")
  {
    vtkSmartPointer<vtkPNMReader> ppmreader = 
      vtkSmartPointer<vtkPNMReader>::New();
    reader = ppmreader; 
  }
  reader->SetFileName(fname.c_str());
  reader->Update();
  return reader->GetOutput();
}

Triangulator* jcompress_triangulate ( const char* fname, 
                                      unsigned int nSamp,
                                      unsigned int nRuns,
                                      unsigned int mtarget,
                                      bool useMultiChannel,
                                      bool viz )
{
  // jacobi poly
  double a = 0.5, b = 0.5, c = 0.5;
  // triangulation quad settings
  unsigned int N = 55;
  unsigned int m = 6;
  //unsigned int N = 496;
  //unsigned int m = mtarget;

  // Read the image
  std::string fWithExt(fname);
  std::filesystem::path p(fWithExt);
  std::string fout = p.stem().string();
  std::string fWithoutExt = p.parent_path().string() + "/" + fout;
  vtkSmartPointer<vtkImageData> imagedata = readImage(fWithoutExt.c_str(), p.extension().c_str()); 
  // View the image
  if (viz)
  {
    vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
    vtkSmartPointer<vtkImageViewer2> imageViewer = vtkSmartPointer<vtkImageViewer2>::New();
    imageViewer->SetInputData(imagedata);
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    imageViewer->SetupInteractor(renderWindowInteractor);
    imageViewer->Render();
    imageViewer->GetRenderer()->ResetCamera();
    imageViewer->GetRenderer()->SetBackground(
        colors->GetColor3d("DarkSlateGray").GetData());
    imageViewer->GetRenderWindow()->SetWindowName(fname);
    imageViewer->Render();
    renderWindowInteractor->Start();
  }
  // create cubic interpolator on image
   vtkSmartPointer<vtkImageInterpolator> interpolator = 
    vtkSmartPointer<vtkImageInterpolator>::New();
  interpolator->Initialize(imagedata);
  interpolator->SetInterpolationModeToCubic();
  
  // triangulate
  Triangulator* T = new Triangulator  ( N, m, a, b, c, nSamp, nRuns, 
                                        "xtri_N55_n9_M91_m12.txt",
                                        "ytri_N55_n9_M91_m12.txt",            
                                        "wtri_N55_n9_M91_m12.txt",
                                        imagedata, interpolator, 
                                        useMultiChannel );
   
  //Triangulator* T = new Triangulator  ( N, m, a, b, c, nSamp, nRuns, 
  //                                      "xtri_N496_n30_M1378_m51.txt",                   
  //                                      "ytri_N496_n30_M1378_m51.txt",                               
  //                                      "wtri_N496_n30_M1378_m51.txt",                   
  //                                      imagedata, interpolator, 
  //                                      useMultiChannel );
  //Triangulator* T = new Triangulator  ( N, m, a, b, c, nSamp, 5, 
  //                                      "xtri_N55_n9_M91_m12.txt",
  //                                      "ytri_N55_n9_M91_m12.txt",            
  //                                      "wtri_N55_n9_M91_m12.txt",
  //                                      imagedata, interpolator);
  T->Mtarget = static_cast<unsigned int>(0.5 * (mtarget + 1) * (mtarget + 2));
  T->run();
  return T;
}

double jcompress  ( Triangulator* T,
                    const char* fname, 
                    unsigned int nSamp,
                    unsigned int nRuns,
                    unsigned int order,
                    bool useMultiChannel,
                    bool viz )
{
  // compress
  Compressor* C = new Compressor(T, order);
  C->run();
  double totalBytes; 
  // write grids with compressed data
  std::string fWithExt(fname);
  std::filesystem::path p(fWithExt);
  std::string fout = p.stem().string();
  std::string fWithoutExt = p.parent_path().string() + "/" + fout;
  if (useMultiChannel)
  {
    std::string fout1 = fout + "_" + std::to_string(order) + "_" + "1.vtp";
    std::string fout2 = fout + "_" + std::to_string(order) + "_" + "2.vtp";
    std::string fout3 = fout + "_" + std::to_string(order) + "_" + "3.vtp";
    writeVTP(T->polytri1, fout1.c_str());
    writeVTP(T->polytri2, fout2.c_str());
    writeVTP(T->polytri3, fout3.c_str());
  }
  else
  {
    std::string fout1 = fout + "_" + std::to_string(nSamp) + "_" 
                                   + std::to_string(nRuns) + "_"
                                   +  std::to_string(order) + ".vtp";
    std::string fout2 = fout + "_" + std::to_string(nSamp) + "_" 
                                   + std::to_string(nRuns) + "_"
                                   + std::to_string(order) + ".stl";
    std::string fout3 = fout + "_" + std::to_string(nSamp) + "_" 
                                   + std::to_string(nRuns) + "_"
                                   + std::to_string(order) + ".coeffs";
    std::string fout4 = fout + "_" + std::to_string(nSamp) + "_" 
                                   + std::to_string(nRuns) + "_"
                                   + std::to_string(order) + ".bench";
    writeVTP(C->polytri, fout1.c_str());
    writeSTL(C->polytri, fout2.c_str());
    writeCoeffs(C->polytri, fout3.c_str());

  
    std::string command = "./lz.sh " + fout4 + " " + fout2 + " " + fout3;

    FILE* fp = popen(command.c_str(), "r");
    if (fp == NULL)
    {
      printf("Failed to run script\n");
    } 
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), fp) != NULL){}
    if (buffer[strlen(buffer) - 1] == '\n')
    {
      buffer[strlen(buffer) - 1] == '\0';
    } 
    totalBytes = atof(buffer);
    pclose(fp);
    // double totalBytes = C->totalBytes;
  }
   
  // clean
  delete C;
  
  return totalBytes;
}



double jcompress  ( const char* fname, 
                    unsigned int nSamp,
                    unsigned int nRuns,
                    unsigned int order,
                    bool useMultiChannel,
                    bool viz )
{
  // triangulate
  Triangulator* T = jcompress_triangulate ( fname, nSamp, nRuns, order,
                                            useMultiChannel, false );

  // compress
  Compressor* C = new Compressor(T, order);
  C->run();
  double totalBytes; 
  // write grids with compressed data
  std::string fWithExt(fname);
  std::filesystem::path p(fWithExt);
  std::string fout = p.stem().string();
  std::string fWithoutExt = p.parent_path().string() + "/" + fout;
  if (useMultiChannel)
  {
    std::string fout1 = fout + "_" + std::to_string(order) + "_" + "1.vtp";
    std::string fout2 = fout + "_" + std::to_string(order) + "_" + "2.vtp";
    std::string fout3 = fout + "_" + std::to_string(order) + "_" + "3.vtp";
    writeVTP(T->polytri1, fout1.c_str());
    writeVTP(T->polytri2, fout2.c_str());
    writeVTP(T->polytri3, fout3.c_str());
  }
  else
  {
    std::string fout1 = fout + "_" + std::to_string(order) + ".vtp";
    std::string fout2 = fout + "_" + std::to_string(order) + ".pts";
    std::string fout5 = fout + "_" + std::to_string(order) + ".cells";
    std::string fout3 = fout + "_" + std::to_string(order) + ".coeffs";
    std::string fout4 = fout + "_" + std::to_string(nSamp) + "_" 
                                   + std::to_string(nRuns) + "_"
                                   + std::to_string(order) + ".bench";

    std::string fout00 = fout + "_" + std::to_string(order);
    writeVTP(C->polytri, fout1.c_str());
    writeTri(T->polytri, fout00.c_str());
    writeCoeffs(C->polytri, fout3.c_str()); 

    std::string command = "./lz.sh " + fout4 + " " + fout2 + " " + fout3 + " " + fout5;

    FILE* fp = popen(command.c_str(), "r");
    if (fp == NULL)
    {
      printf("Failed to run script\n");
    } 
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), fp) != NULL){}
    if (buffer[strlen(buffer) - 1] == '\n')
    {
      buffer[strlen(buffer) - 1] == '\0';
    } 
    totalBytes = atof(buffer);
    std::cout << "Compressed all channels into " << totalBytes << " bytes\n";
    pclose(fp);
  }


   
  // clean
  delete T;
  delete C;
  
  return totalBytes;
}

void jdecompress  ( bool useMultiChannel,
                    const char* fmt, 
                    const char* channel1,
                    const char* channel2,
                    const char* channel3 )
{
  Decompressor* D;
  if (useMultiChannel)
  {
    D = new Decompressor ( useMultiChannel, channel1, channel2, channel3 );
  }
  else
  {
    D = new Decompressor ( useMultiChannel, channel1 ); 
  }
  D->run();
  unsigned int order = D->mmax;
  std::string fWithExt(channel1);
  std::filesystem::path p(fWithExt);
  std::string ext(fmt); ext = "." + ext;
  std::string fout = p.stem().string() + "_deco";
  D->writeImage(fout, ext);
  delete D;
}

double ssim ( vtkSmartPointer<vtkImageData> img1,
              vtkSmartPointer<vtkImageData> img2 )
{

  vtkSmartPointer<vtkImageSSIM> ssim = vtkSmartPointer<vtkImageSSIM>::New();
  ssim->SetInputToRGB();
  ssim->SetInputData(img1); 
  ssim->SetImageData(img2);  
  ssim->Update();
  
  vtkSmartPointer<vtkImageData> ssim_output =  ssim->GetOutput();
  double tup[3], ave = 0;
  int numtuples = ssim_output->GetPointData()->GetNumberOfTuples();
  vtkSmartPointer<vtkDataArray> scalars = ssim_output->GetPointData()->GetScalars();
  for (unsigned int i = 0; i < numtuples; ++i)
  {
    scalars->GetTuple(i, tup);
    ave += (tup[0] + tup[1] + tup[2]) / 3;
  }
  ave /= numtuples; 
 
  std::cout << "ssim average: " << ave << std::endl;
  return ave;
} 

double ssim ( const char* F1WithExt, const char* F2WithExt )
{
  std::string f1WithExt(F1WithExt);
  std::filesystem::path p1(f1WithExt);
  std::string f1 = p1.stem().string();
  std::string ext1 = p1.extension();
  std::string f1WithoutExt;
  if (not p1.parent_path().string().empty())
  {
    f1WithoutExt = p1.parent_path().string() + "/" + f1; 
  }
  else
  {
    f1WithoutExt = f1;
  }

  std::string f2WithExt(F2WithExt);
  std::filesystem::path p2(f2WithExt);
  std::string f2 = p2.stem().string();
  std::string ext2 = p2.extension();
  std::string f2WithoutExt;
  if (not p2.parent_path().string().empty())
  {
    f2WithoutExt = p2.parent_path().string() + "/" + f2;
  }
  else
  {
    f2WithoutExt = f2;
  }
  return ssim ( readImage(f1WithoutExt, ext1),  
                readImage(f2WithoutExt, ext2) );

}


/************************************ 

  m-adaptive compression code 
   
    while (m <= mmax)
    {
      unsigned int M = static_cast<unsigned int>(0.5 * (m + 1) * (m + 2));
      unsigned int Mprev = static_cast<unsigned int>(0.5 * m * (m + 1));
      Pm->computeCoeffsM(interpc, R, S, W, cimg, m, hasweights);
      hasweights = true;
      blknorm = cblas_dnrm2(M-Mprev, cimg+Mprev, 1);
      if (blknorm < blknormtol)
      { 
        for (unsigned int j = 0; j < M; ++j)
        {
          coeffs->InsertNextValue(cimg[j]);
        }
        offsets->SetValue(i, offset);
        orders->SetValue(i, m);
        offset += M;
        break;
      }
      m += 2; 
    }
    // if blknormtol never statisfied
    if (m > mmax)
    {
      for (unsigned int j = 0; j < Mmax; ++j)
      {
        coeffs->InsertNextValue(cimg[j]);
      }
      offsets->SetValue(i, offset);
      orders->SetValue(i, mmax);
      offset += Mmax;
    }


  pix to tri map checks

      if (rspix[i] < 0 || rspix[i] > 1 || rspix[i+npix] < 0 || rspix[i+npix] > 1-rspix[i])
      {
        std::cerr << "outside of triangle" << std::endl;
        std::cerr << rspix[i] << " " << rspix[i+npix] << " " 
                  << ret << std::setprecision(17) 
                  << rspix[i+npix]-(1-rspix[i]) << std::endl;
      }
  
  dgelsd smoothing checks
      double res = 0;
      for (unsigned int i = Mmax; i < N+3; ++i)
      {
        res += cRHS[i] * cRHS[i];
      }
      std::cout << res << std::endl;
      smoothCoeffs_help  ( 0, icell, nleg, lPm, bPm, hPm,
                           cRHS, cimg11, imgbnd, imgbnd1,
                           pcoords10, pcoords11,
                           offset11, subid, wts, dist2,
                           cellPtIds, neighborCellIds,
                           neighbors, 2, true);             

  
        // populate cRHS for channel 0
        cblas_dgemv ( CblasColMajor, CblasNoTrans,
                      N, Mmax, 1.0, Pm->V, N,
                      cimg1, 1, 0.0, cRHS, 1 );      
        cRHS[N] = sum1;
        cRHS[N+1] = sum2;
        cRHS[N+2] = sum3;
        cblas_dcopy ( (N+3)*Mmax, cMat, 1, cMat_cpy, 1 );
        if (LAPACKE_dgelsd  ( LAPACK_COL_MAJOR, N+3, Mmax, 1, cMat_cpy, 
                              N+3, cRHS, N+3, S, -1.0, rank ) )
        {
          std::cerr << "ERROR: Lapack dgelsd: Pseudoinverse" << std::endl;
        }
        else
        {
          for (unsigned int i = 0; i < Mmax; ++i)
          {
            coeffs->SetComponent(offset1+i, 0, cRHS[i]);
          }
        }
        // populate cRHS for channel 1
        cblas_dgemv ( CblasColMajor, CblasNoTrans,
                      N, Mmax, 1.0, Pm->V, N,
                      cimg2, 1, 0.0, cRHS, 1 );      
        cRHS[N] = sum11;
        cRHS[N+1] = sum21;
        cRHS[N+2] = sum31;
        cblas_dcopy ( (N+3)*Mmax, cMat, 1, cMat_cpy, 1 );
        if (LAPACKE_dgelsd  ( LAPACK_COL_MAJOR, N+3, Mmax, 1, cMat_cpy, 
                              N+3, cRHS, N+3, S, -1.0, rank ) )
        {
          std::cerr << "ERROR: Lapack dgelsd: Pseudoinverse" << std::endl;
        }
        else
        {
          for (unsigned int i = 0; i < Mmax; ++i)
          {
            coeffs->SetComponent(offset2+i, 1, cRHS[i]);
          }
        }
        // populate cRHS for channel 2
        cblas_dgemv ( CblasColMajor, CblasNoTrans,
                      N, Mmax, 1.0, Pm->V, N,
                      cimg3, 1, 0.0, cRHS, 1 );      
        cRHS[N] = sum12;
        cRHS[N+1] = sum22;
        cRHS[N+2] = sum32;
        cblas_dcopy ( (N+3)*Mmax, cMat, 1, cMat_cpy, 1 );
        if (LAPACKE_dgelsd  ( LAPACK_COL_MAJOR, N+3, Mmax, 1, cMat_cpy, 
                              N+3, cRHS, N+3, S, -1.0, rank ) )
        {
          std::cerr << "ERROR: Lapack *gelsd: Pseudoinverse" << std::endl;
        }
        else
        {
          for (unsigned int i = 0; i < Mmax; ++i)
          {
            coeffs->SetComponent(offset3+i, 2, cRHS[i]);
          }
        }

struct ArrayHash {
    std::size_t operator()(const std::array<double, 2>& arr) const {
        std::size_t h1 = std::hash<double>{}(arr[0]);
        std::size_t h2 = std::hash<double>{}(arr[1]);
        return h1 ^ (h2 << 1); // Simple way to combine hashes
    }
};

// Define a custom equality operator for std::array<double, 2>
struct ArrayEqual {
    bool operator()(const std::array<double, 2>& lhs, const std::array<double, 2>& rhs) const {
        return lhs[0] == rhs[0] && lhs[1] == rhs[1];
    }
};


*************************************/
