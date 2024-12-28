#include <jcodec.hh>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkInformation.h>
#include <vtkCellTreeLocator.h>
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
#include <map>
#include <vector>
#include <timer.hh>
#include <filesystem>

vtkSmartPointer<vtkDelaunay2D> triangulateUniform ( int* dims, 
                                                    double* origin,
                                                    unsigned int nSamp )
{
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
      ptId += 1;  
    }
  }
 
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
    // read the image
    this->imagedata = _imagedata;
    int dims[3]; imagedata->GetDimensions(dims);
    double origin[3]; imagedata->GetOrigin(origin);
    // get the interpolator
    this->interpolator = _interpolator;
  
    // triangulate uniform subsampled grid on image
    vtkSmartPointer<vtkDelaunay2D> triangulator = 
      triangulateUniform(dims, origin, nSamp);
    polytri = triangulator->GetOutput();
    std::cout << "Num cells on subsampled grid: " 
              << polytri->GetNumberOfCells() << std::endl;
    if (useMultiChannel)
    { 
      this->polytri1 = vtkSmartPointer<vtkPolyData>::New();
      this->polytri2 = vtkSmartPointer<vtkPolyData>::New();
      this->polytri3 = vtkSmartPointer<vtkPolyData>::New();
      this->polytri1->DeepCopy(polytri);
      this->polytri2->DeepCopy(polytri);
      this->polytri3->DeepCopy(polytri);
    }
  }

Triangulator::~Triangulator()
{
  if (Pm) { delete Pm; }
  if (Pmx) { delete Pmx; }
  if (Pmy) { delete Pmy; }
  if (Habc) { free(Habc); }
  if (Ha1bc1) { free(Ha1bc1); }
  if (Hab1c1) { free(Hab1c1); }
  if (Dx) { free(Dx); }
  if (Dy) { free(Dy); }
  if (X) { free(X);  } 
  if (Y) { free(Y);  } 
  if (W) { free(W);  }
  if (cdimg) { free(cdimg);  }
  if (dimgr) { free(dimgr);  }
  if (dimgt) { free(dimgt);  } 
  if (interpc) { free(interpc); }
}

double Triangulator::triangulateEntropy_help  ( double* intgn, unsigned int channel,
                                                vtkSmartPointer<vtkPolyData> polytri )
{
  double v1t[3], v2t[3], v3t[3], tmpwts[3];
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  double invIxe[4], gradnorm;
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
    // inverse of deformation map 
    invIxe[0] = v3t[1] - v1t[1];
    invIxe[1] = -(v2t[1] - v1t[1]);
    invIxe[2] = -(v3t[0] - v1t[0]); 
    invIxe[3] = v2t[0] - v1t[0];
    // coefficients of interpolated data in channel
    Pm->computeCoeffs(interpc, X, Y, W, cimg);
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
      gradnorm += (dimgt[j] * dimgt[j] + dimgt[j+N] * dimgt[j+N]) * W[j] / 2.0;
    }
    intgn[i] = gradnorm;
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

  double ave = 0;//, ave1 = 0, ave2 = 0;
  #pragma omp simd reduction(+:ave)
  for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  {
    ave += intgn[i];
  } 
  
  ave /= polytri->GetNumberOfCells();

  return ave;
}

vtkSmartPointer<vtkDelaunay2D> Triangulator::triangulateEntropy ( double* intgn, unsigned int channel  )
{
  vtkSmartPointer<vtkPolyData> polytri;
  double ave;
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
    ave = triangulateEntropy_help ( intgn, channel, polytri );
  }
  else
  {
    polytri = this->polytri;
    double* intgn1 = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
    double* intgn2 = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
    ave = ( triangulateEntropy_help ( intgn,  0, polytri ) +
            triangulateEntropy_help ( intgn1, 1, polytri ) +
            triangulateEntropy_help ( intgn2, 2, polytri )  ) / 3.0;
    for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
    {
      intgn[i] += intgn1[i] + intgn2[i];
      intgn[i] /= 3.0;
    }
    free(intgn1);
    free(intgn2);
  }

  vtkSmartPointer<vtkIncrementalOctreePointLocator> locator = 
    vtkSmartPointer<vtkIncrementalOctreePointLocator>::New();
  
  locator->SetDataSet(polytri);
  locator->BuildLocator();
  double v1r[3] = {0.5, 0, 0};
  double v2r[3] = {0.5, 0.5, 0};
  double v3r[3] = {0, 0.5, 0};
  double v1t[3], v2t[3], v3t[3], tmpwts[3];
  vtkIdType ptid; int subid;
  for (int icell = 0; icell < polytri->GetNumberOfCells(); ++icell)
  {
    if (intgn[icell] > ave)
    {
      polytri->GetCell(icell)->EvaluateLocation(subid, v1r, v1t, tmpwts);
      polytri->GetCell(icell)->EvaluateLocation(subid, v2r, v2t, tmpwts);
      polytri->GetCell(icell)->EvaluateLocation(subid, v3r, v3t, tmpwts);
      locator->InsertUniquePoint(v1t, ptid);  
      locator->InsertUniquePoint(v2t, ptid);  
      locator->InsertUniquePoint(v3t, ptid);
    } 
  } 
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New(); 
  polydata->SetPoints(locator->GetLocatorPoints());
  vtkSmartPointer<vtkDelaunay2D> triangulator = vtkSmartPointer<vtkDelaunay2D>::New();
  triangulator->SetInputData(polydata);
  triangulator->Update();
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
      std::cout << "Refinement step: " << iRun << std::endl;
      std::cout << "Num cells on each grid: "
                << polytri1->GetNumberOfCells() << ", "
                << polytri2->GetNumberOfCells() << ", "
                << polytri3->GetNumberOfCells() << "\n\n";
    }
  }
  else
  {
    for (unsigned int iRun = 0; iRun < nRuns; ++iRun)
    {
      double* intgn = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
      vtkSmartPointer<vtkDelaunay2D> triangulator = triangulateEntropy ( intgn, 0 );
      polytri = triangulator->GetOutput();
      free(intgn);
      std::cout << "Refinement step: " << iRun << std::endl;
      std::cout << "Num cells on grid: "
                << polytri->GetNumberOfCells() << "\n\n";
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
  }
  this->interpolator = T->interpolator;
  this->R = (double*) calloc(this->N, sizeof(double));
  this->S = (double*) calloc(this->N, sizeof(double));
  this->W = (double*) calloc(this->N, sizeof(double));
  this->interpc = (double*) calloc(this->N, sizeof(double));
  this->a = T->a; this->b = T->b; this->c = T->c;
  this->Pm = new jPoly<double>(N, morder, a, b, c, nthreads); 
  this->Mmax = Pm->Np;
  // storage buffer for img coefs
  this->cimg = (double*) calloc(Mmax, sizeof(double)); 
  readQuad(trix, triy, triw, N, R, S, W);
}

Compressor::~Compressor()
{
  if (Pm) { delete Pm; Pm = 0; }
  if (R) { free(R); R = 0; }
  if (S) { free(S); S = 0; }
  if (W) { free(W); W = 0; }
  if (interpc) { free(interpc); interpc = 0; }
  if (cimg) { free(cimg); cimg = 0; }
}

void Compressor::compressChannel_help ( unsigned int channel,
                                        vtkSmartPointer<vtkPolyData> polytri,
                                        vtkSmartPointer<vtkDoubleArray> coeffs,
                                        vtkSmartPointer<vtkIntArray> offsets,
                                        vtkSmartPointer<vtkUnsignedIntArray> orders )
{
  double v1t[3], v2t[3], v3t[3], tmpwts[3];
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  double invIxe[4], gradnorm;
  
  vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
  unsigned int m = this->morder; 
  int offset = 0;
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
    // get tri verts
    polytri->GetCellPoints(i, ptids);
    polytri->GetPoints()->GetPoint(ptids->GetId(0), v1t);
    polytri->GetPoints()->GetPoint(ptids->GetId(1), v2t);
    polytri->GetPoints()->GetPoint(ptids->GetId(2), v3t);
    unsigned int M = static_cast<unsigned int>(0.5 * (m + 1) * (m + 2));
    Pm->computeCoeffs(interpc, R, S, W, cimg);
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
      for (unsigned int j = 0; j < M; ++j)
      {
        coeffs->InsertComponent(offset+j, channel, cimg[j]);
      }
      offsets->SetComponent(i, channel, offset);
      orders->SetComponent(i, channel, m);
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
    compressChannel_help ( channel, polytri, coeffs, offsets, orders );
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
    compressChannel_help ( 0, polytri, coeffs, offsets, orders );
    compressChannel_help ( 1, polytri, coeffs, offsets, orders );
    compressChannel_help ( 2, polytri, coeffs, offsets, orders );
  }
  
  polytri->GetFieldData()->AddArray(coeffs);
  polytri->GetCellData()->AddArray(offsets); 
  polytri->GetCellData()->AddArray(orders); 
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
    unsigned int M = static_cast<unsigned int>(0.5 * (morder+1)*(morder+2));
    /* (ntri*(mcoeff floats/channel)*(3 channels)*(4bytes/float) 
      + ntri*(3 indices)*(4 bytes/int) 
      + npts * (2 floats for x,y)*(4 bytes/float) 
    */
    totalBytes =  polytri->GetNumberOfCells() * (M * 12 + 12) +
                  polytri->GetNumberOfPoints() * 8; 
    std::cout << "compressed all channels into "
              << totalBytes / 1e6 << " MB\n";
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
  pixels->SetNumberOfPoints(dims[0] * dims[1]);
  int ptid = 0; double pt[3]; pt[2] = 0;
  for (unsigned int iy = 0; iy < dims[1]; ++iy)
  {
    pt[1] = (double) iy;
    for (unsigned int ix = 0; ix < dims[0]; ++ix)
    {
      pt[0] = (double) ix; 
      pixels->SetPoint(ptid, pt);
      ptid += 1;
    }
  }

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
  // build locator on image triangulation
  vtkSmartPointer<vtkCellTreeLocator> triloc =
    vtkSmartPointer<vtkCellTreeLocator>::New();
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
  double* img = (double*) calloc(maxpix, sizeof(double));
  // large buffer jacobi interp op
  jPoly<double>* Pm = new jPoly<double>(maxpix, mmax, a, b, c, 1); 
  
  // iterate over triangles
  double v1t[3], v2t[3], v3t[3], tmpwts[3], r[3], dist2;
  int offset, npix, subid;
  unsigned int m, M;
  unsigned short color;
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
      polytri->GetCell(it->first)->EvaluatePosition(pixels->GetPoint(it->second[i]), nullptr, subid, r, dist2, tmpwts);
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
          cimg[i] = coeffs->GetComponent(offset+i, ichannel);
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
    
    if (!(it->first % 200)) 
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


void writeVTP(vtkSmartPointer<vtkPolyData> polytri, const char* ofname)
{
  // write mesh
  vtkNew<vtkXMLPolyDataWriter> writer;
  writer->SetFileName(ofname);
  writer->SetInputData(polytri);  
  writer->SetDataModeToBinary();
  writer->Write();
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
                                      unsigned int order,
                                      bool useMultiChannel,
                                      bool viz )
{
  // jacobi poly
  double a = 0.5, b = 0.5, c = 0.5;
  // triangulation quad settings
  unsigned int N = 55;
  unsigned int m = 6;

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
    writeVTP(C->polytri, fout1.c_str());
  }

  double totalBytes = C->totalBytes;
   
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
                                            useMultiChannel, viz );

  // compress
  Compressor* C = new Compressor(T, order);
  C->run();
  
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
    writeVTP(C->polytri, fout1.c_str());
  }

  double totalBytes = C->totalBytes;
   
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

*************************************/
