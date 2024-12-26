#include<jcodec.hh>
#include<vtkCellData.h>
#include<vtkFieldData.h>

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

void readQuad(const std::string& trix, 
              const std::string& triy, 
              const std::string& triw,
              unsigned int N,
              double* X, double* Y, 
              double*  W)
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
                              unsigned int _nthreads )
  : N(_N), m(_m), a(_a), b(_b), c(_c), nSamp(_nSamp), nRuns(_nRuns),
    trix(_trix), triy(_triy), triw(_triw), 
    nthreads(_nthreads)
  {
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
    vtkSmartPointer<vtkPolyData> polytri = triangulator->GetOutput();
    std::cout << "Num cells on subsampled grid: " 
              << polytri->GetNumberOfCells() << std::endl;
    this->polytri1 = vtkSmartPointer<vtkPolyData>::New();
    this->polytri2 = vtkSmartPointer<vtkPolyData>::New();
    this->polytri3 = vtkSmartPointer<vtkPolyData>::New();
    this->polytri1->DeepCopy(polytri);
    this->polytri2->DeepCopy(polytri);
    this->polytri3->DeepCopy(polytri);
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

vtkSmartPointer<vtkDelaunay2D> Triangulator::triangulateEntropy ( double* intgn, unsigned int channel  )
{
  vtkSmartPointer<vtkPolyData> polytri;
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
    // inverse of incidence matrix 
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

  vtkSmartPointer<vtkIncrementalOctreePointLocator> locator = 
    vtkSmartPointer<vtkIncrementalOctreePointLocator>::New();
  
  locator->SetDataSet(polytri);
  locator->BuildLocator();
  ptids->Reset();
  double v1r[3] = {0.5, 0, 0};
  double v2r[3] = {0.5, 0.5, 0};
  double v3r[3] = {0, 0.5, 0};
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

Compressor::Compressor(Triangulator* T)
{
  this->polytri1 = T->polytri1; 
  this->polytri2 = T->polytri2; 
  this->polytri3 = T->polytri3;
  this->interpolator = T->interpolator;
  this->nthreads = T->nthreads; 
  this->R = (double*) calloc(this->N, sizeof(double));
  this->S = (double*) calloc(this->N, sizeof(double));
  this->W = (double*) calloc(this->N, sizeof(double));
  this->interpc = (double*) calloc(this->N, sizeof(double));
  this->a = T->a; this->b = T->b; this->c = T->c;
  this->Pm = new jPoly<double>(N, mmax, a, b, c, nthreads); 
  this->Mmax = Pm->Np;
  // storage buffer for img coefs
  this->cimg = (double*) calloc(Mmax, sizeof(double)); 
  readQuad(trix, triy, triw, N, R, S, W);
  this->coeffs1 = vtkSmartPointer<vtkDoubleArray>::New();
  this->coeffs2 = vtkSmartPointer<vtkDoubleArray>::New();
  this->coeffs3 = vtkSmartPointer<vtkDoubleArray>::New();
  this->offsets1 = vtkSmartPointer<vtkIntArray>::New();
  this->offsets2 = vtkSmartPointer<vtkIntArray>::New();
  this->offsets3 = vtkSmartPointer<vtkIntArray>::New();
  
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
void Compressor::compressChannel(unsigned int channel, double blknormtol)
{
  vtkSmartPointer<vtkPolyData> polytri;
  vtkSmartPointer<vtkDoubleArray> coeffs;
  vtkSmartPointer<vtkIntArray> offsets;
  switch (channel)
  { 
    case 0:
    {
      polytri = polytri1;
      coeffs = coeffs1;
      offsets = offsets1;
      break;
    }
    case 1:
    {
      polytri = polytri2;
      coeffs = coeffs2;
      offsets = offsets2;
      break;
    }
    case 2:
    {
      polytri = polytri3;
      coeffs = coeffs3;
      offsets = offsets3;
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
  coeffs->SetName("coeffs");
  coeffs->SetNumberOfComponents(1);

  double v1t[3], v2t[3], v3t[3], tmpwts[3];
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  double invIxe[4], gradnorm;
  
  vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
  unsigned int m;
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
    // starting with linear order, get coeffs 
    // of increasing order and stop when norm of 
    // last block is within blknormtol 
    m = 1; double blknorm; 
    while (m < mmax)
    {
      unsigned int M = static_cast<unsigned int>(0.5 * (m + 1) * (m + 2));
      unsigned int Mprev = static_cast<unsigned int>(0.5 * m * (m + 1));
      Pm->computeCoeffsM(interpc, R, S, W, cimg, m);
      blknorm = cblas_dnrm2(M-Mprev, cimg+Mprev, 1);
      if (blknorm < blknormtol)
      { 
        for (unsigned int j = 0; j < M; ++j)
        {
          coeffs->InsertNextValue(cimg[j]);
        }
        offsets->SetValue(i, offset);
        offset += M; 
        break;
      } 
      m += 1; 
    }
  }
  polytri->GetFieldData()->AddArray(coeffs);
  polytri->GetCellData()->AddArray(offsets); 
}

void Compressor::run(double ctol)
{ 
  compressChannel(0, ctol);
  compressChannel(1, ctol);
  compressChannel(2, ctol);
} 
