#include <vtkImageViewer2.h>
#include <vtkJPEGReader.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkImageToStructuredGrid.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkStructuredGrid.h>
#include <vtkUniformGrid.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkImageInterpolator.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vtkIncrementalOctreePointLocator.h>
#include <vtkTriangle.h>
#include<cblas.h>
#include<lapacke.h>
#include<dMat.hh>
#include<jPoly.hh>

vtkSmartPointer<vtkDelaunay2D> triangulateImage ( int* dims, 
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
              double* X, 
              double* Y, 
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

int main(int argc, char* argv[])
{
  vtkNew<vtkNamedColors> colors;

  // Verify input arguments
  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0]
              << " Filename(.jpeg/jpg) e.g. Pileated.jpg " << std::endl;
    return EXIT_FAILURE;
  }

  // jacobi poly
  double a = 0.5, b = 0.5, c = 0.5;
  unsigned int N = 55;
  unsigned int m = 6;
  unsigned int M = static_cast<unsigned int>(0.5 * (m + 1) * (m + 2));
  unsigned int nthreads = 6;
  jPoly<double>* Pm = new jPoly<double>(N, m, a, b, c, nthreads); 
  jPoly<double>* Pmx = new jPoly<double>(N, m, a+1, b, c+1, nthreads); 
  jPoly<double>* Pmy = new jPoly<double>(N, m, a, b+1, c+1, nthreads); 
  double* Habc  = (double*) calloc((m+2)*(m+2), sizeof(double));
  double* Ha1bc1  = (double*) calloc((m+2)*(m+2), sizeof(double));
  double* Hab1c1 = (double*) calloc((m+2)*(m+2), sizeof(double));
  double* Dx = (double*) calloc(M*M, sizeof(double));
  double* Dy = (double*) calloc(M*M, sizeof(double));
  double* X = (double*) calloc(N, sizeof(double));
  double* Y = (double*) calloc(N, sizeof(double));
  double* W = (double*) calloc(N, sizeof(double));
  sFactors(m+2, a, b, c, Habc);  
  sFactors(m+2, a+1, b, c+1, Ha1bc1);  
  sFactors(m+2, a, b+1, c+1, Hab1c1);  
  dMat(a, b, c, Habc, Ha1bc1, m, 0, Dx);
  dMat(a, b, c, Habc, Hab1c1, m, 1, Dy); 
  readQuad( "xtri_N55_n9_M91_m12.txt",
            "ytri_N55_n9_M91_m12.txt",
            "wtri_N55_n9_M91_m12.txt",
             N, X, Y, W);
  // compute Jacobi interpolation matrices
  Pm->computeV(X,Y);
  Pmx->computeV(X,Y);
  Pmy->computeV(X,Y);
  // storage for computing sizefield in adapative triangulation
  double* cimg = (double*) calloc(M, sizeof(double));
  double* cdimg = (double*) calloc(M*2, sizeof(double));
  double* dimgr = (double*) calloc(N*2, sizeof(double));
  double* dimgt = (double*) calloc(2*N, sizeof(double));
  // Read the image
  vtkNew<vtkJPEGReader> jpegReader;
  jpegReader->SetFileName(argv[1]);
  jpegReader->Update();
  vtkSmartPointer<vtkImageData> imagedata = jpegReader->GetOutput(); 
  int dims[3]; imagedata->GetDimensions(dims);
  double origin[3]; imagedata->GetOrigin(origin);

  // triangulate uniform subsampled grid on image
  unsigned int nSamp = 10;
  vtkSmartPointer<vtkDelaunay2D> triangulator = triangulateImage(dims, origin, nSamp);
  vtkSmartPointer<vtkPolyData> polytri = triangulator->GetOutput();
  std::cout << polytri->GetNumberOfPoints() << std::endl;
  std::cout << polytri->GetNumberOfCells() << std::endl;
  // create cubic interpolator on image
   vtkSmartPointer<vtkImageInterpolator> interpolator = 
    vtkSmartPointer<vtkImageInterpolator>::New();
  interpolator->Initialize(imagedata);
  interpolator->SetInterpolationModeToCubic();
  std::cout << interpolator->GetNumberOfComponents() << std::endl; 

  // interpolated val storage for each channel
  double* interpc1 = (double*) calloc(N, sizeof(double));
  double* interpc2 = (double*) calloc(N, sizeof(double));
  double* interpc3 = (double*) calloc(N, sizeof(double));
  double interpc[3];
  // cell pt storage
  double v1t[3], v2t[3], v3t[3], tmpwts[3];
  double xqtr[3], xqt[3]; xqtr[2] = 0; xqt[2] = 0;
  vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
  // incidence matrix
  double invIxe[4], gradnorm;
  vtkSmartPointer<vtkPolyData> polytri1 = polytri; 
  unsigned int nRuns = 5;
  for (unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    polytri = polytri1; 
    // gradient density
    double* intgn1 = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
    double* intgn2 = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
    double* intgn3 = (double*) calloc(polytri->GetNumberOfCells(), sizeof(double));
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
        interpolator->Interpolate(xqt, interpc);
        interpc1[j] = interpc[0];
        interpc2[j] = interpc[1];
        interpc3[j] = interpc[2];
      }
      // get tri verts
      polytri->GetCellPoints(i, ptids);
      polytri->GetPoints()->GetPoint(ptids->GetId(0), v1t);
      polytri->GetPoints()->GetPoint(ptids->GetId(1), v2t);
      polytri->GetPoints()->GetPoint(ptids->GetId(2), v3t);
      // incidence matrix inverse
      invIxe[0] = v3t[1] - v1t[1];
      invIxe[1] = -(v2t[1] - v1t[1]);
      invIxe[2] = -(v3t[0] - v1t[0]); 
      invIxe[3] = v2t[0] - v1t[0];
      // coefficients of interpolated data in channel
      Pm->computeCoeffs(interpc1, X, Y, W, cimg);
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
        gradnorm += (dimgt[j] * dimgt[j] + dimgt[j+N] * dimgt[j+N]) * W[j] / 2;
      }
      intgn1[i] = gradnorm;

      Pm->computeCoeffs(interpc2, X, Y, W, cimg);
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
        gradnorm += (dimgt[j] * dimgt[j] + dimgt[j+N] * dimgt[j+N]) * W[j] / 2;
      }
      intgn2[i] = gradnorm;
      
      Pm->computeCoeffs(interpc3, X, Y, W, cimg);
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
        gradnorm += (dimgt[j] * dimgt[j] + dimgt[j+N] * dimgt[j+N]) * W[j] / 2;
      }
      intgn3[i] = gradnorm;
      
    }
    
    double sum = 0, sum1 = 0, sum2 = 0;
    #pragma omp simd reduction(+:sum,sum1,sum2)
    for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
    {
      sum += intgn1[i];
      sum1 += intgn2[i];
      sum2 += intgn3[i];
    }
    for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
    {
      intgn1[i] /= sum;
      intgn2[i] /= sum1;
      intgn3[i] /= sum2;
    } 

    double ave = 0, ave1 = 0, ave2 = 0;
    #pragma omp simd reduction(+:sum,sum1,sum2)
    for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
    {
      ave += intgn1[i];
      ave1 += intgn2[i];
      ave2 += intgn3[i];
    } 
    
    ave /= polytri->GetNumberOfCells();
    ave1 /= polytri->GetNumberOfCells();
    ave2 /= polytri->GetNumberOfCells();

    vtkSmartPointer<vtkIncrementalOctreePointLocator> locator = 
      vtkSmartPointer<vtkIncrementalOctreePointLocator>::New();
    locator->SetDataSet(polytri);
    locator->BuildLocator();
    ptids->Reset();
    double v1r[3] = {0.5, 0, 0};
    double v2r[3] = {0.5, 0.5, 0};
    double v3r[3] = {0, 0.5, 0};
    vtkIdType ptid; int subid;
    std::cout << "num pts in loc: " << locator->GetNumberOfPoints() << std::endl;
    for (int icell = 0; icell < polytri->GetNumberOfCells(); ++icell)
    {
      if (intgn1[icell] > ave)
      {
        polytri->GetCell(icell)->EvaluateLocation(subid, v1r, v1t, tmpwts);
        polytri->GetCell(icell)->EvaluateLocation(subid, v2r, v2t, tmpwts);
        polytri->GetCell(icell)->EvaluateLocation(subid, v3r, v3t, tmpwts);
        locator->InsertUniquePoint(v1t, ptid);  
        locator->InsertUniquePoint(v2t, ptid);  
        locator->InsertUniquePoint(v3t, ptid);
      } 
    } 
    std::cout << "num pts in loc: " << locator->GetNumberOfPoints() << std::endl;
    vtkSmartPointer<vtkPolyData> polydata1 = vtkSmartPointer<vtkPolyData>::New(); 
    polydata1->SetPoints(locator->GetLocatorPoints());
    vtkSmartPointer<vtkDelaunay2D> triangulator1 = vtkSmartPointer<vtkDelaunay2D>::New();
    triangulator1->SetInputData(polydata1);
    triangulator1->Update();
    polytri1 = triangulator1->GetOutput();
    free(intgn1);
    free(intgn2);
    free(intgn3);
  }
  // write mesh
  std::cout << "NEW num points " << polytri1->GetNumberOfPoints() << std::endl;
  vtkNew<vtkXMLPolyDataWriter> writer1;
  writer1->SetFileName("test1.vtp");
  writer1->SetInputData(polytri1);  
  //writer->SetDataModeToBinary();
  writer1->SetDataModeToAscii();
  writer1->Write();
  
  // Visualize
  vtkNew<vtkImageViewer2> imageViewer;
  imageViewer->SetInputData(imagedata);
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  imageViewer->SetupInteractor(renderWindowInteractor);
  imageViewer->Render();
  imageViewer->GetRenderer()->ResetCamera();
  imageViewer->GetRenderer()->SetBackground(
      colors->GetColor3d("DarkSlateGray").GetData());
  imageViewer->GetRenderWindow()->SetWindowName("JPEGReader");
  imageViewer->Render();

  renderWindowInteractor->Start();
 
  // cleanup 
  delete Pm;
  delete Pmx;
  delete Pmy;
  free(Habc);
  free(Ha1bc1);
  free(Hab1c1);
  free(Dx);
  free(Dy);
  free(X); 
  free(Y); 
  free(W);
  free(cdimg);
  free(dimgr);
  free(dimgt); 

  return EXIT_SUCCESS;
}
  

/* Random test codes */

  //for (int icell = 0; icell < polytri->GetNumberOfCells(); ++icell)
  //{

  //  polytri->GetCellPoints(icell, ptids);
  //  polytri->GetPoints()->GetPoint(ptids->GetId(0), v1);
  //  polytri->GetPoints()->GetPoint(ptids->GetId(1), v2);
  //  polytri->GetPoints()->GetPoint(ptids->GetId(2), v3);
  //
  //} 



  // testing coef expansions and derivatives
  //double* F = (double*) calloc(N, sizeof(double));
  //double* cF = (double*) calloc(M, sizeof(double));
  //double* cdxF = (double*) calloc(M, sizeof(double));
  //double* cdyF = (double*) calloc(M, sizeof(double));
  //for (unsigned int i = 0; i < N; ++i) { F[i] = std::pow(X[i] + Y[i], 6); }
  //Pm->computeCoeffs(F, X, Y, W, cF);
  //cblas_dgemv( CblasColMajor, CblasNoTrans,  M, M, 1.0, Dx, M, cF, 1, 0.0, cdxF, 1 );
  //cblas_dgemv( CblasColMajor, CblasNoTrans,  M, M, 1.0, Dy, M, cF, 1, 0.0, cdyF, 1 );
  //for (unsigned int i = 0; i < M; ++i) 
  //{
  //  printf("%.17g\n", cF[i]);
  //}
  //printf("\n");
  //for (unsigned int i = 0; i < M; ++i) 
  //{
  //  printf("%.17g\n", cdxF[i]);
  //}
  //printf("\n");
  //for (unsigned int i = 0; i < M; ++i) 
  //{
  //  printf("%.17g\n", cdyF[i]);
  //}
  //printf("\n");
  //free(F); free(cF);
  //free(cdxF); free(cdyF);

  //std::cout << imagedata->GetNumberOfPoints() << std::endl;
  //std::cout << imagedata->GetNumberOfCells() << std::endl;
  //double x[3];
  //imagedata->GetPoint(24160255, x);
  //std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
  //double origin[3];
  //imagedata->GetOrigin(origin);
  //std::cout << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
  //std::cout << imagedata->GetActualMemorySize() * 1024 << std::endl;
  //int dims[3];
  //imagedata->GetDimensions(dims);
  //std::cout << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;



  //double xqtr[3], xqt[3]; 
  //xqtr[2] = 0; xqt[2] = 0;
  //double interpwts[3];

  //for (unsigned int j = 0; j < N; ++j)
  //{
  //  // get quad point in ref
  //  xqtr[0] = X[j]; xqtr[1] = Y[j]; 
  //  // map quad point into each cell
  //  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  //  {
  //    int subid;
  //    polytri->GetCell(i)->EvaluateLocation(subid, xqtr, xqt, interpwts);
  //    points->InsertNextPoint(xqt); 
  //  } 
  //}



  //vtkDoubleArray* gnint = vtkDoubleArray::New();
  //gnint->SetName("Integral(GradNorm)");
  //gnint->SetNumberOfComponents(3);
  //gnint->SetNumberOfTuples(polytri->GetNumberOfCells());
  //double sums[3];
  //for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  //{
  //  sums[0] = intgn1[i]; sums[1] = intgn2[i]; sums[2] = intgn3[i];
  //  gnint->SetTuple(i, sums);
  //}

  //polytri->GetCellData()->AddArray(gnint); 
  //vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  //points->DeepCopy(polytri->GetPoints());


  //std::cout << "NEW num points " << polytri->GetNumberOfPoints() << std::endl;
  //vtkNew<vtkXMLPolyDataWriter> writer;
  //writer->SetFileName("test.vtp");
  //writer->SetInputData(polytri);  
  ////writer->SetDataModeToBinary();
  //writer->SetDataModeToAscii();
  //writer->Write();
